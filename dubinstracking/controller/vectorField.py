"""
vector field controller that forces a Dubins vehicle onto a morphing circular orbit
"""
import numpy as np
from scipy.integrate import solve_ivp
from dubinstracking.pathplanning import MorphingOrbit


def angleDiff(a: float, b: float):
    """
    difference between to angles a - b

    Parameters
    ----------
    a: float
        first angle
    b: float
        second angle

    Returns
    -------
    float [-pi, pi]
    """
    c = a - b
    c = (c + np.pi) % (2 * np.pi) - np.pi
    return c


class VectorFieldController:
    """
    vector field controller that forces a Dubins vehicle onto a morphing circular orbit
    """
    def __init__(self, beta, k, maxcurv):
        """
        create the vector field controller

        Parameters
        ----------
        beta: float
            parameter that controller the rate of attraction towards the orbit
        k: float
            proportinal gain that forces the Dubins vehicle heading to the desired heading
        maxcurv: float
            the maximum curvature for the Dubins vehicle 1/r_min
        """
        self.beta = beta
        self.k = k
        self.maxcurv = maxcurv

    def makeU(self, orbit: MorphingOrbit):
        """
        creates a function the is the guidance vector field for the Dubins vehicle
        
        Parameters
        ----------
        orbit: MorphingOrbit
            the current orbit to attract to
        
        Returns
        -------
        (np.ndarray, float) -> np.ndarray\\
        ([x, y], t) -> (u_x, u_y)
        """
        speed = orbit.speed
        def u(x: np.ndarray, t: float) -> np.ndarray:
            vel = orbit.g_dot(t)
            g_v = np.linalg.norm(vel)
            g = orbit.g(t)
            theta = np.arctan2(x[1] - g[1], x[0] - g[0])
            radius_dot = orbit.radius_dot(theta, t)
            radius = orbit.radius(theta, t)
            r = np.sqrt((x[0] - g[0]) ** 2 + (x[1]- g[1]) ** 2)
            e_r = np.array([(x[0] - g[0]), (x[1] - g[1])]) / r
            e_theta = np.array([-x[1] + g[1], x[0] - g[0]]) / r

            G = (speed - g_v - np.abs(radius_dot)) * 2 / np.pi * np.arctan(self.beta * (r - radius))
            H = np.sqrt(speed ** 2 - (-G + np.dot(vel, e_r) + radius_dot) ** 2)

            return (-G + np.dot(vel, e_r) + radius_dot) * e_r + H * e_theta
        return u

    def makeSys(self, orbit: MorphingOrbit):
        """
        creates a function that returns the derivative of the current state for the Dubins vehicle, the function for scipy.integrate.solve_ivp.

        Parameters
        ----------
        orbit: MorphingOrbit
            the orbit to attract the Dubins Vehilce to
        
        Returns
        -------
        (np.ndarray, float) -> np.ndarray\\
        ([x, y, psi], t) -> (x_dot, y_dot, psi_dot)
        """
        speed = orbit.speed
        u = self.makeU(orbit)
        def psi_dot(x, t):
            vel = orbit.g_dot(t)
            g_v = np.linalg.norm(vel)
            g = np.array([vel[0]  * t, vel[1] * t])
            theta = np.arctan2((x[1]- g[1]), (x[0]- g[0]))
            radius = orbit.radius(theta, t)
            radius_dot = orbit.radius_dot(theta, t)
            r = np.sqrt((x[0] - g[0]) ** 2 + (x[1]- g[1]) ** 2)
            e_theta = np.array([-x[1] + g[1], x[0] - g[0]]) / r
            u_theta = np.dot(e_theta, u(x, t))
            theta_dot = u_theta / r - np.dot(vel, e_theta) / r
            dot_phi_p = -4 / np.pi ** 2 * self.beta * np.arctan(self.beta * (r - radius)) * ((speed - g_v - np.abs(radius_dot)) ** 2) / (1 + self.beta ** 2 * (r - radius) ** 2)
            return dot_phi_p / u_theta + r * theta_dot ** 2 / (u_theta)
        
        def invert(x, t):
            vel = orbit.g_dot(t)
            g = orbit.g(t)
            theta = np.arctan2(x[1]- g[1],x[0] - g[0] )
            g = np.array([vel[0]  * t, vel[1] * t])
            radius = orbit.radius(theta, t)
            radius_dot = orbit.radius_dot(theta, t)
            r = np.sqrt((x[0] - g[0]) ** 2 + (x[1]- g[1]) ** 2)
            e_theta = np.array([-x[1] + g[1], x[0] - g[0]]) / r
            e_r = np.array([x[0] - g[1], x[1] - g[1]]) / r
            u_theta = np.dot(e_theta, u(x, t))
            
            psi_d = np.arctan2(u(x, t)[1], u(x, t)[0])
            psi = x[2]
            psi_diff = (psi - psi_d + np.pi) % (2 * np.pi) - np.pi
            a = (np.dot(g, e_r) + radius_dot) * (1 - np.cos(psi_diff))/(psi_diff)
            b = u_theta * np.sin(psi_diff) / psi_diff
            c = (self.beta * 2 / np.pi * np.arctan(self.beta * (r - radius)) / (1 + self.beta ** 2 * (r - radius) ** 2))
            return (a + b) * c

        def sys(t, state):

            xy = u(state[:2], t)
            psi = state[2]

            
            psiff = psi_dot(state[:2], t)
            psi_d = np.arctan2(xy[1], xy[0])
            dpsi = self.k * angleDiff(psi_d, psi) + psiff + invert(state, t)
            dpsi = np.clip(dpsi, -speed * self.maxcurv, speed * self.maxcurv)
            return [
                speed * np.cos(psi),
                speed * np.sin(psi),
                dpsi
            ]
        return sys

    def solve(self, orbit: MorphingOrbit, t_start: float, t_final: float, x_0: 'list', t_step: float):
        """
        solve for the trajectory for a Dubins vehicle attracting to an orbit

        Parameters
        ----------
        orbit: MorphingOrbit
            the orbit to attract to
        t_start: float
            the time at the start of the trajectory
        t_end: float
            the time at the end of the trajectory
        x_0: list[float]
            the initial configuration of the Dubins vehicle
        t_step: float
            the maximum time step along the trajectory
        
        Returns
        -------
        np.ndarray, np.ndarray, np.ndarray, np.ndarray\\
        time [n], Dubins vehicle state (x, y, psi, x_dot, y_dot, psi_dot) [n, 6], 
        target position (g_x, g_y) [n: 2], guidance vector field (u_x, u_y) [n, 2]
        """
        sys = self.makeSys(orbit)
        u = self.makeU(orbit)
        r = solve_ivp(sys, [t_start, t_final], x_0, max_step=t_step)
        vel = [sys(t, x) for t, x in zip(r.t, r.y.T)]
        return r.t, np.column_stack([r.y.T, vel]), np.array([orbit.g(t) for t in r.t]), np.array([u(x, t) for x, t in zip(r.y.T, r.t)])

    def control(self, orbit: MorphingOrbit, state: np.ndarray, t: float):
        """
        get the current control

        Parameters
        ----------
        orbit: MorphingOrbit
            the orbit to attract to
        state: np.ndarray
            the Dubins vehicle state [x, y, psi]
        t: float
            current time
        
        Returns
        -------
        float: u_psi the control
        """
        sys = self.makeSys(orbit)
        ds = sys(t, state)
        return ds[2]
    
    def controlU(self, orbit: MorphingOrbit, state: np.ndarray, t: float):
        """
        get the current guidance

        Parameters
        ----------
        orbit: MorphingOrbit
            the orbit to attract to
        state: np.ndarray
            the xy position [x, y]
        t: float
            current time
        
        Returns
        -------
        np.ndarray: (u_x, u_y) the guidance vector
        """
        u = self.makeU(orbit)
        ds = u(state, t)
        return ds