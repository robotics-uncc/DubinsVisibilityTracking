"""
    an orbit that translates and grows/shrinks
"""
from typing import Callable
from dubinstracking.util import fit, fit2deriv, fit2deriv2, fit2func
import numpy as np
from enum import IntEnum


class Value(IntEnum):
    """
    parameter to get a value from the orbit
    """
    UNKNOWN = 0
    XY = 1
    THETA_D = 2
    XY_D = 3
    THETA_DD = 4
    XY_DD = 5
    CURVATURE = 6
    RADIUS = 7
    RADIUS_DOT = 8


class OrbitFactory:
    """
        creates orbits
    """
    @staticmethod
    def make(points, g_dot_0, g_0, speed, ccw=True):
        """
        make an orbit from a set of points

        Parameters
        ----------
        points: np.ndarray [n, 2] (theta, r)
            orbit points to fit a spline to
        g_dot_0: np.ndarray
            translating vector
        g_0: np.ndarray
            starting poitn for the orbit
        ccw:
            direction of rotation about the center of the orbit
        """
        coefficients = fit(points, type=2)
        radius = fit2func(points, coefficients)
        radius_prime = fit2deriv(points, coefficients)
        radius_pprime = fit2deriv2(points, coefficients)

        return Orbit(
            points,
            coefficients,
            g_0,
            g_dot_0,
            radius,
            radius_prime,
            radius_pprime,
            speed,
            ccw=ccw
        )

    @staticmethod
    def fromJson(json):
        """
        create an orbit from a json file
        """
        return OrbitFactory.make(np.array(json['points']), np.array(json['g_dot_0']), np.array(json['g_0']), json['speed'])


class Orbit:
    """
    An orbit made of a spline and a translating rate.
    """
    def __init__(
        self,
        points: np.ndarray,
        coefficients: np.ndarray,
        g_0: np.ndarray,
        g_dot_0: np.ndarray,
        radius: Callable[[float], float],
        radius_prime: Callable[[float], float],
        radius_pprime: Callable[[float], float],
        speed: float,
        ccw: bool = True
    ):
        """
        Parameters
        ----------
        points: np.ndarray [n, 2] (theta, r)
            orbit points to fit a spline to
        coefficients: np.ndarray [4n]
            a vector of cubic spline coefficients
        g_0: np.ndarray
            starting point for the orbit
        g_dot_0: np.ndarray
            translating velocity for the orbit
        radius: (theta) -> r 
            gets the radius at a specific theta
        radius_prime: (theta) -> dr/dtheta
            gets the radius's partial derivative with respect to theta 
        radius_pprim: (theta) -> d^2r/dtheta^2
            gets the radius's second derivative with respect to theta
        speed: float
            speed of the Dubins vehicle
        ccw: bool
            true for counter clockwise rotation

        """
        self.points = points
        self.coefficients = coefficients
        self.g_0 = g_0
        self.g_dot_0 = g_dot_0
        self.g_ddot_0 = np.zeros_like(self.g_dot_0)
        self.radius = radius
        self.radius_prime = radius_prime
        self.radius_pprime = radius_pprime
        self.speed = speed
        self.ccw = ccw

    def find_value(self, theta, value: Value):
        """
        finds a value of the orbit at a specific theta

        Parameters
        ----------
        theta: float
            angle to find the value at
        value: Value
            enum describing what value to find
        """
        if value == Value.UNKNOWN:
            raise Exception('Unkown Value')

        e_r = np.array([np.cos(theta), np.sin(theta)])
        e_theta = np.array([-np.sin(theta), np.cos(theta)])
        R = self.radius(theta)
        g_0 = self.g_0
        g_dot = self.g_dot_0
        g_ddot = self.g_ddot_0

        if value == Value.RADIUS:
            return R

        if value == Value.XY:
            return R * e_r + g_0

        lambda_dot = 0
        lambda_prime = self.radius_prime(theta)
        a = lambda_prime ** 2 + R ** 2
        b = lambda_prime * np.dot(g_dot, e_r) + R * np.dot(g_dot, e_theta) + lambda_prime * lambda_dot
        c = np.dot(g_dot, g_dot) - self.speed ** 2 + 2 * lambda_dot * np.dot(g_dot, e_r) + lambda_dot ** 2

        if b ** 2 < a * c:
            raise Exception('Radius changing too fast')

        theta_dot = (-b + np.sqrt(b ** 2 - a * c)) / (a) if self.ccw else (-b - np.sqrt(b ** 2 - a * c)) / (a)

        if value == Value.THETA_D:
            return theta_dot

        x_dot = g_dot + (theta_dot * lambda_prime + lambda_dot) * e_r + R * theta_dot * e_theta

        if value == Value.XY_D:
            return x_dot

        lambda_pprime = self.radius_pprime(theta)
        lambda_dot_prime = 0

        f1 = b ** 2 - a * c

        b_dot = (
            (lambda_pprime * theta_dot + lambda_dot_prime) * np.dot(g_dot, e_r) + lambda_prime * np.dot(g_ddot, e_r)
            + lambda_prime * theta_dot * np.dot(g_ddot, e_theta) + (lambda_prime * theta_dot + lambda_dot) * np.dot(g_dot, e_theta)
            + R * np.dot(g_ddot, e_theta) - R * theta_dot * np.dot(g_dot, e_r)
            + (lambda_pprime * theta_dot + lambda_dot_prime) * lambda_dot + lambda_prime * lambda_dot_prime * theta_dot
        )

        a_dot = (
            2 * lambda_prime * (lambda_pprime * theta_dot + lambda_dot_prime) +
            2 * R * (lambda_prime * theta_dot + lambda_dot)
        )

        c_dot = (
            2 * np.dot(g_dot, g_ddot) + 2 * lambda_dot_prime * theta_dot * np.dot(g_dot, e_r) + 2 * lambda_dot * np.dot(g_ddot, e_r) +
            2 * lambda_dot * theta_dot * np.dot(g_dot, e_theta) + 2 * lambda_dot * lambda_dot_prime * theta_dot
        )
        f1_dot = 2 * b * b_dot - a * c_dot - a_dot * c
        f2 = -b
        f2_dot = -b_dot

        theta_ddot = (a * (f2_dot + f1_dot / (2 * np.sqrt(f1))) - a_dot * (f2 + np.sqrt(f1))) / \
            (a ** 2) if self.ccw else (a * (f2_dot - f1_dot / (2 * np.sqrt(f1))) - a_dot * (f2 - np.sqrt(f1))) / (a ** 2)

        if value == Value.THETA_DD:
            return theta_ddot

        x_ddot = (
            g_ddot + (theta_ddot * lambda_prime + lambda_pprime * theta_dot ** 2 + 2 * lambda_dot_prime * theta_dot - R * theta_dot ** 2) * e_r
            + (lambda_prime * theta_dot + 2 * lambda_dot * theta_dot + lambda_prime * theta_dot ** 2 + R * theta_ddot) * e_theta
        )

        if value == Value.XY_DD:
            return x_ddot

        curvature = np.cross(x_dot, x_ddot) / self.speed ** 3
        if value == Value.CURVATURE:
            return curvature

        raise Exception('Unknown Value')

    def to_dict(self):
        """
        chagne the object into a dict for export to json
        """
        return {
            'points': self.points.tolist(),
            'g_0': self.g_0.tolist(),
            'g_dot_0': self.g_dot_0.tolist(),
            'speed': self.speed,
            'ccw': self.ccw
        }

    def theta_dot(self, theta):
        """
        gets the polar angle rate about the orbit

        Parameters
        ----------
        theta: float
            angle about the point g_0
        
        Returns
        float
        """
        return self.find_value(theta, Value.THETA_D)

    def theta_ddot(self, theta):
        """
        gets the polar angle acceleration about the orbit

        Parameters
        ----------
        theta: float
            angle about the point g_0
        
        Returns
        float
        """
        return self.find_value(theta, Value.THETA_DD)

    def curvature(self, theta):
        """
        gets the curvature of the orbit at theta

        Parameters
        ----------
        theta: float
            angle about the point g_0
        
        Returns
        float
        """
        return self.find_value(theta, Value.CURVATURE)

    def xy_dot(self, theta):
        """
        gets the xy velocity on the orbit at theta

        Parameters
        ----------
        theta: float
            angle about the point g_0
        
        Returns
        np.ndarray
        """
        return self.find_value(theta, Value.XY_D)

    def xy_ddot(self, theta):
        """
        gets the xy acceleration on the orbit at theta

        Parameters
        ----------
        theta: float
            angle about the point g_0
        
        Returns
        np.ndarray
        """
        return self.find_value(theta, Value.XY_DD)

    def xy(self, theta):
        """
        gets the xy position on the orbit at theta

        Parameters
        ----------
        theta: float
            angle about the point g_0
        
        Returns
        np.ndarray
        """
        return self.find_value(theta, Value.XY)


class MorphingOrbit:
    """
    the orbit creates by morphing from one orbit to another
    """
    def __init__(self, orbit1: Orbit, orbit2: Orbit, t1: float, t2: float):
        """
        Parameters
        ----------
        orbit1: Orbit
            orbit at t1
        orbit2: Orbit:
            orbit at t2
        t1: float
            start time starts at orbit1
        t2: float
            end time arrives at orbit2
        """
        self.orbit1 = orbit1
        self.orbit2 = orbit2
        self.t1 = t1
        self.t2 = t2
        self.delta = t2 - t1
        self.g_0 = self.orbit1.g_0
        self.g_dot_0 = self.orbit1.g_dot_0
        self.g_ddot_0 = self.orbit1.g_ddot_0
        self.speed = self.orbit1.speed
        self.ccw = self.orbit1.ccw

    def g(self, t):
        """
        gets the center of the orbit

        Parameters
        ----------
        t: float
            current time t1 <= t <= t2
        """
        return self.g_0 + (t - self.t1) * self.g_dot_0

    def g_dot(self, t):
        """
        gets the translational velocity of the orbit

        Parameters
        ----------
        t: float
            current time t1 <= t <= t2
        """
        return self.g_dot_0

    def g_ddot(self, t):
        """
        gets the translational acceleration of the orbit

        Parameters
        ----------
        t: float
            current time t1 <= t <= t2
        """
        return self.g_ddot_0

    def find_value(self, theta, t, value: Value):
        if value == Value.UNKNOWN:
            raise Exception('Unkown Value')

        e_r = np.array([np.cos(theta), np.sin(theta)])
        e_theta = np.array([-np.sin(theta), np.cos(theta)])
        R1 = self.orbit1.radius(theta)
        R2 = self.orbit2.radius(theta)
        R = R1 * (self.t2 - t) / self.delta + R2 * (t - self.t1) / self.delta
        g_t = self.g(t)
        g_dot = self.g_dot(t)
        g_ddot = self.g_ddot(t)

        if value == Value.RADIUS:
            return R

        if value == Value.XY:
            return R * e_r + g_t

        lambda_dot = 1 / self.delta * (R2 - R1)

        if value == Value.RADIUS_DOT:
            return lambda_dot
    
        lambda_prime = self.orbit1.radius_prime(theta) * (self.t2 - t) / self.delta + self.orbit2.radius_prime(theta) * (t - self.t1) / self.delta
        a = lambda_prime ** 2 + R ** 2
        b = lambda_prime * np.dot(g_dot, e_r) + R * np.dot(g_dot, e_theta) + lambda_prime * lambda_dot
        c = np.dot(g_dot, g_dot) - self.speed ** 2 + 2 * lambda_dot * np.dot(g_dot, e_r) + lambda_dot ** 2

        if b ** 2 < a * c:
            raise Exception('Radius changing too fast')

        theta_dot = (-b + np.sqrt(b ** 2 - a * c)) / (a) if self.ccw else (-b - np.sqrt(b ** 2 - a * c)) / (a)

        if value == Value.THETA_D:
            return theta_dot

        x_dot = g_dot + (theta_dot * lambda_prime + lambda_dot) * e_r + R * theta_dot * e_theta

        if value == Value.XY_D:
            return x_dot

        lambda_pprime = self.orbit1.radius_pprime(theta) * (self.t2 - t) / self.delta + self.orbit2.radius_pprime(theta) * (t - self.t1) / self.delta
        lambda_dot_prime = 1 / self.delta * (self.orbit2.radius_prime(theta) - self.orbit1.radius_prime(theta))

        f1 = b ** 2 - a * c

        b_dot = (
            (lambda_pprime * theta_dot + lambda_dot_prime) * np.dot(g_dot, e_r) + lambda_prime * np.dot(g_ddot, e_r)
            + lambda_prime * theta_dot * np.dot(g_ddot, e_theta) + (lambda_prime * theta_dot + lambda_dot) * np.dot(g_dot, e_theta)
            + R * np.dot(g_ddot, e_theta) - R * theta_dot * np.dot(g_dot, e_r)
            + (lambda_pprime * theta_dot + lambda_dot_prime) * lambda_dot + lambda_prime * lambda_dot_prime * theta_dot
        )

        a_dot = (
            2 * lambda_prime * (lambda_pprime * theta_dot + lambda_dot_prime) +
            2 * R * (lambda_prime * theta_dot + lambda_dot)
        )

        c_dot = (
            2 * np.dot(g_dot, g_ddot) + 2 * lambda_dot_prime * theta_dot * np.dot(g_dot, e_r) + 2 * lambda_dot * np.dot(g_ddot, e_r) +
            2 * lambda_dot * theta_dot * np.dot(g_dot, e_theta) + 2 * lambda_dot * lambda_dot_prime * theta_dot
        )
        f1_dot = 2 * b * b_dot - a * c_dot - a_dot * c
        f2 = -b
        f2_dot = -b_dot

        theta_ddot = (a * (f2_dot + f1_dot / (2 * np.sqrt(f1))) - a_dot * (f2 + np.sqrt(f1))) / \
            (a ** 2) if self.ccw else (a * (f2_dot - f1_dot / (2 * np.sqrt(f1))) - a_dot * (f2 - np.sqrt(f1))) / (a ** 2)

        if value == Value.THETA_DD:
            return theta_ddot

        x_ddot = (
            g_ddot + (theta_ddot * lambda_prime + lambda_pprime * theta_dot ** 2 + 2 * lambda_dot_prime * theta_dot - R * theta_dot ** 2) * e_r
            + (lambda_prime * theta_dot + 2 * lambda_dot * theta_dot + lambda_prime * theta_dot ** 2 + R * theta_ddot) * e_theta
        )

        if value == Value.XY_DD:
            return x_ddot

        curvature = np.cross(x_dot, x_ddot) / self.speed ** 3
        if value == Value.CURVATURE:
            return curvature

        raise Exception('Unknown Value')

    def theta_dot(self, theta, t):
        """
        gets the angular rate about g at theta, t
    
        Parameters
        ----------
        theta: float
            angle about g
        t: float
            time t1 <= t <= t2
        """
        return self.find_value(theta, t, Value.THETA_D)

    def theta_ddot(self, theta, t):
        """
        gets the angular acceleration about g at theta, t
    
        Parameters
        ----------
        theta: float
            angle about g
        t: float
            time t1 <= t <= t2
        """
        return self.find_value(theta, t, Value.THETA_DD)

    def curvature(self, theta, t):
        """
        gets the curvature of the orbit at theta, t
    
        Parameters
        ----------
        theta: float
            angle about g
        t: float
            time t1 <= t <= t2
        """
        return self.find_value(theta, t, Value.CURVATURE)

    def xy_dot(self, theta, t):
        """
        gets the velocity at theta, t
    
        Parameters
        ----------
        theta: float
            angle about g
        t: float
            time t1 <= t <= t2
        """
        return self.find_value(theta, t, Value.XY_D)

    def xy_ddot(self, theta, t):
        """
        gets the acceleration at theta, t
    
        Parameters
        ----------
        theta: float
            angle about g
        t: float
            time t1 <= t <= t2
        """
        return self.find_value(theta, t, Value.XY_DD)

    def xy(self, theta, t):
        """
        gets the position at theta, t
    
        Parameters
        ----------
        theta: float
            angle about g
        t: float
            time t1 <= t <= t2
        """
        return self.find_value(theta, t, Value.XY)

    def radius(self, theta, t):
        """
        gets the radius at t
    
        Parameters
        ----------
        theta: None
            ignored
        t: float
            time t1 <= t <= t2
        """
        return self.find_value(theta, t, Value.RADIUS)
    
    def radius_dot(self, theta, t):
        """
        gets the radial rate at t
    
        Parameters
        ----------
        theta: None
            ignored
        t: float
            time t1 <= t <= t2
        """
        return self.find_value(theta, t, Value.RADIUS_DOT)
