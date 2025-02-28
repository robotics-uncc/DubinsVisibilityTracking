from targettracking.cli.subapplication import Subapplication
from targettracking.pathplanning import OrbitFactory, Orbit
from argparse import ArgumentParser, FileType
import numpy as np
import os
from shapely.geometry import Polygon, Point
from shapely.affinity import affine_transform
import json
import logging
from tqdm import tqdm
from scipy.optimize import basinhopping


CUTOFF = 1e-2


def folder(item):
    """
    is an item a folder

    Parameters
    ----------
    item: str
        to test if folder
    
    Returns
    -------
    str: the absolute path to the folder

    Raises
    ------
    ValueError: when not a folder
    """
    if os.path.isdir(item):
        return os.path.abspath(item) + '/'
    raise ValueError(f'{item} is not a directory')


def getCurvature(R, g_dot, g_ddot, speed, ccw):
    """
    get a function that returns the curvature of an orbit

    Parameters
    ----------
    R: float
        radius
    g_dot: np.ndarray
        translational rate of the orbit
    g_ddot: np.ndarray
        translational acceleration of the orbit
    speed: float
        speed of the agent on the orbit
    ccw: bool
        direction of the orbit true = ccw
    """
    def curvature(theta, *args):
        theta = (theta + 2 * np.pi) % (2 * np.pi)
        e_r = np.array([np.cos(theta), np.sin(theta)]).flatten()
        e_theta = np.array([-np.sin(theta), np.cos(theta)]).flatten()

        a = R ** 2
        b = R * np.dot(g_dot, e_theta)
        c = np.dot(g_dot, g_dot) - speed ** 2

        if b ** 2 < a * c:
            raise Exception('Target Moving Too Fast')

        theta_dot = (-b + np.sqrt(b ** 2 - a * c)) / (a) if ccw else (-b - np.sqrt(b ** 2 - a * c)) / (a)

        x_dot = g_dot + R * theta_dot * e_theta

        f1 = b ** 2 - a * c

        b_dot = (
            R * np.dot(g_ddot, e_theta) - R * theta_dot * np.dot(g_dot, e_r)
        )

        c_dot = (
            2 * np.dot(g_dot, g_ddot)
        )
        f1_dot = 2 * b * b_dot - a * c_dot
        f2_dot = -b_dot

        theta_ddot = (f2_dot + f1_dot / (2 * np.sqrt(f1))) / a if ccw else (f2_dot - f1_dot / (2 * np.sqrt(f1))) / a

        x_ddot = (
            g_ddot + (- R * theta_dot ** 2) * e_r
            + (R * theta_ddot) * e_theta
        )

        return np.cross(x_dot, x_ddot) / speed ** 3
    return curvature


def outside(R: float, polygon: Polygon, n: int = 200):
    """
    is an orbit with center (0, 0) inside or outside of a polygon

    Parameters
    ----------
    R: float
        radius of the orbit
    polygon: Polygon
        polygon to check if the orbit is inside of
    n: int
        number of points to check on the orbit
    
    Returns
    -------
    bool
    """
    theta = np.linspace(0, 2 * np.pi, n)
    points = R * np.array([np.cos(theta), np.sin(theta)])
    for i in range(points.shape[1]):
        if not polygon.contains(Point(*points[:, i])):
            return True
    return False


class Circles(Subapplication):
    """
    cli subapplication that creates visibility orbits inside a set of visibility volumes
    """
    def __init__(self):
        self.logger = logging.getLogger(__name__)
        super().__init__('circles')

    def modifyParser(self, parser: ArgumentParser):
        parser.add_argument('--volumes', dest='volumes', type=FileType('r'), required=True)
        parser.add_argument('--out', dest='out', type=folder, required=True)
        parser.add_argument('--tspeed', dest='tspeed', type=float, required=False, default=5.0)
        parser.add_argument('--uavspeed', dest='uavspeed', type=float, required=False, default=20.0)
        parser.add_argument('--uavturn', dest='uavturn', type=float, required=False, default=10.0)
        super().modifyParser(parser)

    def run(self, args):
        logger = logging.getLogger(__name__)
        if args.tspeed > args.uavspeed:
            raise Exception(f'UAV is not fast enough to keep up with target, condition {args.tspeed} < {args.uavspeed}.')
        volumes = json.load(args.volumes)
        polygons = [Polygon(volume['polygon']) for volume in volumes]
        points = np.array([volume['point'] for volume in volumes])
        radii = []
        r_dot_cutoff = (args.uavspeed - args.tspeed)
        t = 0
        for i, polygon in tqdm(enumerate(polygons), 'orbit'):
            newPoly = affine_transform(polygon, [1, 0, 0, 1, -points[i, 0], -points[i, 1]])
            r = args.uavturn * .1
            rmin = 0
            rmax = np.inf
            for k in range(200):
                r = (rmin + rmax) / 2 if not np.isinf(rmax) else r * 2
                inside = (not outside(r, newPoly))
                if inside:
                    rmin = r
                else:
                    rmax = r
                if abs(rmax - rmin) < CUTOFF:
                    break

            if i < len(polygons) - 1:
                direction = points[i + 1, :] - points[i, :]
            else:
                direction = points[i, :] - points[i - 1, :]
            distance = np.linalg.norm(direction)
            direction /= distance
            g_dot_0 = direction[:2] * args.tspeed
            time = distance / args.tspeed

            radii.append((r, t, g_dot_0))
            t += time

            curv = getCurvature(r, g_dot_0, np.zeros_like(g_dot_0), args.uavspeed, True)
            minCurv = basinhopping(curv, np.pi)
            maxCurv = basinhopping(lambda x: -curv(x), np.pi)
            if abs(1 / minCurv.fun) < args.uavturn:
                raise Exception(f'Orbit with radius {r} violates curvature constraint |k| < {1/args.uavturn} with a value of {minCurv.fun} at {minCurv.x}')

            if abs(1 / maxCurv.fun) < args.uavturn:
                raise Exception(f'Orbit with radius {r} violates curvature constraint |k| < {1/args.uavturn} with a value of {maxCurv.fun} at {maxCurv.x}')

        for i in range(len(radii) - 1):
            delta = (radii[i + 1][1] - radii[i][1])
            r_dot = (radii[i + 1][0] - radii[i][0]) / delta
            if r_dot > r_dot_cutoff:
                logger.warning(f'R is changing too fast reducing radius at {i + 1}')
                radii[i + 1] = (radii[i][0] + r_dot_cutoff * delta, radii[i + 1][1], radii[i + 1][2])

        for i in range(len(radii) - 2,  -1, -1):
            delta = (radii[i + 1][1] - radii[i][1])
            r_dot = (radii[i + 1][0] - radii[i][0]) / delta
            if r_dot < -r_dot_cutoff:
                logger.warning(f'R is changing too fast reducing radius at {i}')
                radii[i] = (radii[i + 1][0] + r_dot_cutoff * delta, radii[i][1], radii[i][2])

        thetas = np.linspace(0, 2 * np.pi, 16)
        orbits: 'list[Orbit]' = []
        for i, radius in enumerate(radii):
            q = np.array([thetas, radius[0] * np.ones_like(thetas)])
            orbits.append(OrbitFactory.make(q.T, radius[2], points[i, :2], args.uavspeed))

        with open(args.out + 'circles.json', 'w') as f:
            json.dump([orbit.to_dict() for orbit in orbits], f)

        logging.info('Found valid orbits for view polygons')
