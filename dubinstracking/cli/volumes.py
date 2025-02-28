from dubinstracking.cli.subapplication import Subapplication
from dubinstracking.util.meshSlicing import polygonFromBody
from argparse import ArgumentParser, FileType
import numpy as np
import os
import subprocess
import json
import logging

try:
    import pyvista as pv
except ModuleNotFoundError:
    class pv:
        def __getattribute__(self, name):
            raise NotImplementedError()
        def __setattribute__(self, name):
            raise NotImplementedError()


def folder(item):
    if os.path.isdir(item):
        return os.path.abspath(item) + '/'
    raise ValueError(f'{item} is not a directory')


def file(item):
    if os.path.isfile(item):
        return os.path.abspath(item)
    raise ValueError(f'{item} is not a file')


class Volumes(Subapplication):
    """
    cli subapplication that creates visibility volumes along a trajectory
    """
    def __init__(self):
        self.logger = logging.getLogger(__name__)
        super().__init__('volumes')

    def modifyParser(self, parser: ArgumentParser):
        parser.add_argument('--csv', dest='csv', type=FileType('r'), required=True)
        parser.add_argument('--dist', dest='dist', type=float, required=False, default=20.0)
        parser.add_argument('--out', dest='out', type=folder, required=True)
        parser.add_argument('--map', dest='map', type=file, required=True)
        parser.add_argument('--radius', dest='radius', type=float, required=True)
        parser.add_argument('--cutoff', dest='cutoff', required=True, type=float)
        parser.add_argument('--alt', dest='alt', required=True, type=float)
        super().modifyParser(parser)

    def run(self, args):
        path = np.genfromtxt(args.csv, delimiter=',')
        points = []
        for i in range(path.shape[1] - 1):
            points.append(path[i, :])
            v = path[i + 1] - path[i]
            d = np.linalg.norm(v)
            n = v / d
            num = np.floor(d / args.dist)
            adj = (d - num * args.dist) / 2
            for j in range(1, int(num)):
                points.append(path[i] + n * (j * args.dist + adj))
        points.append(path[-1, :])
        points = np.array(points)
        with open(args.out + 'points.json', 'w') as f:
            json.dump(points.tolist(), f)

        cmdArgs = [
            'blender', '-b', '--python', 'dubinstracking/cli/helpers/adaptiveSampling.py',
            '--',
            '--out', args.out,
            '--points', args.out + 'points.json',
            '--volumes', args.out + 'volumes.json',
            '--map', args.map,
            '--radius', str(args.radius),
            '--cutoff', str(args.cutoff)
        ]
        try:
            result = subprocess.run(
                cmdArgs,
                check=True
            )
            result.check_returncode()
        except subprocess.CalledProcessError as e:
            self.logger.error(e.stdout.decode('utf-8'))
            self.logger.error(e.stderr.decode('utf-8'))

        with open(args.out + 'volumes.json') as f:
            volumes = json.load(f)

        for volume in volumes:
            if volume['polygon'] != '':
                continue
            file = volume['volume']
            vol = pv.read(file)
            vol = vol.transform(np.array([
                [1, 0, 0, 0],
                [0, 0, -1, 0],
                [0, 1, 0, 0],
                [0, 0, 0, 1]
            ]))
            polygon = polygonFromBody(args.alt, vol)
            if polygon is None:
                self.logger.warning(f'{file} not sliced at {args.alt}')
                continue
            volume['polygon'] = np.array(polygon.exterior.coords).tolist()

        with open(args.out + 'volumes.json', 'w') as f:
            volumes = json.dump(volumes, f)
