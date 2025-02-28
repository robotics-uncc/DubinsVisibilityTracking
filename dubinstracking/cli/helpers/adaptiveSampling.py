"""
run in blender
"""
import os
import bpy
import json
from argparse import ArgumentParser
import subprocess
import sys


def clear():
    # clear workspace
    for item in bpy.data.objects:
        bpy.data.objects.remove(item)


def decimate(file, faces):
    # decimate view volumes and change file type
    clear()
    bpy.ops.wm.obj_import(filepath=file)
    viewVolume = bpy.context.selected_objects[0]

    ratio = faces / len(viewVolume.data.polygons)
    decimate = viewVolume.modifiers.new(name='decimate', type='DECIMATE')
    decimate.ratio = ratio
    decimate.use_collapse_triangulate = True
    viewVolume.select_set(True)
    with bpy.context.temp_override(object=viewVolume):
        bpy.ops.object.modifier_apply(modifier=decimate.name)

    bpy.ops.export_mesh.ply(filepath=file.replace('.obj', '.ply'))
    if os.path.exists(file.replace('.obj', '.ply')):
        os.remove(file)


def viewRegionDifference(a: str, b: str, aPoints: 'list', bPoints: 'list', remesh: float):
    clear()

    # import view volume
    bpy.ops.wm.obj_import(filepath=a)
    aVolume = bpy.context.selected_objects[0]
    bpy.ops.transform.translate(value=[-1 * p for p in aPoints])
    bpy.ops.wm.obj_import(filepath=b)
    bVolume = bpy.context.selected_objects[0]
    bpy.ops.transform.translate(value=[-1 * p for p in bPoints])
    # boolean difference
    # select imported object
    bool = aVolume.modifiers.new(name='difference', type='BOOLEAN')
    bool.object = bVolume
    bool.operation = 'DIFFERENCE'
    bool.use_self = True
    aVolume.select_set(True)
    with bpy.context.temp_override(object=aVolume):
        bpy.ops.object.modifier_apply(modifier=bool.name)

    # remesh intersection
    re = aVolume.modifiers.new(name='remesh', type='REMESH')
    re.voxel_size = remesh
    aVolume.select_set(True)
    with bpy.context.temp_override(object=aVolume):
        bpy.ops.object.modifier_apply(modifier=re.name)
    aVolume.select_set(False)
    bVolume.select_set(True)
    bpy.ops.object.delete()
    aVolume.select_set(True)
    # save
    bpy.context.view_layer.objects.active = aVolume
    bpy.ops.rigidbody.object_add()
    bpy.ops.rigidbody.mass_calculate(material='Custom', density=1)
    # mass of density 1 is volume
    return aVolume.rigid_body.mass


def metric(aFile: str, bFile: str, aPoints: 'list', bPoints: 'list', remesh: float):
    # compute difference between visibility volumes
    a_volume = viewRegionDifference(aFile, bFile, aPoints, bPoints, remesh)
    b_volume = viewRegionDifference(bFile, aFile, aPoints, bPoints, remesh)
    return a_volume + b_volume


def main(args):
    # run view volume generation program
    if not os.path.exists(args.out):
        os.mkdir(args.out)
    # load in view volume data
    if args.volumes is not None and os.path.exists(args.volumes):
        with open(args.volumes) as f:
            volumes = json.load(f)
    elif args.points is not None and os.path.exists(args.points):
        with open(args.points) as f:
            points = json.load(f)
        volumes = [
            {
                'volume': None,
                'point': points[i],
                'polygon': ''
            }
            for i in range(len(points))
        ]
    else:
        raise Exception('--volumes or --points must exist')
    # create volumes
    volumes = makeVolumes(volumes, args.out, args.map, args.radius)
    # make queue of adjacent volumes
    queue = [(volumes[i], volumes[i + 1]) for i in range(len(volumes) - 1)]
    order = []
    while len(queue) > 0:
        # place volumes in the correct output sequence
        a, b = queue.pop(0)
        if a not in order and b not in order:
            order.append(a)
            order.append(b)
        elif a in order and b in order:
            pass
        elif a in order:
            i = order.index(a)
            order.insert(i + 1, b)
        elif b in order:
            i = order.index(b)
            order.insert(i, a)

        # make visibility volumes if they don't exist
        if a['volume'] is None or b['volume'] is None:
            volumes = makeVolumes(volumes, args.out, args.map, args.radius)
        # recursively call adaptive sampling
        queue += adaptive(a, b, volumes, args.cutoff)

    # enusre all visibility volumes are created
    volumes = makeVolumes(volumes, args.out, args.map, args.radius)

    # calculate the metric along the veiw volumes
    metricOut = [metric(order[k]['volume'], order[k + 1]['volume'], order[k]['point'], order[k + 1]['point'], 2.5) for k in range(len(volumes) - 1)]
    with open(args.out + 'metric.json', 'w') as f:
        json.dump(metricOut, f)

    # for volume in volumes:
    #     if volume['volume'].endswith('.obj'):
    #         decimate(volume['volume'], 500)
    #         volume['volume'] = volume['volume'].replace('.obj', '.ply')

    # save visibility volume data
    with open(args.out + 'volumes.json', 'w') as f:
        json.dump(order, f)
    points = [volume['point'] for volume in volumes]
    with open(args.out + 'points.json', 'w') as f:
        json.dump(points, f)


def makeVolumes(volumes: 'list', out: str, worldMap: str, radius: float):
    # read input file templates
    with open('data/template/worldTemplate.txt') as f:
        worldTemplate = f.read()
    with open('data/template/volumeTemplate.txt') as f:
        volumeTemplate = f.read()
    outfile = out + f'out_{os.getpid()}.yaml'
    # write input file 
    k = 0
    with open(outfile, 'w') as f:
        f.write(worldTemplate.format('world', os.path.abspath(worldMap)))
        for i, volume in enumerate(volumes):
            if volume['volume'] is not None or volume['volume'] == '':
                continue
            fname = f'{out}v_{i:03d}.obj'
            point = volume['point']
            f.write(volumeTemplate.format(i, point[0], -point[1], point[2], radius, os.path.abspath(fname)))
            volume['volume'] = fname
            k += 1
    # call view volume generator as subprocess
    if k > 0:
        result = subprocess.run(
            ['./ogl_depthrenderer', '-c', os.path.abspath(outfile)],
            cwd=os.path.abspath('subs/OpenGLDepthRenderer/build/bin/'),
            capture_output=True
        )
        result.check_returncode()
    os.remove(outfile)

    return volumes


def adaptive(aVolume: dict, bVolume: dict, volumes: 'list', cutoff: float):
    aFile = aVolume['volume']
    bFile = bVolume['volume']
    # evaluate metric on two volumes
    value = metric(aFile, bFile, aVolume['point'], bVolume['point'], 2.5)
    # if metric is met return
    if value < cutoff:
        print(aFile, bFile, f'{value} < {cutoff}')
        return []
    print(aFile, bFile, f'{value} > {cutoff}')
    newVolume = {
        'volume': None,
        'polygon': '',
        'point': [.5 * (aVolume['point'][i] + bVolume['point'][i]) for i in range(len(aVolume['point']))]
    }
    volumes.append(newVolume)
    # returm new set of volumes to evaluate
    return [(aVolume, newVolume), (newVolume, bVolume)]


if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('--out', dest='out', required=True, type=str)
    parser.add_argument('--volumes', dest='volumes', default=None, type=str)
    parser.add_argument('--points', dest='points', default=None, type=str)
    parser.add_argument('--map', dest='map', required=True, type=str)
    parser.add_argument('--radius', dest='radius', required=True, type=float)
    parser.add_argument('--cutoff', dest='cutoff', required=True, type=float)
    i = sys.argv.index('--')
    args = parser.parse_args(sys.argv[i + 1:])
    try:
        main(args)
    except Exception as e:
        print('error', e)
