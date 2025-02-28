import argparse
from targettracking.configuration import ConfigurationFactory
from targettracking.cli import Volumes, Circles
import logging
from logging.handlers import RotatingFileHandler
import sys
import os


def main():
    subCommands = [
        Volumes(),
        Circles()
    ]
    parser = argparse.ArgumentParser(prog='targettracking')
    parser.add_argument('-c', '--config', type=str, default='./targettracking/config.yaml')
    subparsers = parser.add_subparsers()
    for subCommand in subCommands:
        subparser = subparsers.add_parser(subCommand.key)
        subCommand.modifyParser(subparser)
    args = parser.parse_args()
    config = ConfigurationFactory.getInstance(args.config)['logging']
    logger = logging.getLogger('targettracking')
    for log in config['handlers']:
        if log['out'] == 'stdout':
            stream = sys.stdout
            handler = logging.StreamHandler(stream)
        elif log['out'] == 'stderr':
            stream = sys.stderr
            handler = logging.StreamHandler(stream)
        else:
            handler = RotatingFileHandler(log['out'], log['mode'], maxBytes=100 * 2 ** 20, backupCount=3)
        form = logging.Formatter(log['format'])
        handler.setLevel(log['level'])
        handler.setFormatter(form)
        logger.addHandler(handler)
    logger.setLevel('DEBUG')
    try:
        args.application(args)
    except Exception as e:
        logger.critical('Critial Error', exc_info=True)


if __name__ == '__main__':
    main()
