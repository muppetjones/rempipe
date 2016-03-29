'''A series of functions for setting up an argparse parser

Functions for setting up and adding to a single parser.

'''

import argparse

from libpipe.util import path

import logging
log = logging.getLogger(__name__)

#
#   Custom argparse actions
#


class ProtectAbsPathArg(argparse.Action):

    '''Protect path string'''

    def __call__(self, parser, namespace, dir_str, option_string):
        dir_str = path.protect(dir_str)
        setattr(namespace, self.dest, dir_str)


class ProtectAbsPathList(argparse._AppendAction):

    '''Protect path strings'''

    def __call__(self, parser, namespace, values, option_string):
        values = path.protect(values)
        super(ProtectAbsPathList, self).__call__(
            parser, namespace, values, option_string)


class ParseSummaryArg(argparse.Action):

    def __call__(self, parser, namespace, summary_file, option_string):
        summary_file = path.protect(summary_file)
        summary = {}
        with open(summary_file, 'r') as fh:
            rows = [line.lstrip().rstrip().split() for line in fh]
            summary = {col[0]: col[1].split(';') for col in rows}

        setattr(namespace, self.dest, summary)


#
#   Define input argparse objects
#


def pipe_parser(parser=None):

    if not parser:
        parser = argparse.ArgumentParser()

    __add_pipe_group(parser)

    return parser


def __add_pipe_group(parser):

    grp = parser.add_argument_group('Pipe options')
    grp.add_argument(
        '--project', metavar='NAME', dest='project',
        action='store',
        help='The name of the project. Also the name of the output directory.'
    )

    grp.add_argument(
        '--root', metavar='PATH', dest='root',
        action=ProtectAbsPathArg,
        help='The root output directory. Data will be stored in root/project.'
    )

    grp.add_argument(
        '--summary', metavar='FILE', dest='summary',
        action=ParseSummaryArg,
        help=' '.join([
            'Path to file with project summary.',
            'Expects first two columns in format:',
            '<name> <comma-separated filenames>.',
        ])
    )

    grp.add_argument(
        '--data', metavar='PATH', dest='data',
        action=ProtectAbsPathArg,
        help='The dir where the data files given in --summary are stored.'
    )

    grp.add_argument(
        '--genome', metavar='INDEX', dest='genome',
        action=ProtectAbsPathArg,
        help=('The name of the alignment index, including path.'),
    )

    grp.add_argument(
        '--filter', metavar='INDEX', dest='filter_list',
        action=ProtectAbsPathList,
        help=('Other genome indices to use to filter unwanted reads.'),
    )

    grp.add_argument(
        '--debug', dest='debug',
        action='store_true',
        default=False,
        help=('Run in debug mode. Prevent execution of pipe if set'),
    )

    return
