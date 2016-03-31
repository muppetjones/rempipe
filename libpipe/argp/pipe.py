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
    '''Custom argparse.Action for reading a summary file.

    [Deprecated] Read in a summary file and set `summary` with
    a dict of the results.

    CRITICAL: This is what NOT to do. It works, but it's poor design.
        It will be difficult to maintain, and this design pattern
        will result in MANY specialized actions.
    '''

    def __call__(self, parser, namespace, summary_file, option_string):
        raise(DeprecationWarning('Do not use this argparse.Action'))
        summary_file = path.protect(summary_file)
        summary = {}
        with open(summary_file, 'r') as fh:
            rows = [line.lstrip().rstrip().split() for line in fh]
            summary = {col[0]: col[1].split(';') for col in rows}

        setattr(namespace, self.dest, summary)


#
#   Define pipe argparse objects
#


def pipe_parser(parser=None):

    if not parser:
        parser = argparse.ArgumentParser()

    __add_pipe_group(parser)

    return parser


def build_args(args):
    '''Perform additional parsing of given args.

    Rather than create a dozen specialized argparse.Action objects,
    of which ParseSummaryArg is an example, separate that functionality
    into individual functions.

    `build_args` should provide a single interface to these specialized
    functions, removing the need to call each separately.

    Specializations:
    *   Update `file_list`: Read the given directories and append all
        of the found files to the file_list attribute.
    *   Update `summary`: Read the given summary file and create a dict
        of the files with the format:
            {'name': [file, ...]}
    '''

    # TODO(sjbush): `file_list` is not defined here, so build_args
    #   should not handle this! Make a registration.

    args = __read_summary(args)

    return args

#
#   Private functions
#


def __add_pipe_group(parser):

    grp = parser.add_argument_group('Pipe options')

    # BASIC
    grp.add_argument(
        '--project', metavar='NAME', dest='project',
        action='store',
        help='The name of the project. Also the name of the output directory.'
    )

    # PATH
    grp.add_argument(
        '--root', metavar='PATH', dest='root',
        action=ProtectAbsPathArg,
        help='The root output directory. Data will be stored in root/project.'
    )

    grp.add_argument(
        '--summary', metavar='FILE', dest='summary',
        action=ProtectAbsPathArg,
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

    # PATH LIST
    grp.add_argument(
        '--genome', metavar='INDEX', dest='genome_list',
        action=ProtectAbsPathList,
        help=('The name of one or more alignment indices, including path.'),
    )

    grp.add_argument(
        '--filter', metavar='INDEX', dest='filter_list',
        action=ProtectAbsPathList,
        help=('Other genome indices to use to filter unwanted reads.'),
    )

    # FLAGS
    grp.add_argument(
        '--debug', dest='debug',
        action='store_true',
        default=False,
        help=('Run in debug mode. Prevent execution of pipe if set'),
    )

    grp.set_defaults(build_args=build_args)

    return


#
#   Arg handling
#

def __read_summary(args):
    setattr(args, 'summary_file', args.summary)
    with open(args.summary_file, 'r') as fh:
        rows = [line.lstrip().rstrip().split() for line in fh]
        summary = {col[0]: col[1].split(';') for col in rows}
    setattr(args, 'summary', summary)
    return args