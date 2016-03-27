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


class AddLeadingPeriodList(argparse._AppendAction):

    '''Ensure leading period'''

    def __call__(self, parser, namespace, val, option_string):
        val = val if val.startswith('.') else '.{}'.format(val)
        super().__call__(
            parser, namespace, val, option_string)


class ProtectRelPathList(argparse._AppendAction):

    '''Protect path strings'''

    def __call__(self, parser, namespace, values, option_string):
        values = path.protect(values, abspath=False)
        super(ProtectRelPathList, self).__call__(
            parser, namespace, values, option_string)


class ProtectAbsPathList(argparse._AppendAction):

    '''Protect path strings'''

    def __call__(self, parser, namespace, values, option_string):
        values = path.protect(values)
        super(ProtectAbsPathList, self).__call__(
            parser, namespace, values, option_string)


class RemoveEnclosingQuotesList(argparse._AppendAction):

    '''Ensure leading period'''

    def __call__(self, parser, namespace, val, option_string):
        val = val.lstrip('"\'').rstrip('"\'')
        super().__call__(
            parser, namespace, val, option_string)


class UniqueStore(argparse.Action):

    '''Ensure argument is accepted only once.

    http://stackoverflow.com/questions/23032514/
    argparse-disable-same-argument-occurences
    '''

    def __call__(self, parser, namespace, values, option_string):
        if getattr(namespace, self.dest, self.default) is not None:
            parser.error(option_string + " appears several times.")
        setattr(namespace, self.dest, values)


#
#   Define input argparse objects
#

# TODO(sjbush): Make dir feed directly into file_list

def input_parser(parser=None, accept_dirs=True, accept_files=True):

    if not parser:
        parser = argparse.ArgumentParser()

    if accept_files:
        __add_file_group(parser)

    if accept_dirs:
        __add_dir_group(parser)

    return parser


def __add_file_group(parser):
    grp = parser.add_argument_group('File input options')
    grp.add_argument(
        '-f', '--file', metavar='FILE', dest='file_list',
        action=ProtectAbsPathList,
        help='One or more input files.'
    )
    return


def __add_dir_group(parser):
    grp = parser.add_argument_group('Directory input options')
    grp.add_argument(
        '-d', '--dir', metavar='DIR', dest='dir_list',
        action=ProtectAbsPathList,
        help='One or more input directories.'
    )

    grp.add_argument(
        '--extn', metavar='EXTN', dest='extn_list',
        action=AddLeadingPeriodList,
        help='Defines required extension for files in -d DIR.'
    )

    grp.add_argument(
        '--pattern', metavar='REGEX', dest='pattern_list',
        action=RemoveEnclosingQuotesList,
        help='Defines regex patterns for filtering files in -d DIR.'
    )

    grp.set_defaults(find_files=__append_dir_files_to_file_list)
    return


#
#   output argparse functions
#

def output_parser(parser=None, accept_files=True, accept_dirs=False):
    '''Define output arguments on an ArgumentParser

    Adds '-o' and '--out' arguments for file output and '--odir' argument
    for directory output.

    Arguments:
        parser: ArgumentParser object. If None (default), creates new object.
        accept_files: Bool indicating whether or not to define
            an out file argument. Default=True.
        accept_dirs: Bool indicating whether or not to define
            an out dir argument. Default=False,
    '''

    if not parser:
        parser = argparse.ArgumentParser()

    par_grp = parser.add_argument_group('Output options')

    if accept_files:
        __add_outgroup_file(par_grp)

    if accept_dirs:
        __add_outgroup_dir(par_grp)

    return parser


def __add_outgroup_file(parser):

    parser.add_argument(
        '-o', '--outfile', dest='ofile',
        metavar='FILE', action=UniqueStore,
        help='The output file.',
    )
    return


def __add_outgroup_dir(parser):

    parser.add_argument(
        '--odir', dest='odir',
        metavar='DIR', action=UniqueStore,
        help='The directory to store output in.',
    )
    return


#
#   Helper functions
#


def __append_dir_files_to_file_list(args):
    '''Append files found within given directories to args.file_list

    Finds files in each given directory and appends those files to the
    file_list attribute. If extensions or patterns are defined, the
    list is filtered accordingly
    '''

    if not args.dir_list:
        return

    extra_files = []
    for _dir in args.dir_list:
        found_files = path.walk_file(
            _dir,
            extension=args.extn_list,
            pattern=args.pattern_list,
        )
        extra_files.extend(found_files)

    # add all discovered files to file_list
    try:
        args.file_list.extend(extra_files)
    except AttributeError:
        setattr(args, 'file_list', extra_files)

    # delete the directories--we don't want re-add!
    del args.dir_list
    return args.file_list
