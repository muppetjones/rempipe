
import os.path
import sys
import argparse

from libpipe.utility import path

import logging
log = logging.getLogger(__name__)

# =============================================================================
# Declare parsers
# =============================================================================


class DummyParser(argparse.ArgumentParser):

    def error(self, message):
        log.error('error: %s\n\n' % message)
        self.print_help()
        sys.exit(message)


class RemScripted(object):

    parser = DummyParser(
        description="Scripting interface for libpiperary",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    subparsers = parser.add_subparsers(
        title='Commands',
        description='Subcommands available for accessing remsci',
        dest='command',
        help="additional help",
    )

    @classmethod
    def __reset__(cls):

        cls.parser = DummyParser(
            description="Scripting interface for libpiperary",
            formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        )
        cls.subparsers = cls.parser.add_subparsers(
            title='Commands',
            description='Subcommands available for accessing remsci',
            dest='command',
            help="additional help",
        )


# =============================================================================
# Custom argparse classes
# =============================================================================


class UniqueStore(argparse.Action):

    '''http://stackoverflow.com/questions/23032514/
        argparse-disable-same-argument-occurences'''

    def __call__(self, parser, namespace, values, option_string):
        # log.debug(namespace)
        # log.debug(self.dest)
        # log.debug(self.default)
        if getattr(namespace, self.dest, self.default) is not None:
            parser.error(option_string + " appears several times.")
        setattr(namespace, self.dest, values)


class ProtectFileList(argparse._AppendAction):

    '''Protect file strings'''

    def __call__(self, parser, namespace, values, option_string):
        values = path.protect(values)
        super(ProtectFileList, self).__call__(
            parser, namespace, values, option_string)


# =============================================================================
# Parser access
# =============================================================================


def get_parser():
    return RemScripted.parser


def get_subparser():
    return RemScripted.subparsers


def input_file_parser():
    '''Create a parser with options for file input

    Options:
        -f, --file <file>
            File(s) for input. Can be given multiple times.
        --header, --no-header
            States whether or not a column header is expected (Default: True).
        --delim
            Expected column delimiter.

    Examples:
        > myprog.py --file file1.txt --file file2.txt
    '''

    input_parser = argparse.ArgumentParser(
        description="File input options",
        add_help=False)

    file_group = input_parser.add_argument_group(
        'File input options')
    file_group.add_argument('-f', '--file', dest='file_list',
                            metavar='FILE',
                            action=ProtectFileList,
                            help='The input file path(s).')
    file_group.add_argument(
        '--header', dest='header', action='store_true',
        help='Denotes whether the file contains column names'
    )
    file_group.add_argument('--no-header', dest='header', action='store_false')
    file_group.set_defaults(header=True)
    file_group.add_argument(
        '--delim', dest='delim',
        default=argparse.SUPPRESS,
        help='Column delimiter in file.',
    )
    return input_parser


def input_directory_parser():
    '''Parser with options for directory input_parser

    Search a given directory for files matching a given pattern, extension,
    or both. Patterns and extensions may be used together, but be careful
    that they do not cancel each other out.

    Options:
        -d, --directory (string): Directory where input files are stored (depth
            is irrelevant if you use the `get_file_list_from_directory`
            function).
        -e, --extension (multiple, string): Input file extension(s).
        --pattern (multiple, regex string): Input file pattern(s).

    Examples:
        # search for fastq files from controls in a given directory
        > myprog.py --directory "~/data" --pattern 'control' -e '.fq'
    '''

    input_parser = argparse.ArgumentParser(
        description="Directory input options",
        add_help=False)
    dir_group = input_parser.add_argument_group(
        'Directory input options')
    dir_group.add_argument('-d', '--directory', dest='dir_path',
                           metavar='directory'.upper(),
                           action=UniqueStore,
                           help='The directory containing input files')
    dir_group.add_argument('-e', '--extension', dest='extension_list',
                           metavar='extension'.upper(),
                           default=[],
                           action='append',
                           help='The desired file extension',
                           )
    dir_group.add_argument('--pattern', dest='pattern_list',
                           metavar='pattern'.upper(),
                           default=[],
                           action='append',
                           help='A substring to look for in files',
                           )
    dir_group.set_defaults(find_files=_get_file_list_from_directory)
    return input_parser


def output_file_parser():
    '''Return parser with options for output files'''
    output_parser = argparse.ArgumentParser(
        description="Output file options",
        add_help=False)

    group = output_parser.add_argument_group(
        'Output options')
    group.add_argument('-o', '--outfile', dest='outfile',
                       metavar='file'.upper(),
                       action=UniqueStore,
                       help='The output file',
                       )
    return output_parser


def output_directory_parser():
    output_parser = argparse.ArgumentParser(
        description="Output directory options",
        add_help=False)
    group = output_parser.add_argument_group(
        'Output', 'options for output')
    group.add_argument('-o', '--outdir', dest='outdir',
                       metavar='directory'.upper(),
                       action=UniqueStore,
                       help='The output directory',
                       )
    return output_parser

# =============================================================================
# Scripting functions
# =============================================================================


def _get_file_list_from_directory(args):
    '''Return a list of files from the command line options

    NOTE: This might be redundant...perhaps make remove and let each function
          do it for themselves?
    '''

    if not args.dir_path:
        return

    if not args.extension_list and not args.pattern_list:
        msg = "Use of '-d' also requires either '--extension' or '--pattern'"
        raise argparse.ArgumentError(None, msg)
    file_list = path.walk_file(args.dir_path,
                               extension=args.extension_list,
                               pattern=args.pattern_list,
                               )

    log.info("Working with {0} files in \n   {1}".format(
        len(file_list),
        os.path.commonprefix(file_list),
    ))

    setattr(args, 'file_list', file_list)
    del args.dir_path
    del args.extension_list
    del args.pattern_list
    return file_list


def _get_files(args):
    '''Return a list of files from the command line options'''
    try:
        return [path.protect(f) for f in args.file_path]
    except AttributeError:
        if not args.file_extension and not args.file_pattern:
            raise ValueError('Require file extension or pattern')
        file_list = path.walk_file(args.dir_path,
                                   extension=args.file_extension,
                                   # level=5,
                                   # patterns=[r'\.bam'],
                                   pattern=args.file_pattern,
                                   )
    except:
        raise

    log.info("Working with {0} files in \n   {1}".format(
        len(file_list),
        path.common_directory(file_list),
    ))

    return file_list
