
import abc
import argparse
import os.path
import re

import remsci.scripted.base as base
import remsci.lib.utility.path as path
from remsci.scripted.interface import SubparserBase

import logging
log = logging.getLogger(__name__)


class BasePipeParser(SubparserBase, metaclass=abc.ABCMeta):

    '''Creates a subparser with several predefined arguments for pipes

    Arguments:
        --paired-end [PATTERN]  Specify that the data is paired-end.
            Optionally provde a pattern for identifying pairs.
        --project CHAR  Specify the project name. Must be unique.
        --extend        Add to existing project (allows non-unique names).
        --summary FILE  Summary file with sample names and paths.
        --root PATH     Path to root directory.
        --data PATH     Specify location of summary file and data.
                        May omit if summary is an absolute path.
        --force         Force all commands to re-run (overwrite existing data).
        --debug         Write files, but do not run.

    Child must define:
        pipe attribute: The pipe class to use.
        subcmd attribute: The name of the subcommand to use when calling.
        setup method: CALL SUPER! Must implement to differenciate between
            children parsers.

    '''

    def __init__(self, subparser=None):
        self.subparser = subparser.add_parser(
            self.subcmd,
            parents=[
                base.input_file_parser(),
                base.input_directory_parser(),
                base.output_file_parser(),
            ],
            formatter_class=argparse.ArgumentDefaultsHelpFormatter,
            help='Run fastqc on a fastq file',
        )

    @abc.abstractmethod
    def setup(self):
        super().setup()

        self.subparser.add_argument(
            '--paired_end', dest='paired_end',
            nargs='?', const=r'_2\.f', type=str,
            metavar='PATTERN',
            help=(
                'Specify paired end data.'
                + 'Optional: Pattern for identification of 2nd pair.'),
        )

        # unique project name
        self.subparser.add_argument(
            '--project', dest='project',
            default='new_project',
            metavar='NAME',
            help='Name of project. Must be unique unless "--extend" given.',
        )

        # add files to project (allows for non-unique name)
        # -- do we actually use this?
        self.subparser.add_argument(
            '--extend', dest='extend',
            action='store_true',
            help='Flag to indicate that files should be added to project',
        )

        self.subparser.add_argument(
            '--summary', dest='summary_file',
            metavar='FILE',
            help='Path to file with project summary',
        )

        self.subparser.add_argument(
            '--root', dest='root_dir',
            help='root directory',
        )

        self.subparser.add_argument(
            '--data', dest='data_dir',
            help='data directory. Should contain fastq files and summary.',
        )

        self.subparser.add_argument(
            '--force', dest='force',
            action='store_true',
            default=False,
            help=('Force pipe to rerun steps. ' +
                  'Default: Skip step if output found.'),
        )

        self.subparser.add_argument(
            '--debug', dest='debug',
            action='store_true',
            default=False,
            help=('Do not run pipe if set'),
        )

        # only fastq files by default
        self.subparser.set_defaults(
            pattern_list=[r'\.f(?:ast)?q', r'fastq\.gz'])

        return

    def run(self, args):
        # protect filenames
        protect = ['summary_file', 'root_dir', 'data_dir']
        for p in protect:
            try:
                setattr(args, p, path.protect(getattr(args, p)))
            except TypeError:
                pass  # arg might not be set (is None)

        # read summary file
        if args.summary_file is not None:
            setattr(args, 'summary', self.read_summary(args))
        else:
            args.find_files()  # sets file_list attribute

        # add project directory attribute
        project_dir = os.path.join(args.root_dir, args.project, 'samples')
        setattr(args, 'project_dir', project_dir)

        # add pipe attribute
        # -- use to create pipes!
        setattr(args, 'pipe', self.pipe)

        return args

    def read_summary(self, args):
        '''Reads a summary file and returns a list of names and files

        Expected format:
                <name>  <file>  [<file> ...]    [<data> ...]
                <name>  <file>[;<file>]         [<data> ...]

        The first column should be the name of the file. The files,
        should be FASTA or FASTQ format and listed in the following
        columns (or, preferably, semi-colon separated in a single
        column). The columns should be space, tab, or comma separated.
        Semi-colons should only be used within a single column.

        Multiple, consecutive delimiters will be treated as one delimiter.
        Only files to be included should be listed; multiple files should
        denote paired-end or mate-pair reads.
        '''

        with open(args.summary_file, 'r') as fh:
            rows = [re.split(r'[;,\s]+', line.rstrip()) for line in fh]
            data_dir = (args.data_dir
                        if args.data_dir else os.path.dirname(fh.name))

        # compile regex
        rx_seq = re.compile(r'\.f(?:ast)?[aq]$')

        # add and protect path to second (and third) columns in rows
        for row in rows:
            row[1:] = [
                path.protect(os.path.join(data_dir, col))
                for col in row[1:]
                if rx_seq.search(col.lower())
            ]

        return rows
