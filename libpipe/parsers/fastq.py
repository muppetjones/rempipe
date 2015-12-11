import argparse
import os.path
import re

import remsci.scripted.base as base
from remsci.scripted.interface import SubparserBase
from remsci.lib.tree import RemTree
from remsci.lib.utility import path
from remsci.lib.decorators import file_or_handle

import logging
log = logging.getLogger(__name__)


class FastqScripted(SubparserBase):

    def __init__(self, subparser=None):
        # create the subparser
        self.subparser = subparser.add_parser(
            'fastq',
            parents=[
                base.input_file_parser(),
                base.input_directory_parser(),
                base.output_file_parser(),
            ],
            formatter_class=argparse.ArgumentDefaultsHelpFormatter,
            help='Run fastqc on a fastq file',
        )

    def setup(self):
        super().setup()

        self.subparser.add_argument(
            '--paired_end', dest='paired_end',
            nargs='?', const=r'_2\.f', type=str,
            metavar='PATTERN',
            help='Specify paired end data. Optional: Pattern for identification of 2nd pair.',
        )

        self.subparser.add_argument(
            '--project', dest='project',
            metavar='NAME',
            help='Name of project. Must be unique unless "--extend" given.',
        )

        self.subparser.add_argument(
            '--summary', dest='summary',
            metavar='FILE',
            help='Path to file with project summary',
        )

        self.subparser.add_argument(
            '--extend', dest='extend',
            action='store_true',
            help='Flag to indicate that files should be added to project',
        )

        self.subparser.add_argument(
            '--genome', dest='genome',
            help='Hisat genome',
        )

        self.subparser.add_argument(
            '--root', dest='root',
            help='root directory',
        )

        # only fastq files by default
        self.subparser.set_defaults(
            pattern_list=[r'\.f(?:ast)?q', r'fastq\.gz'])

        return

    def run(self, args):
        if args.summary is not None:
            return self._get_files_from_summary(args.summary)

        if args.file_list is None:
            raise ValueError('No fastq files found')

        # separate paired end
        try:
            rx_pe = re.compile(args.paired_end)
        except TypeError:
            raise NotImplementedError('Use --paired_end option')
        else:
            r1 = [f for f in args.file_list if rx_pe.search(f) is None]
            r2 = [f for f in args.file_list if f not in r1]
            pe_files = list(zip(r1, r2))
            if not pe_files:
                raise ValueError('Input is not paired end.')
            log.debug([(os.path.basename(f[0]), os.path.basename(f[1]))
                       for f in pe_files])
            return pe_files

    @file_or_handle
    def _get_files_from_summary(self, fh):
        '''Reads a summary file and returns a list of names and files

        Expected format:
                <name>  <file>  [<file> ...]
            The first two columns should be the name of the sample
            and the file. Paired-end should be listed on the same line.
        '''

        files = []
        rows = [line.rstrip().split() for line in fh]
        names = [cols[0] for cols in rows]
        r1 = [path.protect(cols[1]) for cols in rows]

        # check for valid paired end in third column
        try:
            r2 = [path.protect(cols[2]) for cols in rows]
            r2_0 = os.path.join(os.path.dirname(fh.name), r2[0])
            if not os.path.isfile(r2[0]) and not os.path.isfile(r2_0):
                raise IndexError()
        except IndexError:
            pass

        try:
            return list(zip(names, r1, r2))
        except UnboundLocalError:
            return list(zip(names, r1))
