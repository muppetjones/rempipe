import argparse
import os.path
import re

from argparse import ArgumentError

import remsci.scripted.base as base
from remsci.scripted.interface import SubparserBase
from libpipe.utility import path
from libpipe.decorators import file_or_handle

import logging
log = logging.getLogger(__name__)


class AggregateScripted(SubparserBase):

    def __init__(self, subparser=None):

        if not subparser:
            subparser = base.get_subparser()

        # create the subparser
        self.subparser = subparser.add_parser(
            'aggregate',
            parents=[
                base.input_file_parser(),
                base.input_directory_parser(),
                base.output_file_parser(),
            ],
            formatter_class=argparse.ArgumentDefaultsHelpFormatter,
            help='Aggregate output summaries',
        )

        self.setup()

    def setup(self):
        super().setup()

        self.subparser.add_argument(
            '--fastqc', dest='fastqc',
            action='store_true',
            default=False,
            help=('Aggregate FastQC output'),
        )

        self.subparser.add_argument(
            '--align', dest='alignment_summary',
            action='store_true',
            default=False,
            help=('Aggregate alignment output'),
        )

        self.subparser.add_argument(
            '--count', dest='count_summary',
            action='store_true',
            default=False,
            help=('Aggregate count output'),
        )

        self.subparser.set_defaults(
            pattern_list=[
                r'\.zip',  # fastqc summary archives
                r'\.count',  # count files
                r'hisat2?\.log',  # hisat output
            ])

        return

    def run(self, args):

        # only use given patterns (i.e., erase the default 3)
        if len(args.pattern_list) > 3:
            del args.pattern_list[0:3]  # remember, up to, not including

        # parse directory file list
        # (requires input_directory_parser to be included)
        try:
            args.find_files(args)
        except ArgumentError as e:
            print(e.message)
            exit()
        except AttributeError:
            pass

        return args.file_list
