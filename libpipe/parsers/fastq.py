# import argparse
import os.path
import re

# import libpipe.parsers.pipe as base
# from libpipe.parsers.subparser import SubparserBase
# from libpipe.utility import path
# from libpipe.decorators import file_or_handle

from libpipe.parsers.pipe import BasePipeParser
from libpipe.pipes.assemble import WgsPipe
from libpipe.pipes.align import RnaSeqPipe

import logging
log = logging.getLogger(__name__)


class WgsPipeParser(BasePipeParser):

    subcmd = 'wgs'
    pipe = WgsPipe

    def setup(self):
        super().setup()

        self.subparser.add_argument(
            '--filter', dest='filter_genomes',
            action='append',
            help='Align to genome and use unaligned sequences',
        )

        self.subparser.add_argument(
            '--reference', dest='reference_genome',
            help='The reference genome.',
        )

    def run(self, args):
        return super().run(args)


class RnaseqPipeScripted(BasePipeParser):

    subcmd = 'rnaseq'
    pipe = RnaSeqPipe

    def setup(self):
        super().setup()

        self.subparser.add_argument(
            '--genome', dest='genome',
            action='append',
            help='Hisat genome',
        )

    def run(self, args):
        if args.summary is not None:
            return self._get_files_from_summary(args.summary, args.data_dir)

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
