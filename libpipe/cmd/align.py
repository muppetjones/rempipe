
import os.path

from libpipe.cmd.attr import CmdAttributes
from libpipe.cmd.base import CmdBase
from libpipe.util import path

import logging
log = logging.getLogger(__name__)


class Hisat2Cmd(CmdBase):

    attr = CmdAttributes(
        invoke='hisat2',
        args=[
            ('-x', 'FILE', 'Hisat reference genome index (base name)'),
            ('-1', 'FILE[,FILE]',
             'comma separated list of paired-end 1 files'),
            ('-2', 'FILE[,FILE]',
             'comma separated list of paired-end 2 files'),
            ('-U', 'FILE[,FILE]', 'comma separated list of unpaired reads'),
            ('-S', 'FILE', 'Output sam file (defaults to read prefix)'),
            ('-p', 'INT', 'number of processors'),
            ('-I', 'INT', 'minimum fragment length. Default = 0.'),
            ('-X', 'INT', 'aximum fragment length. Default = 500.'),
            ('--un-conc', 'PATH',
             'Path to write unaligned, paired-end reads to.'),
            ('--phred33', None, 'Illumina 1.9+ encoding'),
            ('--phred64', None, 'Illumina 1.8 and earlier encoding'),
            ('--fr', None, 'Upstream downstream mate orientations'),
            ('-q', None, 'Reads are FASTQ files'),
            ('-f', None, 'Reads are FASTA files'),
        ],
        defaults={
            '-p': 3,  # "$(wc -l < $PBS_NODEFILE)",
            '-I': 0,
            '-X': 500,
        },

        req_kwargs=['-x', ('-1', '-2'), ['-1', '-U']],
        req_args=0,
        req_types=[
            [('-1', '-2'), ('.fastq', '.fq', '.fastq', '.fa'), True],
            [('-U', ), ('.fastq', '.fq', '.fastq', '.fa'), True],
            [('-S', ), ('.sam', )],
        ],
    )

    #
    #   Cmd Interface
    #

    def output(self):
        return [self.kwargs['-S'], ]

    #
    #   Overrides
    #

    def _match_input_with_args(self):
        super()._match_input_with_args()

        unused_input = [
            _input for _input in self.input()
            if _input not in self.args and _input not in self.kwargs.values()
        ]

        for unused in unused_input:
            try:
                self._check_for_index_files(unused)
            except ValueError:
                pass  # not a valid index
            else:
                self.kwargs['-x'] = unused
                break  # only need one!

    def _pre_cmd(self):

        if '-S' not in self.kwargs:
            try:
                sam_file = '{}.sam'.format(
                    os.path.commonprefix(
                        [self.kwargs['-1'], self.kwargs['-2']]).rstrip('.')
                )
            except KeyError:
                sam_file = '{}.sam'.format(
                    os.path.splitext(self.kwargs['-U'])[0]
                )
            self.kwargs['-S'] = sam_file

    #
    #   Class methods
    #

    @classmethod
    def _check_for_index_files(
            cls, index_name, expected_file_count=8, extension='.ht2'):
        # ensure the index exists
        genome_dir, genome_base = os.path.split(index_name)
        if not genome_dir:
            genome_dir = './'

        index_pattern = r'{}\..*{}'.format(genome_base, extension)
        index_files = path.walk_file(genome_dir, pattern=index_pattern)

        if len(index_files) != expected_file_count:
            msg = 'Fewer index files ({}) than expected ({})'.format(
                len(index_files), expected_file_count)
            raise ValueError(msg)
