
import os.path

from libpipe.cmd.attr import CmdAttributes
from libpipe.cmd.base import CmdBase
from libpipe.type import index as _index
from libpipe.util import path

import logging
log = logging.getLogger(__name__)

#
#   Indices
#

Hisat2Index = _index.IndexType.factory({'.ht2': 8}, name='Hisat2Index')

#
#   Commands
#


class Hisat2Cmd(CmdBase):

    '''Hisat2 Sequence Aligner

    Command to run the hisat2 aligner.

    TODO(sjbush): Add stricter formatting to output (sam, log, etc.) by
        adding timestamp, job_name, etc.
    TODO(sjbush): Enable caching of check index results.
    TODO(sjbush): Add AlignIndex class and factory.

    Command usage:
        hisat2 [options]* -x <bt2-idx> {-1 <m1> -2 <m2> | -U <r>} -S <sam>

    '''

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
            [('-x', ), (Hisat2Index, )]
        ],
    )

    #
    #   Cmd Interface
    #

    def output(self):
        _output = [self.kwargs['-S'], ]
        try:
            _output.append(self.kwargs['--un'])
        except KeyError:
            prefix, extn = os.path.splitext(self.kwargs['--un-conc'])
            _output.extend([
                '{}.{}{}'.format(prefix, i + 1, '.fastq')
                for i in range(2)
            ])
        return _output

    #
    #   Overrides
    #

    def _match_input_with_args(self):
        print('^' * 500)
        try:
            super()._match_input_with_args()
        except Exception as e:
            print('*' * 500)
            log.debug(str(e))
            raise
        else:
            print('V' * 500)

        unused_input = [
            _input for _input in self.input()
            if _input not in self.args and _input not in self.kwargs.values()
        ]
        log.debug(unused_input)

        # for unused in unused_input:
        #     try:
        #         _index = self.index(unused)
        #         self._check_for_index_files(unused)
        #     except ValueError:
        #         pass  # not a valid index
        #     else:
        #         log.debug(_index)
        #         log.debug(type(_index))
        #         self.kwargs['-x'] = _index  # unused
        #         break  # only need one!

    def _pre_cmd(self):

        if '-S' not in self.kwargs:
            try:
                sam_file = '{}.sam'.format(
                    os.path.commonprefix(
                        [self.kwargs['-1'], self.kwargs['-2']]).rstrip('._')
                )
            except KeyError:
                sam_file = '{}.sam'.format(
                    os.path.splitext(self.kwargs['-U'])[0]
                )
            self.kwargs['-S'] = sam_file

        # basic log file name
        log_file = self.kwargs['-S'].replace('.sam', '.log')
        log_file = os.path.join(
            os.path.dirname(log_file),
            'hisat_{}'.format(os.path.basename(log_file))
        )
        self.redirect = ('2>&1', '|', 'tee -a', log_file)

        # ensure unaligned reads are written to a file
        # -- tag with _unal label (but make sure there's only one)
        unal_key = '--un' if '-U' in self.kwargs else '--un-conc'
        unal_base = os.path.splitext(self.kwargs['-S'])[0].replace('_unal', '')
        unal = unal_base + '_unal.fastq'
        self.kwargs.update({unal_key: unal})

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
