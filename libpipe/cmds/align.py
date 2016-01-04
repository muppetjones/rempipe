
import os.path
import time

from libpipe.cmds.base import BaseCmd

import logging
log = logging.getLogger(__name__)


class HisatCmd(BaseCmd):

    NAME = 'hisat'
    INVOKE_STR = 'hisat'

    ARGUMENTS = [
        ('-x', 'FILE', 'Hisat reference genome index (base name)'),
        ('-1', 'FILE[,FILE]', 'comma separated list of paired-end 1 files'),
        ('-2', 'FILE[,FILE]', 'comma separated list of paired-end 2 files'),
        ('-U', 'FILE[,FILE]', 'comma separated list of unpaired reads'),
        ('-S', 'FILE', 'output sam file'),
        ('-p', 'INT', 'number of processors'),
        ('-I', 'INT', 'minimum fragment length. Default = 0.'),
        ('-X', 'INT', 'aximum fragment length. Default = 500.'),
        ('--un-conc', 'PATH', 'Path to write unaligned, paired-end reads to.'),
        ('--phred33', None, 'Illumina 1.9+ encoding'),
        ('--phred64', None, 'Illumina 1.8 and earlier encoding'),
        ('--fr', None, 'Upstream downstream mate orientations'),
        ('-q', None, 'Reads are FASTQ files'),
        ('-f', None, 'Reads are FASTA files'),
    ]
    DEFAULTS = {
        '-p': 3,  # "$(wc -l < $PBS_NODEFILE)",
        '-I': 0,
        '-X': 500,
    }

    REQ_KWARGS = ['-x', ('-1', '-2'), ['-1', '-S']]
    REQ_ARGS = 0
    REQ_TYPE = [
        [('-1', '-2'), ('.fastq', '.fq', '.fastq', '.fa'), False],
        [('-U', ), ('.fastq', '.fq', '.fastq', '.fa'), False],
        [('-S', ), ('sam', )],
    ]

    #
    #   Magic methods
    #

    def __init__(
            self, *args,
            encoding='--phred33', orientation='--fr', format='-q',
            **kwargs):

        try:
            super().__init__(*args, **kwargs)
        except ValueError:
            raise   # requirements failure; pass it on up

        # update flags
        self.flags.extend([
            encoding, orientation
        ])

        # set the timestamp if not done already
        if not self.timestamp:
            self.timestamp = time.strftime("%y%m%d-%H%M%S")

    #
    #   "Public" methods
    #

    def output(self):
        return [self.kwargs['-S']]

    #
    #   "Private" methods
    #

    # def _input(self):
    #     '''Match linked output to sequence inputs'''
    #
    #     args = self._get_input()
    #     if not args:
    #         return
    #
    #     fq_types = ('.fastq', '.fq', '.fastq', '.fa')
    #     filtered = self._filter_by_type(args, fq_types)
    #
    #     try:
    #         self.kwargs['-1'], self.kwargs['-2'] = filtered
    #     except ValueError:
    #         try:
    #             self.kwargs['-U'], = filtered
    #         except ValueError:
    #             msg = 'Unknown input from link: {} ({})'.format(
    #                 args, self.input.__self__.name)
    #             raise self.CmdLinkError(msg)

    def _prepcmd(self):
        '''Prep for hisat cmd

        > parse log file name and set for redirect
        > ensure unaligned reads are output
        '''

        # parse log file name (based on given file info)
        # sample_name: the basename of the given file
        try:  # single end reads
            log_dir = os.path.dirname(self.kwargs['-U'])
            sample_name = self._trubase(self.kwargs['-U'])
            unal_key = '--un'
        except KeyError:  # paired-end reads
            log_dir = os.path.dirname(self.kwargs['-1'])
            sample_name = self._trubase(self.kwargs['-1'])
            unal_key = '--un-conc'

        # parse the genome name
        genome_name = os.path.basename(self.kwargs['-x'])
        self.id = '_'.join(
            [sample_name, genome_name, self.timestamp, self.name])
        log_path = os.path.join(log_dir, self.id + '.log')

        # setup stdout redirect
        self.redirect = '2>&1 | tee -a {}'.format(log_path)

        # ensure unaligned reads are written to a file
        unal = os.path.splitext(self.kwargs['-S'])[0] + '.unal.fastq'
        self.kwargs.update({unal_key: unal})


class Bowtie2Cmd(HisatCmd):

    '''Bowtie 2 Aligner

    Current version uses the same parameters as hisat
    '''

    NAME = 'bowtie2'
    INVOKE_STR = 'bowtie2'
