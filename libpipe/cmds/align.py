
import os.path
import time

from libpipe.cmds.base import BaseCmd


class HisatCmd(BaseCmd):

    NAME = 'hisat'
    INVOKE_STR = 'hisat'

    HELP_LIST = [
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
        '-p': "$(wc -l < $PBS_NODEFILE)",
        '-I': 0,
        '-X': 500,
    }

    REQ_KWARGS = ['-x', ('-1', '-2'), '-S']
    REQ_ARGS = 0

    def __init__(
            self, *args,
            encoding='--phred33', orientation='--fr', format='-q',
            timestamp=None,
            **kwargs):

        try:
            super().__init__(*args, **kwargs)
        except ValueError:
            raise   # requirements failure; pass it on up

        self.flags.extend([
            encoding, orientation
        ])

        # parse log file name (based on given file info)
        try:
            log_dir = os.path.dirname(self.kwargs['-U'])
            sample_name = self._trubase(self.kwargs['-U'])
        except AttributeError:
            log_dir = os.path.dirname(self.kwargs['-1'])
            sample_name = self._trubase(self.kwargs['-1'])

        # sample_name: the basename of the given file
        genome_name = os.path.basename(self.kwargs['-x'])
        timestamp = timestamp if timestamp else time.strftime("%y%m%d-%H%M%S")

        log_file = '_'.join(
            [sample_name, genome_name, timestamp, self.name]) + '.log'
        log_path = os.path.join(log_dir, log_file)

        # setup stdout redirect
        self.redirect = '2>&1 | tee -a {}'.format(log_path)

        # ensure unaligned reads are written to a file
        un_conc = os.path.splitext(self.kwargs['-S'])[0] + '.unal.fastq'
        self.kwargs.update({'--un-conc': un_conc})

    @property
    def output(self):
        return [self.kwargs['-S']]


class Bowtie2Cmd(HisatCmd):

    '''Bowtie 2 Aligner

    Current version uses the same parameters as hisat
    '''

    NAME = 'bowtie2'
    INVOKE_STR = 'bowtie2'
