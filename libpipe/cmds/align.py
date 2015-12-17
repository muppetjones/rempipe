
import os.path
import time

from libpipe.cmds.base import BaseCmd


class HisatCmd(BaseCmd):

    bin_name = 'hisat'

    attributes = {
        '-x': 'FILE\tHisat reference genome index (base name)',
        '-1': 'FILE[,FILE]\tcomma separated list of paired-end mate 1 files',
        '-2': 'FILE[,FILE]\tcomma separated list of paired-end mate 2 files',
        '-U': 'FILE[,FILE]\tcomma separated list of unpaired reads',
        '-S': 'FILE\toutput sam file',
        '-p': 'INT\tnumber of processors',
        '-I': 'INT\tminimum fragment length. Default = 0.',
        '-X': 'INT\taximum fragment length. Default = 500.',
        '--phred33': 'Illumina 1.9+ encoding',
        '--phred64': 'Illumina 1.8 and earlier encoding',
        '--fr': 'Upstream downstream mate orientations',
        '--un-conc': 'PATH\tPath to write unaligned, paired-end reads to.',
        '-q': 'Reads are FASTQ files',
        '-f': 'Reads are FASTA files',
    }
    defaults = {
        '-p': "$(wc -l < $PBS_NODEFILE)",
        '-I': 0,
        '-X': 500,
    }

    required_kwargs = ['-x', ('-1', '-2'), '-S']
    required_args = 0

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
            log_dir = os.path.dirname(self.kwargs['-S'])
        except AttributeError:
            log_dir = os.path.dirname(self.kwargs['-1'])

        log_sample = os.path.splitext(os.path.basename(self.kwargs['-1']))[0]
        log_genome = os.path.basename(self.kwargs['-x'])
        log_time = timestamp if timestamp else time.strftime("%y%m%d-%H%M%S")

        log_file = '_'.join(
            [self.bin_name, log_sample, log_genome, log_time]) + '.log'
        log_path = os.path.join(log_dir, log_file)

        # setup stdout redirect
        self.redirect = '2>&1 | tee -a {}'.format(log_path)

        # ensure unaligned reads are written to a file
        un_conc = os.path.splitext(self.kwargs['-S'])[0] + '.unal.fastq'
        self.kwargs['--un-conc'] = un_conc

    @property
    def output(self):
        return [self.kwargs['-S']]


class Bowtie2Cmd(HisatCmd):

    '''Bowtie 2 Aligner

    Current version uses the same parameters as hisat
    '''

    bin_name = 'bowtie2'
