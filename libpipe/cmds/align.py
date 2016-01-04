
import os.path
import time

from libpipe.cmds.base import BaseCmd

import logging
log = logging.getLogger(__name__)


class HisatCmd(BaseCmd):

    '''HISAT command setup

    Command usage:
        hisat [options]* -x <bt2-idx> {-1 <m1> -2 <m2> | -U <r>} -S <sam>

    '''

    NAME = 'hisat'
    INVOKE_STR = 'hisat'

    ARGUMENTS = [
        ('-x', 'FILE', 'Hisat reference genome index (base name)'),
        ('-1', 'FILE[,FILE]', 'comma separated list of paired-end 1 files'),
        ('-2', 'FILE[,FILE]', 'comma separated list of paired-end 2 files'),
        ('-U', 'FILE[,FILE]', 'comma separated list of unpaired reads'),
        ('-S', 'FILE', 'Output sam file (defaults to read prefix)'),
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

    REQ_KWARGS = ['-x', ('-1', '-2'), ['-1', '-U']]
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
        out_list = [self.kwargs['-S'], ]
        try:
            out_list.append(self.kwargs['--un'])
        except KeyError:
            out_list.append(self.kwargs['--un-conc'])
        return out_list

    #
    #   "Private" methods
    #

    def _prepcmd(self):
        '''Prep for hisat cmd

        > parse log file name and set for redirect
        > ensure unaligned reads are output
        '''
        # parse the genome name
        genome_name = os.path.basename(self.kwargs['-x'])

        # ensure we have an output file
        try:
            out_dir = os.path.dirname(self.kwargs['-S'])
            run_name = self._trubase(self.kwargs['-S'])
        except KeyError:
            try:
                # unpaired sequence file
                out_dir = os.path.dirname(self.kwargs['-U'])
                run_name = self._trubase(self.kwargs['-U'])
            except KeyError:
                # paired-end sequence file
                out_dir = os.path.dirname(self.kwargs['-1'])
                run_name = os.path.commonprefix(
                    self.kwargs['-1'], self.kwargs['-2'])

                # ensure common prefix includes some of the base name
                if run_name == out_dir:
                    run_name = self._trubase(self.kwargs['-1'])
                else:
                    run_name = os.path.basename(run_name)
            finally:
                # generated output name should contain genome name, too
                if genome_name not in run_name:
                    run_name = '_'.join([run_name, genome_name])

                # ensure we have '-S' set
                self.kwargs['-S'] = os.path.join(out_dir, run_name + '.sam')

        # set log file name
        self.id = '_'.join(
            [run_name, genome_name, self.timestamp, self.name])
        log_path = os.path.join(out_dir, self.id + '.log')

        # setup stdout redirect
        self.redirect = '2>&1 | tee -a {}'.format(log_path)

        # ensure unaligned reads are written to a file
        unal_key = '--un' if '-U' in self.kwargs else '--un-conc'
        unal = os.path.splitext(self.kwargs['-S'])[0] + '.unal.fastq'
        self.kwargs.update({unal_key: unal})


class Bowtie2Cmd(HisatCmd):

    '''Bowtie 2 Aligner

    Current version uses the same parameters as hisat
    '''

    NAME = 'bowtie2'
    INVOKE_STR = 'bowtie2'
