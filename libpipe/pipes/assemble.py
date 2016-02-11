import os

from libpipe.pipes.base import PresetPipe
from libpipe.pipes.qc import TrimPipe
from libpipe.cmds import (
    Hisat2Cmd,
    SamtoolsSortCmd, SamtoolsIndexCmd,
    HtseqCountCmd,
)

from libpipe.cmds.assemble import (
    VelvetkCmd, VelvethCmd, VelvetgCmd, AbacasCmd, ContigSummaryCmd,
)
from libpipe.cmds.utility import (
    InsertSizesCmd, ParseInsertSizesCmd,
)


import logging
log = logging.getLogger(__name__)


class AssemblePipe(PresetPipe):

    '''Assemble a bacterial genome

    NOTE: This pipe is intended for use as a subpipe.

    Velvet Pipeline
        1. Velvetk
        2. Velveth
        3. Velvetg
    Post processing
        4. Summarize contig sizes (BASH one-liner)
        5. abacas (order contigs according to reference genome)
        6. Summarize contig sizes (BASH one-liner)

    NOTE: We could move the post processing into another pipe, but
        both of these commands are fairly crucial to each assembly.
    '''

    REQ_PARAM = ['reference']

    def _setup(self, *pipe_args, input_list=[], reference='', **pipe_kwargs):

        # check for ins_len and ins_len_sd

        # 1. velvetk
        args = input_list[:]
        kwargs = {'--genome': reference}
        velvetk = VelvetkCmd(*args, wrap='k', **kwargs)

        # 2. velveth
        args = []
        kwargs = {'k': '${k}'}
        velveth = VelvethCmd(*args, **kwargs)

        # 3. velvetg
        args = []
        kwargs = {}
        kwargs.update({
            k: v for k, v in pipe_kwargs.items()
            if k in ['-ins_length', '-ins_length_sd']
        })
        velvetg = VelvetgCmd(*args, **kwargs)

        # 4. contig size summary
        args = []
        kwargs = {}
        contig_sizes = ContigSummaryCmd(*args, **kwargs)

        # 5. abacas
        args = []
        kwargs = {'-r': reference}
        abacas = AbacasCmd(*args, **kwargs)

        # 6. contig size summary
        args = []
        kwargs = {}
        contig_sizes = ContigSummaryCmd(*args, **kwargs)

        # add 'em to the pipe
        self.add(velvetk, velveth, velvetg, contig_sizes, abacas, )


class WgsPipe(PresetPipe):

    '''Assemble a bacterial genome, from sequencer to annotation

    Pipeline:
        1. TrimPipe
            a. FastQC (reads)
            b. Skewer (reads)
            c. FastQC (skewered reads)
        2. (optional) Filter reads through reference genome (usually human)
            a. Align to filter genome (will use unaligned output)
        3. Determine insertion length (skip)
            a. Remove unmapped reads
            b. Calculate insert length
        4. Run through Velvet
            a. estimate k (velvetk.pl)
            b. velveth
            c. velvetg
        5. Assembly statistics
            > MUST FIRST implement batch handling


    '''

    REQ_PARAM = ['input_list', 'odir', 'reference', ]

    def _setup(self, *pipe_args, **pipe_kwargs):
        # NOTE: TrimPipe uses 'input_list' to set FastQC input.

        # get cpu count
        # NOTE: this does NOT work with qsub -- needs update!
        cpu_count = os.cpu_count()

        # 1. TrimPipe
        trim = TrimPipe(*pipe_args, **pipe_kwargs)

        # 2. (optional) Filter through reference genome
        # > input: link from previous (fastq)
        # > output: sam files to given directory (+ the unaligned fastq)
        hisat_list = []
        try:
            for filter_genome in pipe_kwargs['filter']:
                args = []
                kwargs = {'-x': filter_genome, '-p': cpu_count}
                if kwargs['-p'] is None:
                    del kwargs['-p']
                else:
                    kwargs['-p'] = kwargs['-p'] - 1
                hisat = Hisat2Cmd(*args, **kwargs)
                hisat_list.append(hisat)
        except KeyError:
            pass  # no filter genomes -- no problem!

        # 2a -- sort sam file
        # > input: link from previous (sam)
        # > output: bam files to given directory
        args = []
        kwargs = {}
        st_sort = SamtoolsSortCmd(*args, **kwargs)

        # 2b -- index bam file
        # > input: link from previous (bam)
        # > output: bam index files to given directory
        args = []
        kwargs = {}
        st_index = SamtoolsIndexCmd(*args, **kwargs)

        # 3. Determine Insertion length

        # 3a. Calculate the insert length
        # > input: sam|bam file
        # > output: metrics file + prevous
        args = []
        kwargs = {
            # 'fall_through': True,
        }
        picard_ins_len = InsertSizesCmd(*args, **kwargs)

        # 3b. Parse the mean insert length
        # > input: metrics file
        # > output: previous
        args = ['mean_ins']
        kwargs = {
            'fall_through': True,
            'wrap': 'ins_len_mean'
        }
        mean_ins_len = ParseInsertSizesCmd(*args, **kwargs)

        # 3c. parse the insert length standard deviation
        # > input: metrics file
        # > output: previous
        args = ['stdev']
        kwargs = {
            'fall_through': True,
            'wrap': 'ins_len_stdev'
        }
        stdev_ins_len = ParseInsertSizesCmd(*args, **kwargs)

        # 4. Assemble
        # > input: link from previous (unaligned fastq)
        # > output: static files in child directory (contigs.fa)
        # -- remove input_list: don't accidently pass on untrimmed reads!
        del pipe_kwargs['input_list']
        pipe_kwargs.update({
            '-ins_length': '${ins_len_mean}',
            '-ins_length_sd': '${ins_len_stdev}',
        })
        assemble = AssemblePipe(
            *pipe_args,
            **pipe_kwargs
        )

        # add all valid commands
        sam_stuff = [st_sort, st_index]
        ins_len = [picard_ins_len, mean_ins_len, stdev_ins_len]
        cmd_list = filter(
            None, [trim] + hisat_list + sam_stuff + ins_len + [assemble])

        # add and link all commands
        self.add(*cmd_list)

        # circumvent the extra stuff
        # -- either the trimmed or filtered reads
        start_index = 0 + len(hisat_list)
        self.cmds[start_index].link(self.cmds[-1])
