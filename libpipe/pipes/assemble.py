import os

from libpipe.pipes.base import PresetPipe
from libpipe.pipes.qc import TrimPipe
from libpipe.cmds import (
    Hisat2Cmd,
    SamtoolsSortCmd, SamtoolsIndexCmd,
    HtseqCountCmd,
)

from libpipe.cmds.assemble import (
    VelvetkCmd, VelvethCmd, VelvetgCmd, ContigSummaryCmd
)

from libpipe.cmds.utility import (
    AbacasCmd,
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

    NOTE: We could move the post processing into another pipe, but
        both of these commands are fairly crucial to each assembly.
    '''

    REQ_PARAM = ['input_list', 'reference']

    def _setup(self, *pipe_args, input_list=[], reference='', **pipe_kwargs):

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
        velvetg = VelvetgCmd(*args, **kwargs)

        # 4. contig size summary
        args = []
        kwargs = {}
        contig_sizes = ContigSummaryCmd(*args, **kwargs)

        # 5. abacas
        args = []
        kwargs = {'-r': reference}
        abacas = AbacasCmd(*args, **kwargs)

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

        # check for filter genome input

        # 1. TrimPipe
        trim = TrimPipe(*pipe_args, **pipe_kwargs)

        # 2. (optional) Filter through reference genome
        # > input: link from previous (fastq)
        # > output: sam files to given directory (+ the unaligned fastq)
        if 'filter' not in pipe_kwargs:
            hisat = None
        else:
            args = []
            kwargs = {'-x': pipe_kwargs['filter'], '-p': os.cpu_count()}
            if kwargs['-p'] is None:
                del kwargs['-p']
            else:
                kwargs['-p'] = kwargs['-p'] - 1
            hisat = Hisat2Cmd(*args, **kwargs)

        # 3. Determine Insertion length
        # ** SKIP **

        # 4. Assemble
        # > input: link from previous (unaligned fastq)
        # > output: static files in child directory (contigs.fa)
        assemble = AssemblePipe(*pipe_args, **pipe_kwargs)

        # add all valid commands
        cmd_list = filter(
            None, [trim, hisat, assemble])
        self.add(*cmd_list)
