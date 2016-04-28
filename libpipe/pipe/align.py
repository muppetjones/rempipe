

from libpipe.cmd import align
from libpipe.cmd import count
from libpipe.cmd import samtools
from libpipe.pipe.base import Pipe
from libpipe.type import seq
from libpipe.util import path

import logging
log = logging.getLogger(__name__)


class AlignPipe(Pipe):

    def _setup(self):
        '''Setup a basic RNA-Seq pipeline from alignment to counting reads

        Creates and links commands for a simple RNA-Seq pipeline. Also
        attempts to define an input directory based on FASTQ file input.

        Pipeline:
            1. Hisat2  (cleaned reads, genome)
            2. Samtools sort (aligned SAM file)
            3. Samtools index (sorted BAM file)
            4. Bedtools multicov (sorted BAM file, genome)
        '''

        # genome passed from input
        _align = align.Hisat2Cmd(output_dir=self.output_dir)
        _sort = samtools.SamtoolsSortCmd()
        _index = samtools.SamtoolsIndexCmd()
        _count = count.HtseqCountCmd()

        cmd_list = [_align, _sort, _index, _count]
        self.add(*cmd_list)

        # set the input and output dir
        fastq_files = [f for f in self.input() if isinstance(f, seq.FastqType)]
        self.input_dir = path.common_directory(fastq_files)

        # HACK: Explicitly add index used by _align to _count args
        #   Better than some of the alternatives, at least until
        #   genome pass through is working
        # CRITICAL: No longer needed due to IndexType pass through!
        # self.cmds[0].cmd()
        # self.cmds[-1].args.append(self.cmds[0].kwargs['-x'])

    #
    #   Cmd Interface
    #

    def output(self):
        '''Add hisat2 unaligned fq to output'''
        _output = super().output()
        _unal = [o for o in self.cmds[0].output() if 'unal' in o]
        return _output + _unal
