
import os

from libpipe.pipes.base import PresetPipe
from libpipe.pipes.qc import TrimPipe
from libpipe.cmds import (
    Hisat2Cmd,
    SamtoolsSortCmd, SamtoolsIndexCmd,
    HtseqCountCmd,
)


import logging
log = logging.getLogger(__name__)


class _AlignPipe(PresetPipe):

    '''Execute a basic RNA-Seq pipeline from alignment to counting reads

    NOTE: This pipe is intended for use as a subpipe.

    Pipeline:
        1. Hisat  (skewered reads, genome)
        2. Samtools sort (aligned SAM file)
        3. Samtools index (sorted BAM file)
        4. Bedtools multicov (sorted BAM file, genome)
    '''

    def _setup(self, input_list=[], genome='', odir='', **kwargs):
        if not isinstance(genome, str):

            # could be a size 1 list!
            if len(genome) != 1:
                msg = '{} only accepts a single genome. {}'.format(
                    self.__class__.__name__,
                    'Use "NestedGenomicsPipe" for multiple'
                )
                raise ValueError(msg)
            else:
                genome = genome[0]

        # Step 1 -- align trimmed reads
        # > input: link from previous (fastq)
        # > output: sam files to given directory
        args = []
        kwargs = {'-x': genome, '-p': os.cpu_count() - 1}
        if kwargs['-p'] is None:
            del kwargs['-p']
        hisat = Hisat2Cmd(*args, **kwargs)

        # Standalone usage (vs. linked)
        if input_list:
            hisat.input = lambda: list(input_list)

        # Step 2 -- sort sam file
        # > input: link from previous (sam)
        # > output: bam files to given directory
        args = []
        kwargs = {}
        st_sort = SamtoolsSortCmd(*args, **kwargs)

        # Step 3 -- index bam file
        # > input: link from previous (bam)
        # > output: bam index files to given directory
        args = []
        kwargs = {}
        st_index = SamtoolsIndexCmd(*args, **kwargs)

        # # Step 4 -- count index [DEPRECATED -- USE HtseqCountCmd]
        # # > input: link from previous (bam)
        # # > output: bam index files to given directory
        # args = []
        # kwargs = {'-bed': genome}
        # bt_multicov = BedtoolsMulticovCmd(*args, **kwargs)

        # Step 4 -- count index
        # > input: link from previous (bam)
        # > output: bam index files to given directory
        args = [genome, ]
        kwargs = {}
        count_cmd = HtseqCountCmd(*args, **kwargs)

        self.add(
            hisat,
            st_sort, st_index,
            count_cmd,
        )

        return


class AlignPipe(PresetPipe):

    '''Execute a basic RNA-Seq pipeline from sequencing to counting reads

    Pipeline:
        1. FastQC (reads)
        2. Skewer (reads)
        3. FastQC (skewered reads)
        4. Hisat  (skewered reads, genome)
        5. Samtools sort (aligned SAM file)
        6. Samtools index (sorted BAM file)
        7. Bedtools multicov (sorted BAM file, genome)
    '''

    REQ_PARAM = ['input_list', 'genome', 'odir']

    def _setup(self, *args, **kwargs):
        # NOTE: TrimPipe uses 'input_list' to set FastQC input.
        trim = TrimPipe(*args, **kwargs)
        align = _AlignPipe(*args, **kwargs)

        self.add(trim, align)


class NestedAlignPipe(AlignPipe):

    '''Execute a basic RNA-Seq genomics pipeline with multiple genomes

    Basic genomics pipeline, but with multiple subsequent genome alignments.
    For each genome, use the unaligned output from the previous genome.
    NOTE: Genomes MUST already be in desired order.

    Pipeline:
        1. FastQC (reads)
        2. Skewer (reads)
        3. FastQC (skewered reads)
        4. Hisat  (skewered reads, genome)
        5. Samtools sort (aligned SAM file)
        6. Samtools index (sorted BAM file)
        7. Bedtools multicov (sorted BAM file, genome)

        Then, repeat 4-7 using separate genomes and the unaligned output
        from hisat (step 4).
    '''

    def _setup(self, *args, genome=[], **kwargs):

        if not isinstance(genome, list) or len(genome) < 2:
            msg = 'Multiple genomes must be given. For a single genome,' + \
                'use "AlignPipe"'
            raise ValueError(msg)

        kwargs['genome'] = genome[0]
        secondary_genomes = genome[1:]

        # setup initial trim and alignment
        # CRITICAL: Must use soft_output to ensure the unaligned fastq
        #           files from hisat are available as input
        super()._setup(*args, **kwargs)
        self.cmds[-1].soft_output = True

        # create subpipe for each genome
        # CRITICAL: remove the genome from the kwargs before comprehension!!
        # CRITICAL: each subpipe MUST have soft_output, same as above
        del kwargs['genome']
        secondary_alignments = [
            _AlignPipe(*args, soft_output=True, genome=gen, **kwargs)
            for gen in secondary_genomes
        ]

        self.add(*secondary_alignments)
