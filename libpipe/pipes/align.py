
import os

from libpipe.pipes.base import BasePipe, PresetPipe
from libpipe.pipes.qc import TrimPipe
from libpipe.cmds import (
    Hisat2Cmd,
    SamtoolsSortCmd, SamtoolsIndexCmd, BedtoolsMulticovCmd,
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

    def _setup(self, input_list=[], genome='', odir=''):
        if not isinstance(genome, str):
            msg = '{} only accepts a single genome. {}'.format(
                self.__class__.__name__,
                'Use "NestedGenomicsPipe" for multiple'
            )
            raise TypeError(msg)

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

        # Step 4 -- count index
        # > input: link from previous (bam)
        # > output: bam index files to given directory
        args = []
        kwargs = {'-bed': genome}
        bt_multicov = BedtoolsMulticovCmd(*args, **kwargs)

        self.add(
            hisat,
            st_sort, st_index,
            bt_multicov,
        )

        return


class AlignPipe(BasePipe):

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

    def __init__(self, *args, input_list=[], genome=[], **kwargs):

        # This expects multiple genomes, but Genomics Pipe expects a single
        # genome
        if not genome or not isinstance(genome, list):
            raise TypeError('{} requires more than one genome'.format(
                self.__class__.__name__))
        self.secondary_genome = genome[1:]
        genome = genome[0]

        super().__init__(
            *args, input_list=input_list, genome=genome, **kwargs)
