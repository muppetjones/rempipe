
import os

from libpipe.pipes.base import BasePipe
from libpipe.cmds import (
    SkewerCmd, Hisat2Cmd, FastqcCmd,
    SamtoolsSortCmd, SamtoolsIndexCmd, BedtoolsMulticovCmd,
)

import logging
log = logging.getLogger(__name__)


class GenomicsPipe(BasePipe):

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

    def __init__(self, *args, input_list=[], genome='', odir='', **kwargs):
        super().__init__(*args, **kwargs)

        try:
            self._setup(input_list=input_list, genome=genome, odir=odir)
        except TypeError:
            raise

    def _setup(self, input_list, genome, odir):
        if not isinstance(genome, str):
            msg = '{} only accepts a single genome. {}'.format(
                self.__class__.__name__,
                'Use "NestedGenomicsPipe" for multiple'
            )
            raise TypeError(msg)

        # Step 1 -- assess quality
        # > input: given files
        # > output: to given directory
        args = input_list
        kwargs = {'-o': odir}
        fastqc_raw = FastqcCmd(*args, **kwargs)

        # Step 2 -- trim reads
        # > input: link from previous (fastq)
        # > output: trimmed read files
        args = []
        kwargs = {'-o': os.path.join(odir, self.job_name)}
        skewer = SkewerCmd(*args, **kwargs)

        # Step 3 -- assess quality after trim
        # > input: link from previous (fastq)
        # > output: to given directory
        args = []
        kwargs = {'-o': odir}
        fastqc_trim = FastqcCmd(*args, **kwargs)

        # Step 4 -- align trimmed reads
        # > input: link from previous (fastq)
        # > output: sam files to given directory
        args = []
        kwargs = {'-x': genome, '-p': os.cpu_count() - 1}
        if kwargs['-p'] is None:
            del kwargs['-p']
        hisat = Hisat2Cmd(*args, **kwargs)

        # Step 5 -- sort sam file
        # > input: link from previous (sam)
        # > output: bam files to given directory
        args = []
        kwargs = {}
        st_sort = SamtoolsSortCmd(*args, **kwargs)

        # Step 6 -- index bam file
        # > input: link from previous (bam)
        # > output: bam index files to given directory
        args = []
        kwargs = {}
        st_index = SamtoolsIndexCmd(*args, **kwargs)

        # Step 7 -- count index
        # > input: link from previous (bam)
        # > output: bam index files to given directory
        args = []
        kwargs = {'-bed': genome}
        bt_multicov = BedtoolsMulticovCmd(*args, **kwargs)

        self.add(
            fastqc_raw, skewer, fastqc_trim,
            hisat, st_sort, st_index,
            bt_multicov,
        )

        return


class NestedGenomicsPipe(GenomicsPipe):

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
