
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


class AlignPipe(PresetPipe):

    '''Execute a basic RNA-Seq pipeline from alignment to counting reads

    NOTE: This pipe is intended for use as a subpipe.

    Pipeline:
        1. Hisat  (skewered reads, genome)
        2. Samtools sort (aligned SAM file)
        3. Samtools index (sorted BAM file)
        4. Bedtools multicov (sorted BAM file, genome)
    '''

    def _setup(self, input_list=[], reference='', odir='', **kwargs):
        log.debug(reference)

        if not isinstance(reference, str):

            # could be a size 1 list!
            if len(reference) != 1:
                msg = '{} only accepts a single genome. {}'.format(
                    self.__class__.__name__,
                    'Use "NestedGenomicsPipe" for multiple'
                )
                raise ValueError(msg)
            else:
                reference = reference[0]

        # Step 1 -- align trimmed reads
        # > input: link from previous (fastq)
        # > output: sam files to given directory
        args = []
        kwargs = {'-x': reference, '-p': os.cpu_count() - 1}
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
        args = [reference, ]
        kwargs = {}
        count_cmd = HtseqCountCmd(*args, **kwargs)

        self.add(
            hisat,
            st_sort, st_index,
            count_cmd,
        )

        return


class RnaSeqPipe(PresetPipe):

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

    REQ_PARAM = ['input_list', 'reference', 'odir']

    def _setup(self, *args, **kwargs):
        # NOTE: TrimPipe uses 'input_list' to set FastQC input.
        log.debug(args)
        log.debug(kwargs)

        trim = TrimPipe(*args, **kwargs)
        align = AlignPipe(*args, **kwargs)

        self.add(trim, align)


class NestedRnaSeqPipe(RnaSeqPipe):

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

    REQ_PARAM = ['input_list', 'reference', 'odir']

    def _setup(self, *args, reference=[], **kwargs):

        try:
            _filter = kwargs['filter']
        except KeyError:
            _filter = None  # reference[0]
            # reference = reference[1:]

        # or not isinstance(reference, list) or len(reference) < 2:
        if not _filter:
            msg = 'Multiple genomes must be given. For a single genome,' + \
                'use "RnaSeqPipe"'
            raise ValueError(msg)

        # kwargs['genome'] = kwargs['filter']
        # secondary_genomes = reference[:]

        # setup initial trim and alignment
        # CRITICAL: Must use soft_output to ensure the unaligned fastq
        #           files from hisat are available as input
        _kwargs_filter = dict.copy(kwargs)
        _kwargs_filter['reference'] = _filter
        super()._setup(*args, **_kwargs_filter)
        self.cmds[-1].soft_output = True

        # create subpipe for each genome
        # CRITICAL: remove the genome from the kwargs before comprehension!!
        # CRITICAL: each subpipe MUST have soft_output, same as above
        # del kwargs['genome']
        secondary_alignments = [
            AlignPipe(*args, soft_output=True, reference=ref, **kwargs)
            for ref in reference
        ]

        self.add(*secondary_alignments)
