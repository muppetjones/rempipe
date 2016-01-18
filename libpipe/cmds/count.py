import os.path
from libpipe.cmds.base import BaseCmd

import logging
log = logging.getLogger(__name__)


class HtseqCountCmd(BaseCmd):

    '''HT-Seq command line invocation

    Install:
        pip2 install numpy
        pip2 install scipy
        pip2 install matplotlib
        pip2 install htseq

    Versions:
        numpy 1.10.4
        scipy 0.16.1
        matplotlib 1.5.1
        htseq 0.6.1

    '''

    NAME = 'htseq-count'
    INVOKE_STR = 'htseq-count'

    ARGUMENTS = {
        (None, 'FILE', 'Alignment file (BAM or SAM)'),
        (None, 'FILE', 'Feature file (GFF or GTF)'),
        ('-f', 'SAMTYPE', 'Type of alignment file (sam or bam)'),
        ('-r', 'ORDER', 'Alignment sorting order (pos or name) [pe only]'),
        ('-s', 'STRANDED', 'Specify strand specific (yes, no, reverse)'),
        ('-a', 'INT', 'Minimum quality value (skip lower than)'),
        ('-t', 'CHAR', 'Feature type (3rd col) to use. Default: exon'),
        ('-i', 'CHAR', 'GFF attribute to use as feature id. Default: gene_id'),
        ('-m', 'MODE', 'Mode for feature overlap. Default: union. ' +
            'Choices: union, intersection-strict, intersection-nonempty'),
        ('-o', 'FILE', 'Write new SAM file with feature assignments'),
        ('-q', None, 'Suppress progress report'),
    }

    DEFAULTS = {
        '-r': 'pos',  # set for paired-end data (ignored for single-end)
        '-s': 'no',  # not strand specific
    }
    REQ_KWARGS = []
    REQ_ARGS = 2
    REQ_TYPE = [
        [(0, ), ('.bam', '.sam')],
        [(1, ), ('.gff', '.gtf', '.gff3')],
        [('-o',),  ('.sam', )],
    ]

    def __init__(
            self, *args,
            feature_type=None, feature_id=None, mode=None, **kwargs):
        super().__init__(*args, **kwargs)

        # Simpler interface for common functions
        if feature_type:
            self.kwargs['-t'] = feature_type
        if feature_id:
            self.kwargs['-i'] = feature_id
        if mode:
            self.kwargs['-m'] = mode

    def _prepreq(self):
        super()._prepreq()

        # CRITICAL: We probably have the args in the wrong order (sam first!)
        if self.args[1].endswith('am'):
            self.args.reverse()

        # set correct input format flag
        align_file = os.path.splitext(self.args[0])[1]
        self.kwargs['-f'] = align_file[1:].lower()  # remove leading dot

        # NOTE: From here down is duplicated in BedtoolsMulticovCmd
        #       EXCEPT for the line where the genome is added!!
        # TODO: make this a base function
        def _get_types(flag):
            for rt in self.REQ_TYPE:
                if flag in rt[0]:
                    return rt[1]
            return None

        # identify extension
        bed_file = self.args[1]
        bed_types = _get_types(1)
        try:
            bed_extn = [
                extn for extn in bed_types
                if bed_file.endswith(extn)
            ][0]
        except IndexError:
            # check for presence of file with expected extension
            # update with first found match
            # NOTE: LOGIC ERROR! If multiple exist, e.g., both .bed and .gff,
            #       only the first found will ever be used.
            for extn in bed_types:
                if os.path.isfile(bed_file + extn):
                    bed_extn = extn
                    # self.kwargs['-bed'] = bed_file + extn
                    break
            try:
                self.args[1] = bed_file + bed_extn
            except UnboundLocalError:
                msg = 'Unable to find an annotation file {} for "{}" genome'
                raise FileNotFoundError(msg.format(
                    bed_types, os.path.basename(bed_file)))

        # CRITICAL: GFF files probably will NOT have a gene_id attribute
        #           There's not a good, consistant substitute, except
        #           either 'gene' or 'Name'
        if not self.args[1].endswith('.gtf'):
            self.kwargs['-i'] = 'gene'

        self.redirect = ('>', os.path.splitext(self.args[0])[0] + '.count')

    def output(self):
        return (self.redirect[1])


class BedtoolsMulticovCmd(BaseCmd):

    NAME = 'bedtools_multicov'
    INVOKE_STR = 'bedtools multicov'

    ARGUMENTS = {
        ('-bams', 'FILE', 'The bam file'),
        ('-bed', 'FILE', 'The bed file'),
    }

    DEFAULTS = {}

    REQ_KWARGS = ['-bed', '-bams']
    REQ_ARGS = 0
    REQ_TYPE = [
        [('-bams', ), ('.bam', )],
        [('-bed', ), ('.bed', '.gff', '.gtf')],
    ]

    def _prepreq(self):
        super()._prepreq()

        # We may be given a genome prefix rather than an annotation file.
        # We need to (a) identify such cases and (b) add the appropriate
        # extension, defined by the presence of such a file.

        # TODO: make this a base function
        def _get_types(flag):
            for rt in self.REQ_TYPE:
                if flag in rt[0]:
                    return rt[1]
            return None

        # identify extension
        bed_file = self.kwargs['-bed']
        bed_types = _get_types('-bed')
        try:
            bed_extn = [
                extn for extn in bed_types
                if bed_file.endswith(extn)
            ][0]
        except IndexError:
            # check for presence of file with expected extension
            # update with first found match
            # NOTE: LOGIC ERROR! If multiple exist, e.g., both .bed and .gff,
            #       only the first found will ever be used.
            for extn in bed_types:
                if os.path.isfile(bed_file + extn):
                    bed_extn = extn
                    # self.kwargs['-bed'] = bed_file + extn
                    break
            try:
                self.kwargs['-bed'] = bed_file + bed_extn
            except UnboundLocalError:
                msg = 'Unable to find an annotation file {} for "{}" genome'
                raise FileNotFoundError(msg.format(
                    bed_types, os.path.basename(bed_file)))

    def _prepcmd(self):

        extn = os.path.splitext(self.kwargs['-bed'])[1]
        output = self.kwargs['-bams'].replace('.bam', '.count' + extn)
        self.redirect = ('>', output)

    def output(self):
        return (self.redirect[1],)
