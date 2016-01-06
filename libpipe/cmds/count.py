import os.path
from libpipe.cmds.base import BaseCmd

import logging
log = logging.getLogger(__name__)


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
        self.redirect = '> {}'.format(output)
        self.count_file = output

    def output(self):
        return (self.count_file,)
