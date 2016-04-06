
import os.path

from libpipe.cmd import align
from libpipe.cmd import attr
from libpipe.cmd import base

import logging
log = logging.getLogger(__name__)


class HtseqCountCmd(base.CmdBase):

    '''HT-Seq command line invocation

    TODO(sjbush): Add venv handling.
    TODO(sjbush): Check for gtf or gff (instead of just adding)
    TODO(sjbush): Add genome-specific config file to set -t and -i.
    TODO(sjbush): Add redirect to log file.

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

    attr = attr.CmdAttributes(
        name='htseq-count',
        invoke='htseq-count',

        args=[
            (0, 'FILE', 'Alignment file (BAM or SAM)'),
            (1, 'FILE', 'Feature file (GFF or GTF)'),
            ('-f', 'SAMTYPE', 'Type of alignment file (sam or bam)'),
            ('-r', 'ORDER', 'Alignment sorting order (pos or name) [pe only]'),
            ('-s', 'STRANDED', 'Specify strand specific (yes, no, reverse)'),
            ('-a', 'INT', 'Minimum quality value (skip lower than)'),
            ('-t', 'CHAR', 'Feature type (3rd col) to use. Default: exon'),
            ('-i', 'CHAR',
             'GFF attribute to use as feature id. Default: gene_id'),
            ('-m', 'MODE', 'Mode for feature overlap. Default: union. ' +
             'Choices: union, intersection-strict, intersection-nonempty'),
            ('-o', 'FILE', 'Write new SAM file with feature assignments'),
            ('-q', None, 'Suppress progress report'),
        ],

        defaults={
            '-r': 'pos',  # set for paired-end data (ignored for single-end)
            '-s': 'no',  # not strand specific
        },
        req_kwargs=[],
        req_args=2,
        req_types=[
            [(0, ), ('.bam', '.sam')],
            [(1, ), (align.Hisat2Index, '.gff', '.gtf', '.gff3')],
            [('-o',), ('.sam', )],
        ],
    )

    def __init__(self, *args, genome=None, **kwargs):
        super().__init__(*args, **kwargs)

        if genome:
            self.args.append(genome)

    #
    #   Cmd interface
    #

    def output(self):
        return [self.redirect[1], ]

    #
    #   Overrides
    #

    def _pre_req(self):

        # swap args if the first pos arg is NOT a bam or sam
        # -- naive: does not check anything else, including no. args!
        prefix, extn = os.path.splitext(self.args[0])
        if extn not in ('.bam', '.sam'):
            self.args.reverse()

        # set the format based on the alignment file extension
        # -- even if not BAM or SAM, will catch during check requirements
        # -- re-check extn (may have been swapped!)
        self.kwargs['-f'] = extn[1:].lower()  # remove leading dot

        # add a gtf extension to the second pos arg if missing
        # -- allows init(genome) to be given a bowtie index
        # -- naive: does not check for existence of file
        try:
            extn = os.path.splitext(self.args[1])[1]
            if not extn:
                self.args[1] = self.args[1] + '.gtf'
        except IndexError:
            # CRITICAL: We *should* have 2 args, but we're not checking
            #   here. We'll catch it with 'req_args' later.
            pass

    def _pre_cmd(self):

        # set redirect
        # -- possibly allow user to set in future implementation
        count_file = '{}.count'.format(os.path.splitext(self.args[0])[0])
        self.redirect = ('>', count_file)
