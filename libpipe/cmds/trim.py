
import os.path
import libpipe.templates

from libpipe.cmds.base import BaseCmd, CmdAttributes

import logging
log = logging.getLogger(__name__)


class SkewerCmd(BaseCmd):

    '''Skewer Command

    Command usage:
        skewer [options] <reads.fastq> [paired-reads.fastq]
    '''

    attr = CmdAttributes(
        name='skewer',
        invoke_str='skewer',

        arguments=[
            ('-Q', 'INT',
             'The lowest mean quality value allowed before trimming'),
            ('-q', 'INT', 'Trim 3\' end until INT or higher quality reached'),
            ('-l', 'INT', 'The minimum read length allowed after trimming'),
            ('-m', 'CHAR', 'Trimming mode; \n\t' +
             '1) single-end [head: 5\' end, tail: 3\' end, any: anywhere]\n\t' +
             '2) paired-end [pe: paired-end, mp: mate-pair, ap: amplicon'),
            ('-t', 'INT', 'Number of concurrent threads (processors) [1, 32]'),
            ('-x', 'FILE', 'Adapter file (or sequence)'),
            ('-o', 'PREFIX',
             'Prefix to use for output (PREFIX-trimmed-pair.fq)'),
        ],
        defaults={
            '-Q': 20,
            '-q': 9,
            '-l': 35,
            '-m': 'pe',
            '-t': 8,
            '-x': os.path.join(
                os.path.dirname(
                    libpipe.templates.__file__),
                'adapters_illumina_pe.fa',
            ),
        },

        req_kwargs=[],
        req_args=1,
        req_type=[
            [(0, 1), ('.fastq', '.fq')],
        ],
    )

    def __init__(self, *args, quiet=True, **kwargs):
        super().__init__(*args, **kwargs)

        if quiet:
            self.flags.append('--quiet')

    def _prepcmd(self):
        # ensure the '-o' option is given
        # -- use the common prefix of the two input files
        # -- or the basename of the the first file otherwise
        if '-o' not in self.kwargs:
            prefix = os.path.commonprefix(self.args)
            if not prefix or prefix == self.args[0]:
                prefix = os.path.splitext(self.args[0])[0]
            self.kwargs['-o'] = prefix

        if len(self.args) == 1:
            self.kwargs['-m'] = 'tail'

    def output(self):

        if len(self.args) > 1:
            out_list = [
                '{}-trimmed-pair{}.fastq'.format(self.kwargs['-o'], i + 1)
                for i, v in enumerate(self.args)
            ]
        else:
            # Note that there is NOT a number after trimmed
            out_list = ['{}-trimmed.fastq'.format(self.kwargs['-o'])]
        return out_list
