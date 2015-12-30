
import os.path
import libpipe.templates

from libpipe.cmds.base import BaseCmd


class SkewerCmd(BaseCmd):

    NAME = 'skewer'
    INVOKE_STR = 'skewer'

    DEFAULTS = {
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
    }

    attributes = {
        '-Q': "average quality threshold",
        '-q': "3' trim quality threshold",
        '-l': "minimum seq length",
        '-m': "method (paired end)",
        '-t': "processors",
        '-x': "adapters",
        '-o': "the prefix to use for output (will be <given>-trimmed-pair.fq)",
    }

    REQ_KWARGS = []
    REQ_ARGS = 1

    def __init__(self, *args, quiet=True, **kwargs):
        super().__init__(*args, **kwargs)

        # ensure the '-o' option is given
        # -- use the common prefix of the two input files
        # -- or the basename of the the first file otherwise
        if '-o' not in self.kwargs:
            prefix = os.path.commonprefix(self.args)
            if not prefix or prefix == self.args[0]:
                prefix = os.path.splitext(self.args[0])[0]
            self.kwargs['-o'] = prefix

        if quiet:
            self.flags.append('--quiet')

    @property
    def output(self):

        if len(self.args) > 1:
            out_list = [
                '{}-trimmed-pair{}.fastq'.format(self.kwargs['-o'], i + 1)
                for i, v in enumerate(self.args)
            ]
        else:
            out_list = ['{}-trimmed.fastq'.format(self.kwargs['-o'])]
        return out_list
