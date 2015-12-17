
import os.path
from libpipe.cmds.base import BaseCmd


class SkewerCmd(BaseCmd):

    bin_name = 'skewer'
    defaults = {
        '-Q': 20,
        '-q': 9,
        '-l': 35,
        '-m': 'pe',
        '-t': 8,
        '-x': '/data/sbush/data/adapters/adapters_illumina_truseq.fa',
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

    required_kwargs = []
    required_args = 1

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        # ensure the '-o' option is given
        # -- use the common prefix of the two input files
        # -- or the basename of the the first file otherwise
        if '-o' not in self.kwargs:
            prefix = os.path.commonprefix(self.args)
            if not prefix or prefix == self.args[0]:
                prefix = os.path.splitext(self.args[0])[0]
            self.kwargs['-o'] = prefix

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
