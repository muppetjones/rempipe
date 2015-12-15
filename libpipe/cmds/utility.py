
import os.path
from libpipe.cmds.base import BaseCmd


class FastqcCmd(BaseCmd):

    bin_name = 'fastqc'
    defaults = {}

    attributes = {
        'o': "output directory",
    }

    required_kwargs = []
    required_args = 1

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        print(args)
        # ensure the '-o' option is given
        # -- use the common prefix of the two input files
        # -- or the basename of the the first file otherwise
        if 'o' not in self.kwargs:
            print(self.args)
            prefix = os.path.commonprefix(self.args)
            if not prefix or prefix == self.args[0]:
                prefix = os.path.splitext(self.args[0])[0]
            self.kwargs['o'] = prefix

    def output(self):

        if len(self.args) > 1:
            out_list = [
                '{}-trimmed-pair{}.fastq'.format(self.kwargs['o'], i)
                for i, v in enumerate(self.args)
            ]
        else:
            out_list = ['{}-trimmed.fastq'.format(self.kwargs['o'])]
        return out_list


class SamtoolsSortCmd(BaseCmd):

    bin_name = 'samtools'
    subcmd = 'sort'
    defaults = {}

    attributes = {
        'o': 'Output file.'
        'T': None,
    }

    required_kwargs = []
    required_args = 1

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        # if we were only given a single file, make the prefix
        if 'o' not in self.kwargs:
            if self.args[0].endswith('sam'):
                bam_file = os.path.splitext(self.args[0])[0] + '.bam'
            else:
                bam_file = os.path.splitext(self.args[0])[0] + 's.bam'
            self.kwargs['o'] = bam_file

    def output(self):
        return self.kwargs['o']


class SamtoolsIndexCmd(BaseCmd):

    bin_name = 'samtools'
    subcmd = 'index'
    defaults = {}

    attributes = {}

    required_kwargs = []
    required_args = 1

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        # if we were only given a single file, make the prefix
        if not self.args[0].endswith('.bam'):
            raise ValueError('Samtools index requires a BAM file')

    def output(self):
        return self.args[0]  # return the given bam file
