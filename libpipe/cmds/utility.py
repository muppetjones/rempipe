
import os.path
from libpipe.cmds.base import BaseCmd

import logging
log = logging.getLogger(__name__)


class FastqcCmd(BaseCmd):
    NAME = 'fastqc'
    INVOKE_STR = 'fastqc'

    DEFAULTS = {}

    attributes = {
        '-o': "output directory",
    }

    REQ_KWARGS = []
    REQ_ARGS = 1

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

        files = []
        for f in self.args:
            file_base = os.path.splitext(f)[0] + '_fastqc'
            files.extend([
                file_base + '.html',
                file_base + '.zip',
            ])

        # update the path if 'o' was given
        try:
            files = [
                os.path.join(self.kwargs['-o'], os.path.basename(f))
                for f in files
            ]
        except KeyError:
            pass

        return files


class SamtoolsSortCmd(BaseCmd):
    NAME = 'samtools_sort'
    INVOKE_STR = 'samtools sort'

    DEFAULTS = {}

    attributes = {
        '-o': 'Output file.',
        '-O': 'Output format.',
        '-T': 'Temp file prefix.',
    }

    REQ_KWARGS = []
    REQ_ARGS = 1

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        # if we were only given a single file, make the prefix
        if '-o' not in self.kwargs:
            if self.args[0].endswith('sam'):
                bam_file = os.path.splitext(self.args[0])[0] + '.bam'
            else:
                bam_file = os.path.splitext(self.args[0])[0] + 's.bam'
            self.kwargs['-o'] = bam_file

        if '-T' not in self.kwargs:
            t_kw = os.path.splitext(self.kwargs['-o'])[0] + '.tmp'
            self.kwargs['-T'] = t_kw

    @property
    def output(self):
        return [self.kwargs['-o'], ]


class SamtoolsIndexCmd(BaseCmd):
    NAME = 'samtools_index'
    INVOKE_STR = 'samtools index'

    DEFAULTS = {}

    attributes = {}

    REQ_KWARGS = []
    REQ_ARGS = 1

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        # if we were only given a single file, make the prefix
        if not self.args[0].endswith('.bam'):
            raise ValueError('Samtools index requires a BAM file')

    @property
    def output(self):
        return [self.args[0], ]  # return the given bam file
