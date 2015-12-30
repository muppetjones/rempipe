import os.path
from libpipe.cmds.base import BaseCmd

import logging
log = logging.getLogger(__name__)


class BedtoolsMulticovCmd(BaseCmd):

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
        log.debug(self.args)

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
        return (self.kwargs['-o'],)
