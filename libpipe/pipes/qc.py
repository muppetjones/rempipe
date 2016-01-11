import os

from libpipe.pipes.base import BasePipe
from libpipe.cmds import (
    SkewerCmd, FastqcCmd,
)

import logging
log = logging.getLogger(__name__)


class TrimPipe(BasePipe):

    '''Do QC check, trim, then recheck QC

    Pipeline:
        1: FastQC
        2: Skewer
        3: FastQC
    '''

    def __init__(self, *args, input_list=[], genome='', odir='', **kwargs):
        super().__init__(*args, **kwargs)

        try:
            self._setup(input_list=input_list, genome=genome, odir=odir)
        except TypeError:
            raise

    def _setup(self, input_list, genome, odir):

        # Step 1 -- assess quality
        # > input: given files
        # > output: to given directory
        args = input_list
        kwargs = {'-o': odir}
        fastqc_raw = FastqcCmd(*args, **kwargs)

        # Step 2 -- trim reads
        # > input: link from previous (fastq)
        # > output: trimmed read files
        args = []
        kwargs = {'-o': os.path.join(odir, self.job_name)}
        skewer = SkewerCmd(*args, **kwargs)

        # Step 3 -- assess quality after trim
        # > input: link from previous (fastq)
        # > output: to given directory
        args = []
        kwargs = {'-o': odir}
        fastqc_trim = FastqcCmd(*args, **kwargs)

        self.add(
            fastqc_raw, skewer, fastqc_trim,
        )

        return
