

from libpipe.cmd import align
from libpipe.cmd import samtools
from libpipe.cmd import count
from libpipe.pipe.base import PipeBase


class AlignPipe(PipeBase):

    def _setup(self):

        hisat2 = align.Hisat2Cmd()
        st_sort = samtools.SamtoolsSortCmd()
        st_index = samtools.SamtoolsIndexCmd()
        htseq = count.HtseqCountCmd()

    #
    #   Cmd Interface
    #
