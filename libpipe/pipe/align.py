

from libpipe.cmd import align
from libpipe.cmd import count
from libpipe.cmd import samtools
from libpipe.pipe.base import PipeBase

import logging
log = logging.getLogger(__name__)


class AlignPipe(PipeBase):

    def _setup(self):

        _align = align.Hisat2Cmd()  # genome passed from input
        _sort = samtools.SamtoolsSortCmd()
        _index = samtools.SamtoolsIndexCmd()
        _count = count.HtseqCountCmd()

        cmd_list = [_align, _sort, _index, _count]
        self.add(*cmd_list)

        # HACK: Explicitly add index used by _align to _count args
        #   Better than some of the alternatives, at least until
        #   genome pass through is working
        self.cmds[0].cmd()
        self.cmds[-1].args.append(self.cmds[0].kwargs['-x'])

    #
    #   Cmd Interface
    #

    def output(self):
        '''Add hisat2 unaligned fq to output'''
        _output = super().output()
        _unal = [o for o in self.cmds[0].output() if 'unal' in o]
        return _output + _unal
