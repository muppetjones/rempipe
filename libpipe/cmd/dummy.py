

from libpipe.cmd import CmdAttributes
from libpipe.cmd import CmdBase


class CmdDummy(CmdBase):

    attr = CmdAttributes(**{
        'invoke': None,
        'args': [],
    })

    def output(self):
        pass
