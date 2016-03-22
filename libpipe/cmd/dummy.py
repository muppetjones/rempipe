

from libpipe.cmd import CmdAttributes
from libpipe.cmd import CmdBase
from libpipe.decorators.cmd import universal_fall_through


class CmdDummy(CmdBase):

    '''A dummy Cmd object for starting linking with Pipes

    The Cmd object was designed to automatically figure out which
    input values are matched with any given argument. Unfortunately,
    this is done via the 'input' attribute, which is only set when
    a command is linked to another. In other words, the first cmd
    of a chain must have the input set correctly (as the cmd expects).

    This class is meant to negate that disadvantage. Every Pipe object
    should have a dummy command that will pass through all arguments
    it is given. Linking the dummy to the first command will allow
    the expected automatic argument matching.

    TODO(sjbush): Make the dummy cmd() return various pipe details,
        commented out, of course.

    Notable Methods:
        input: A lambda method that mimics linkage to a previous command.
    '''

    attr = CmdAttributes(**{
        'invoke': None,
        'args': [],
    })

    def __init__(self, *args, **kwargs):

        input_args = list(args)
        input_kwargs = sorted(list(kwargs.values()))
        self.input = lambda: input_args + input_kwargs

    def cmd(self, *args, **kwargs):
        return ''

    @universal_fall_through
    def output(self):
        return []
