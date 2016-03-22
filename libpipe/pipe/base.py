

from libpipe.cmd.dummy import CmdDummy


class PipeBase(object):

    '''Pipe object for storing and handling a series of command-line objs

    PipeBase stores a series of command line objects and can execute those
    commands a number of fashions
        * Via a BASH script
        * Directly through python subprocesses
        * Through a resource manager system, such as Torque

    TODO(sjbush): Add Torque and GridEngine handling (both individually
        and through a pbs script).
    TODO(sjbush): Add direct execution via subprocesses (via cmd.run)
    TODO(sjbush): Allow commands to be skipped if their output is detected.

    Attributes:


    '''

    def __init__(self, **kwargs):

        # check if 'input' was given
        # -- use kwargs b/c input is a builtin function
        try:
            _input = kwargs['input']
        except KeyError:
            _input = []

        self.dummy = CmdDummy(*_input)
        self.cmds = []

    def add(self, *cmds):
        '''Add a new cmd to the pipe and link it to existing commands'''

        self.cmds.extend(cmds)

        self.dummy.link(self.cmds[0])
        for cmd, next_cmd in zip(self.cmds[:-1], self.cmds[1:]):
            cmd.link(next_cmd)

        return self
