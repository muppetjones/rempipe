

from libpipe.cmd.dummy import CmdDummy
from libpipe.cmd.interface import CmdInterface
from libpipe.decorators.cmd import fall_through

import logging
log = logging.getLogger(__name__)


class PipeBase(CmdInterface):

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

        try:
            if kwargs['fall_through']:
                self.output = fall_through(self.output)
        except KeyError:
            pass  # do nothing--output is fine as is

    #
    #   Cmd Interface
    #

    def cmd(self):
        pass

    def link(self, cmd):
        if not self.cmds:
            msg = 'Error: pipe is empty'
            raise ValueError(msg)
        cmd.input = self.output
        return cmd

    @property
    def input(self):
        '''Pass through input to dummy cmd

        Possibly make a property? interface with dummy cmd.
        '''
        return self.dummy.input

    @input.setter
    def input(self, cmd_output_method):
        self.dummy.input = cmd_output_method
        return

    def output(self):
        try:
            return self.cmds[-1].output()
        except IndexError:
            msg = 'Error: pipe is empty'
            raise ValueError(msg)

    #
    #   Public methods
    #

    def add(self, *cmds):
        '''Add a new cmd to the pipe and link it to existing commands

        NOTE: Because it's so easy to forget the '*',
            'add' will also accept a single list argument.

        Arguments:
            *cmds: Any number of cmds
        '''

        # ID10T proof -- just in case you forget the '*'
        if len(cmds) == 1 and isinstance(cmds[0], list):
            cmds = cmds[0]

        self.cmds.extend(cmds)
        self.dummy.link(self.cmds[0])
        for cmd, next_cmd in zip(self.cmds[:-1], self.cmds[1:]):
            cmd.link(next_cmd)

        return self
