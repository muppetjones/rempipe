
import time

from libpipe.cmd.dummy import CmdDummy
from libpipe.cmd.interface import CmdInterface
from libpipe.decorators.cmd import fall_through as fall_through_wrapper

import logging
log = logging.getLogger(__name__)


class PipeBase(CmdInterface):

    '''Pipe object for storing and handling a series of command-line objs

    PipeBase stores a series of command line objects and can execute those
    commands a number of fashions
        * Via a BASH script
        * Directly through python subprocesses
        * Through a resource manager system, such as Torque

    PipeBase also impements a full command-line interface:
        * cmd, link, output

    TODO(sjbush): Add Torque and GridEngine handling (both individually
        and through a pbs script).
    TODO(sjbush): Add direct execution via subprocesses (via cmd.run)
    TODO(sjbush): Allow commands to be skipped if their output is detected.
    TODO(sjbush): Implement pass_through (soft_output on rempipe) to
        collect unique output from all cmds. (as decorator)
    TODO(sjbush): Maintain a symbolic link to current (i.e., self) script
        and log files (named with timestamp).

    Arguments:
        fall_through: A boolean indicating whether input should be
            included in the output.

    Attributes:
        cmds: A list of commands in the pipeline

    Notable Methods:
        add: Accepts command objects, links and synchronizes them, and
            appends them to the list. Accepts any number of commands or
            commands given in a list. Chainable.

    '''

    def __init__(self, timestamp=None, fall_through=False, **kwargs):

        # save a timestamp, given or generated
        self.timestamp = (timestamp if timestamp
                          else time.strftime("%y%m%d-%H%M%S"))

        # check if 'input' was given
        # -- use kwargs b/c input is a builtin function
        try:
            _input = kwargs['input']
        except KeyError:
            _input = []

        self.dummy = CmdDummy(*_input)
        self.cmds = []

        if fall_through:
            self.output = fall_through_wrapper(self.output)

    #
    #   Cmd Interface
    #

    def cmd(self, cmd_sep='\n'):
        if not self.cmds:
            msg = 'Error: pipe is empty'
            raise ValueError(msg)
        return cmd_sep.join(self._get_cmd_strs())

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
            cmd.timestamp = self.timestamp
        self.cmds[-1].timestamp = self.timestamp

        return self

    def write(self, _file):

        with open(_file, 'w') as fh:
            fh.write(self.cmd())

    #
    #   Protected methods
    #

    def _get_cmd_strs(self):
        return [cmd.cmd() for cmd in self.cmds]
