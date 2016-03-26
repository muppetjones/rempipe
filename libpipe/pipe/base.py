
import os
import os.path
import stat
import subprocess
import time


from libpipe.cmd.dummy import CmdDummy
from libpipe.cmd.interface import CmdInterface
from libpipe.decorators.cmd import fall_through as fall_through_wrapper
from libpipe.decorators.io import file_or_handle
from libpipe import templates

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

    TODO(sjbush): Add 'job_name' attribute and use (+ timestamp) for files.
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

    def __init__(
            self,
            cmds=None, fall_through=False, template_path=None, timestamp=None,
            **kwargs):

        # save a timestamp, given or generated
        self.timestamp = (timestamp if timestamp
                          else time.strftime("%y%m%d-%H%M%S"))

        # check if 'input' was given
        # -- use kwargs b/c input is a builtin function
        try:
            _input = kwargs['input']
            del kwargs['input']
        except KeyError:
            _input = []

        self.dummy = CmdDummy(*_input)
        self.script_file = None
        self.cmds = []
        if cmds:
            self.add(*cmds)

        if fall_through:
            self.output = fall_through_wrapper(self.output)

        # initialize template information
        if template_path:
            self.pbs_template_path = template_path
        else:
            self.pbs_template_path = os.path.join(
                os.path.dirname(templates.__file__), 'template.pbs')
        self.pbs_template = None

        # pass all leftover kwargs to setup
        # -- rather, all that weren't deleted
        if not self.cmds:
            try:
                self._setup(**kwargs)
                kwargs = {}  # empty it! The child should handle cleanup.
            except NotImplementedError:
                pass  # no preset cmds...no worries

        # otherwise, complain about unexpected kwargs
        if kwargs:
            raise TypeError('__init__ got unexpected argument(s): {}'.format(
                list(kwargs.keys())))

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
    #   Overrideable
    #

    def _setup(self):
        '''Function used to create preset pipes--Implement in children'''
        raise NotImplementedError('Implement in children')

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

    def run(self, mode='local'):

        # we have a script (via `write`). use it.
        if self.script_file:
            log_file = os.path.splitext(self.script_file)[0] + '.log'
            self.log_file = log_file

            # TODO(sjbush): Make this dynamic using dir
            run_via = {
                'local': self._run_local,
                'pbs': self._run_pbs,
            }
            try:
                run_via[mode]()
            except KeyError:
                msg = 'Unknown run mode "{}". Try: {}'.format(
                    mode, list(run_via.keys()))
                raise ValueError(msg)

        # no script--call each cmd individually
        else:
            for cmd in self.cmds:
                cmd.run()

    @file_or_handle(mode='w')
    def write(self, fh, fmt=None):
        '''Handle write preparation and go to appropriate 'write' method

        Attributes:
            _file: The file to write to
            fmt: The script format. Must use if passing a handle.
                Default: 'pbs'.
        '''

        try:
            # get filename (if file-type handle)
            filename = fh.name
            extn = os.path.splitext(filename)[1]

            # pbs / shell specific processing
            if extn in ['.pbs', '.sh']:
                self._write_pbs_template(fh)

            self.script_file = filename
        except AttributeError:
            # Not a file handle (StringIO or otherwise)
            filename = None
        finally:
            self._write(fh)

        self._update_file_permissions(filename)  # nothing if filename = None
        return

    #
    #   Protected methods
    #

    def _get_cmd_strs(self):
        return [cmd.cmd() for cmd in self.cmds]

    def _load_pbs_template(self):
        if not self.pbs_template:
            with open(self.pbs_template_path, 'r') as ih:
                self.pbs_template = ih.read().rstrip()
        return

    def _run_local(self):
        try:
            subprocess.check_call(
                self.script_file, stdout=self.log_file, stderr=self.log_file)
        except subprocess.CalledProcessError:
            raise

    def _run_pbs(self):
        pass

    def _write(self, fh):
        fh.write(self.cmd())
        return

    def _write_pbs_template(self, fh):
        self._load_pbs_template()
        fh.write(self.pbs_template + "\n")
        return

    #
    #   Class Methods
    #

    @classmethod
    def _update_file_permissions(cls, _file):
        '''Make file executable by all parties

        http://stackoverflow.com/a/12792002/977342
        '''

        if _file:
            file_stat = os.stat(_file)
            st_all_exec = stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH
            os.chmod(_file, file_stat.st_mode | st_all_exec)
