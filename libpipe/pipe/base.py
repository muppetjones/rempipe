
import datetime
import os
import os.path
import re
import stat
import subprocess
import sys
import time


from libpipe.cmd.dummy import CmdDummy
from libpipe.cmd.interface import CmdInterface
from libpipe.decorators.cmd import fall_through as fall_through_wrapper
from libpipe.decorators.io import file_or_handle
from libpipe import templates

import logging
log = logging.getLogger(__name__)


class Pipe(CmdInterface):

    '''Pipe object for storing and handling a series of command-line objs

    Pipe stores a series of command line objects and can execute those
    commands a number of fashions
        * Via a BASH script
        * Directly through python subprocesses
        * Through a resource manager system, such as Torque

    Pipe also impements a full command-line interface:
        * cmd, link, output

    Arguments/Attributes:
        fall_through: A boolean indicating whether input should be
            included in the output.
        cmds: A list of commands in the pipeline
        odir/output_dir: The dir for output. In most cases, the Cmd
            should handle this itself.
        template_path: The path to the script template.
        timestamp: A timestamp string.

    Notable Methods:
        add: Accepts command objects, links and synchronizes them, and
            appends them to the list. Accepts any number of commands or
            commands given in a list. Chainable.

    '''
    # TODO(sjbush): Add 'job_name' attribute and use (+ timestamp) for files.
    # TODO(sjbush): Add Torque and GridEngine handling (both individually
    #   and through a pbs script).
    # TODO(sjbush): Add direct execution via subprocesses (via cmd.run)
    # TODO(sjbush): Allow commands to be skipped if their output is detected.
    # TODO(sjbush): Implement pass_through (soft_output on rempipe) to
    #   collect unique output from all cmds. (as decorator)
    # TODO(sjbush): Maintain a symbolic link to current (i.e., self) script
    #   and log files (named with timestamp).

    def __init__(
            self,
            cmds=None, fall_through=False,
            job_name=None, odir=None,
            template_path=None, timestamp=None,
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

        # set the output dir (if given)
        # -- children must worry about input_dir
        self.input_dir = None
        self.output_dir = odir

        # set an rx dict for general use
        # -- see '_write_dir_shortcuts' for example
        self.rx = {}

        # set the job name
        if job_name:
            self.job_name = job_name
        else:
            self.job_name = self.__class__.__name__.lower()

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

    def _setup(self, *args, **kwargs):
        '''Function used to create preset pipes--Implement in children

        This function is used to add commands to preset pipes.
        To avoid potentially masking true exceptions, the base function
        accepts any parameters, but will raise a NotImplementedError.
        '''
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

    def write(self, _file=None, fmt=None):
        '''Write the pipe commands to a file

        Attributes:
            _file: The file to write to
            fmt: The script format. Use if non-standard file extension
        '''

        # this method serves as a wrapper around the method(s) that
        # actually write the file, so that we can accept
        # -- a file name str
        # -- a file handle
        # -- nothing at all! (we create our own name)
        if not _file:
            try:
                _file = os.path.join(
                    self.output_dir,
                    '{}__{}.pbs'.format(self.job_name, self.timestamp)
                )
            except TypeError:
                # no output_dir given
                msg = '[{}.write] No file given and output_dir not set'.format(
                    self.__class__.__name__)
                raise ValueError(msg)

        self._write_script(_file, fmt=fmt)

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

    def _write_cmds(self, fh):
        cmd = self.cmd(cmd_sep="\n\n")

        # parse in directory shortcuts
        # -- sub output_dir first in case it's a sub dir to input_dir
        # -- CLEANUP: del the rx attributes set by _write_dir_shortcuts
        rx_keys = ['output_dir', 'input_dir']
        try:
            for k in rx_keys:
                try:
                    rx, rep = self.rx[k]
                    cmd = rx.sub(rep, cmd)
                    del self.rx[k]
                except KeyError:
                    pass  # specific rx not set
        except AttributeError:
            pass  # rx set in '_write_dir_shortcuts'

        fh.write(cmd)
        fh.write('\n')  # add trailing new line
        return

    def _write_calling_command(self, fh):
        '''Add the calling command to the script

        Ever wonder what command was used to generate data? (or in this
        case, to generate a script to generate data)
        Never fear! Now you, too, can easily recall what you did when
        you were too tired to think. The write-o-matic pastes the
        information directly into the output. Wow!
        '''

        # split the command by kwarg if longer than 80 char
        cmd = ' '.join(sys.argv)
        if len(cmd) > 80:
            cmd = cmd.replace(' -', ' \\\n#   -')

        # reformat the timestamp to something more readable
        time_fmt = datetime.datetime.strptime(
            self.timestamp, "%y%m%d-%H%M%S").strftime('%Y-%m-%d @ %H:%M:%S')

        # parse together the statement and write it
        blurb = '\n# '.join([
            'Created from dir [{cwd}]',
            'on [{timestamp}] with the command',
            cmd
        ]).format(
            cwd=os.getcwd(),
            timestamp=time_fmt,
        )
        fh.write('# {}\n\n'.format(blurb))
        return

    def _write_pbs_template(self, fh):
        self._load_pbs_template()
        fh.write(self.pbs_template + "\n\n")
        return

    def _write_dir_shortcuts(self, fh):
        '''Write repeated dir paths as BASH variables

        Write the input and output directories to the script as
        BASH variables. Also, compile regex for easier replacement later.
        (These regex should be deleted ASAP!)

        Arguments:
            fh: The file handle.
        Return or Raise:
            None.
        '''

        if self.input_dir:
            _input = 'DATA_DIR={}\n'.format(self.input_dir)
            self.rx['input_dir'] = (
                re.compile(r'{}'.format(self.input_dir)),
                '${DATA_DIR}'
            )
        else:
            _input = ''

        # output_dir should at least be None
        if self.output_dir:
            output = 'OUTPUT_DIR={}\n'.format(self.output_dir)
            self.rx['output_dir'] = (
                re.compile(r'{}'.format(self.output_dir)),
                '${OUTPUT_DIR}'
            )
        else:
            output = ''

        if _input or output:
            fh.write('{}{}\n'.format(_input, output))

    @file_or_handle(mode='w')
    def _write_script(self, fh, fmt=None):
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
            self._write_calling_command(fh)
            self._write_dir_shortcuts(fh)
            self._write_cmds(fh)

        self._update_file_permissions(filename)  # nothing if filename = None
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
