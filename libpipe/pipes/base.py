import abc
import os
import os.path
import stat
import subprocess
import time
from itertools import chain

from remsci.lib.utility import path
from remsci.lib.decorators import file_or_handle


import libpipe.templates
from libpipe.utility.exceptions import RempipeExceptionMeta

import logging
log = logging.getLogger(__name__)


PIPE_PBS_TEMPLATE = ''
DO_RUN = False


class BasePipe(object):

    '''Define a group of commands.

    Should interface JUST like a command.

    '''

    _PBS_TEMPLATE = None

    #
    #   Custom Exceptions
    #

    class QueueSubmissionError(FileNotFoundError):
        pass

    class GenericPipeError(ValueError):
        ERRMSG = {
            'job_name': 'Pipe requires a job_name',
        }

        def __init__(self, key, *args, **kwargs):
            try:
                super().__init__(self.ERRMSG[key], *args, **kwargs)
            except KeyError:
                super().__init__(key, *args, **kwargs)

    #
    #   Magic methods
    #

    def __init__(
        self, *args,
        force=False,
        job_name=None,
        template_path=None,
        timestamp=None,
        **kwargs
    ):

        # replace do_run with dummy function
        if force:
            self._do_run = lambda x: True

        # initialize template information
        if template_path:
            self.pbs_template_path = template_path
        else:
            self.pbs_template_path = os.path.join(
                os.path.dirname(libpipe.templates.__file__), 'template.pbs')
        self.pbs_template = None

        # save a timestamp, if passed
        if timestamp:
            self.timestamp = timestamp
        else:
            self.timestamp = time.strftime("%y%m%d-%H%M%S")

        if job_name:
            self.job_name = job_name
        else:
            # raise self.GenericPipeError('job_name')
            self.job_name = self.__class__.__name__

        # initialize commands
        self.cmds = []

    #
    #   Properties
    #
    @property
    def qid(self):
        try:
            with open(self.output['qid'], 'r') as qid_fh:
                qid = qid_fh.read().lstrip().rstrip()
            return qid
        except OSError:
            return None

    @property
    def name(self):
        try:
            return self.NAME
        except AttributeError:
            return self.__class__.__name__

    def _create_pbs_file_name(self, dir_str):

        try:
            self.pbs_file = os.path.join(
                path.protect(dir_str),
                '_'.join([self.job_name, self.timestamp]) + '.pbs'
            )
        except (TypeError, AssertionError):
            raise AssertionError('Pipe attribute "job_name" is undefined')
        return

    #
    #   Command Handling
    #

    def add(self, *cmds):
        self.cmds.extend(cmds)

        for cmd, next_cmd in zip(self.cmds[:-1], self.cmds[1:]):
            cmd.link(next_cmd)
            cmd.timestamp = self.timestamp
        self.cmds[-1].timestamp = self.timestamp

        return self

    #
    #   Command Interface
    #

    def output(self, from_all=False):
        if not from_all:
            return self.cmds[-1].output()
        else:
            olist_nonuniq = chain.from_iterable(
                cmd.output() for cmd in self.cmds)
            olist_uniq = []
            for output in olist_nonuniq:
                if output not in olist_uniq:
                    olist_uniq.append(output)
            return olist_uniq

    def link(self, cmd):
        cmd.input = self.output
        return cmd

    def cmd(self, verbose=True):
        try:
            return self._get_cmd_str(verbose=verbose)
        except:
            raise

    def _has_output(self):
        '''Match BaseCmd Interface -- always return False'''
        return False

    #
    #   Script Handling
    #

    def _do_run(self, cmd):
        return False if cmd._has_output() else True

    def write_script(self, directory=''):

        # ensure we have the job_name
        try:
            self._create_pbs_file_name(directory)
        except AttributeError:
            raise

        # for now, only support pbs
        self.__load_template()
        try:
            self.__write_pbs()
            self.__update_pbs_permissions()
        except:
            raise

    def __write_pbs(self):

        # write pbs file
        with open(self.pbs_file, 'w') as oh:
            oh.write(self.pbs_template + "\n")
            self._write_commands(oh)
        return

    @file_or_handle(mode='w')
    def _write_commands(self, fh):
        fh.write(self._get_cmd_str())
        return

    def _get_cmd_str(self, verbose=True):
        # static text
        if verbose:
            omit_msg = '# Output found. Skip.\n# '
            comment_str = 'echo "Running {}..."\n'
        else:
            omit_msg = '# '
            comment_str = ''

        cmd_list = []
        for i, cmd in enumerate(self.cmds):
            try:
                cmd_str = cmd.cmd(verbose=verbose)  # must call directly!!
            except:
                raise
            if self._do_run(cmd):
                prefix = comment_str.format(cmd.name)
            else:
                prefix = omit_msg
                cmd_str = cmd_str.replace('\n', '\n#')
            cmd_list.append(prefix + cmd_str)
        return '\n\n'.join(cmd_list)

    def __load_template(self):
        if not self.pbs_template:
            with open(self.pbs_template_path, 'r') as ih:
                self.pbs_template = ih.read()

    def __update_pbs_permissions(self):
        '''Make file executable by all parties
        http://stackoverflow.com/a/12792002/977342
        '''
        file_stat = os.stat(self.pbs_file)
        st_all_exec = stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH
        os.chmod(self.pbs_file, file_stat.st_mode | st_all_exec)

    def run(self, *args, **kwargs):

        # Ensure the PBS script has been created
        # -- the attribute should have been set and the file created in
        #    'write_script'.
        try:
            if not os.path.isfile(self.pbs_file):
                raise OSError(
                    'PBS script {} does not exist'.format(self.pbs_file))
        except AttributeError:
            raise AttributeError('PBS file not set. Write script first.')

        try:
            self.__run_queue(*args, **kwargs)
        except self.QueueSubmissionError:
            # Most likely due to qsub not being installed
            self.__run_local(*args, **kwargs)
        finally:
            self._renew_links()
        return

    def __run_queue(self, log_file=None):
        '''Run given PBS file

        Submits the given pbs_file to queue (TODO: security check).
        Also, creates symbolic links to the log, pbs, and qid files
        based on the job name.

        Arguments:
            job_name    The name to run the job as.
            pbs_file    The pbs file path
            log_file    A log file path. Defaults to pbs file w/ log extension.
        Returns:
            The id of the job.
        Raises:
            FileNotFoundError: Most likely due to missing qsub.
        '''

        qid_file = self.pbs_file.replace('pbs', 'qid')
        if not log_file:
            log_file = self.pbs_file.replace('pbs', 'log')

        # run pbs file
        try:
            qsub_command = [
                'qsub', '-N', self.job_name, '-o', log_file, self.pbs_file, ]
            with open(qid_file, 'w') as qid_fh:
                retcode = subprocess.check_call(
                    qsub_command, stdout=qid_fh, stderr=qid_fh)
        except FileNotFoundError:
            os.remove(qid_file)  # don't keep qsub id file if no job
            raise self.QueueSubmissionError('Unable to submit to qsub')

        # update the link to reflect latest run
        self.output = [qid_file, log_file]

        return

    def __run_local(self, log_file=None):
        '''Run local PBS file.

        Run this only in the case that _run_pbs raises FileNotFoundError.
        Executes the given pbs_file directly from the shell.
        Also, creates symbolic links to the log and pbs based on the job name.

        Arguments:
            job_name    The name to run the job as.
            pbs_file    The pbs file path
            log_file    A log file path. Defaults to pbs file w/ log extension.
        Returns:
            None
        Raises:
            None
        '''

        if not log_file:
            log_file = self.pbs_file.replace('pbs', 'log')

        # run pbs file locally without a resource manager
        run_command = [self.pbs_file, ]  # should be able to run file alone
        with open(log_file, 'w') as lh:
            try:
                retcode = subprocess.check_call(
                    run_command, stdout=lh, stderr=lh)
            except subprocess.CalledProcessError as err:
                lh.write(str(err))

        # update the link to reflect latest run
        self.output = [log_file]

        return

    def _renew_links(self):
        link_base = os.path.join(os.path.dirname(self.pbs_file), self.job_name)

        try:
            to_link = [self.pbs_file] + self.output
        except AttributeError:
            return  # no output; do nothing
        else:
            for f in to_link:
                extn = os.path.splitext(f)[1]
                self._renew_link(link_base + extn, f)
        return

    # def _run_pbs(self, job_name, log_file=''):
    # def _run_local(self, job_name, log_file=''):

    def _renew_link(self, link_name, src_file):
        if os.path.isfile(link_name):
            os.unlink(link_name)
        os.symlink(src_file, link_name)


class PresetPipe(BasePipe):

    '''A pipe object that defines the interface for creating a preset pipe

    Define a group of commands that should be run together.
    Intended for repetative processes, such as align, sort, count.
    The interface should allow it to act like a command and a pipe.
    '''

    __metaclass__ = abc.ABCMeta

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        try:
            self._setup(*args, **kwargs)
        except:
            raise

    @abc.abstractmethod
    def _setup(self, *args, **kwargs):
        pass
