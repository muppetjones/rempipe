
import os
import os.path
import stat
import subprocess


import libpipe.templates

import logging
log = logging.getLogger(__name__)


PIPE_PBS_TEMPLATE = ''
DO_RUN = False


class BasePipe(object):

    _PBS_TEMPLATE = None

    def __init__(self, force=False):

        # replace do_run with dummy function
        if force:
            self._do_run = lambda x: True

        # initialize template information
        self.pbs_template_path = os.path.join(
            os.path.dirname(libpipe.templates.__file__), 'template.pbs')
        self.pbs_template = None

        # initialize commands
        self.cmds = []

    def _do_run(self, cmd):
        return False if self._has_output(cmd) else True

    def write_script(self, *args, **kwargs):

        # for now, only support pbs
        self.pbs_file = self._write_pbs(*args, **kwargs)
        self._chmod_plus_x(self.pbs_file)

    def _write_pbs(self, pbs_file):
        if not self.pbs_template:
            with open(self.pbs_template_path, 'r') as ih:
                self.pbs_template = ih.read()

        # write pbs file
        with open(pbs_file, 'w') as oh:
            oh.write(self.pbs_template + "\n")

            omit_msg = '# Output found. Skip.\n# '
            comment_str = 'echo "Running {}..."\n'

            for cmd in self.cmds:
                cmd_str = str(cmd) + "\n\n"
                if self._do_run(cmd):
                    oh.write(comment_str.format(cmd.name))
                    oh.write(cmd_str)
                else:
                    cmd_str = cmd_str.replace('\n', '\n#')
                    oh.write(omit_msg + cmd_str)

        return pbs_file

    def _chmod_plus_x(self, filepath):
        '''Make file executable by all parties
        http://stackoverflow.com/a/12792002/977342
        '''
        st = os.stat(filepath)
        st_all_exec = stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH
        os.chmod(filepath, st.st_mode | st_all_exec)

    def run(self, *args, **kwargs):

        try:
            return self._run_pbs(*args, **kwargs)
        except FileNotFoundError:
            # Most likely due to qsub not being installed
            return self._run_local(*args, **kwargs)

    def _run_pbs(self, job_name, pbs_file, log_file=''):
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

        qid_file = pbs_file.replace('pbs', 'qid')
        if not log_file:
            log_file = pbs_file.replace('pbs', 'log')

        # run pbs file
        try:
            qsub_command = ['qsub', '-N', job_name, '-o', log_file, pbs_file, ]
            with open(qid_file, 'w') as qid_fh:
                retcode = subprocess.check_call(
                    qsub_command, stdout=qid_fh, stderr=qid_fh)
        except FileNotFoundError:
            os.remove(qid_file)  # don't keep qsub id file if no job
            raise  # expected if qsub not found

        # get the qsub job id
        with open(qid_file, 'r') as qid_fh:
            qid = qid_fh.read()

        # update the link to reflect latest run
        link_base = os.path.join(os.path.dirname(pbs_file), job_name)
        self._renew_link(link_base + '.log', log_file)
        self._renew_link(link_base + '.pbs', pbs_file)
        self._renew_link(link_base + '.qid', qid_file)

        return qid

    def _run_local(self, job_name, pbs_file, log_file=''):
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
            log_file = pbs_file.replace('pbs', 'log')

        # run pbs file locally without a resource manager
        run_command = [pbs_file, ]  # should be able to run file alone
        with open(log_file, 'w') as lh:
            try:
                retcode = subprocess.check_call(
                    run_command, stdout=lh, stderr=lh)
            except subprocess.CalledProcessError as err:
                lh.write(str(err))

        # update the link to reflect latest run
        link_base = os.path.join(os.path.dirname(pbs_file), job_name)
        self._renew_link(link_base + '.log', log_file)
        self._renew_link(link_base + '.pbs', pbs_file)

    def _renew_link(self, link_name, src_file):
        if os.path.isfile(link_name):
            os.unlink(link_name)
        os.symlink(src_file, link_name)

    def _has_output(self, cmd):
        '''Checks for expected output from a command

        Arguments:
            cmd     Command object.
        Returns:
            True if ALL output is found.
            False otherwise
        '''

        for f in cmd.output:
            if not os.path.isfile(f):
                return False
        return True
