
import os
import os.path
import subprocess

import logging
log = logging.getLogger(__name__)

PIPE_PBS_TEMPLATE = ''
DO_RUN = False

import libpipe.templates


class BasePipe(object):

    _PBS_TEMPLATE = None

    def __init__(self, force=False):

        # replace do_run with dummy function
        if force:
            self.do_run = lambda x: True

        # initialize template information
        self.pbs_template_path = os.path.join(
            os.path.dirname(libpipe.templates.__file__), 'template.pbs')
        self.pbs_template = None

        # initialize commands
        self.cmds = []

    def do_run(self, cmd):
        log.debug('in do_run')
        return False if cmd._has_output() else True

    def write_pbs(self, pbs_file):
        if not self.pbs_template:
            with open(self.pbs_template_path, 'r') as ih:
                self.pbs_template = ih.read()

        # write pbs file
        with open(pbs_file, 'w') as oh:
            oh.write(self.pbs_template + "\n")
            oh.write("which samtools\n")

            omit_msg = '# Output found. Skip.\n# '
            comment_str = 'echo "Running {}..."\n'

            for cmd in self.cmds:
                cmd_str = str(cmd) + "\n\n"
                if self.do_run(cmd):
                    oh.write(comment_str.format(cmd.name))
                    oh.write(cmd_str)
                else:
                    oh.write(omit_msg + cmd_str)

    def _run_pbs(self, job_name, pbs_file, log_file=''):

        qid_file = pbs_file.replace('pbs', 'qid')
        if not log_file:
            log_file = pbs_file.replace('pbs', 'log')

        # run pbs file
        qsub_command = ['qsub', '-N', job_name, '-o', log_file, pbs_file, ]
        with open(qid_file, 'w') as qid_fh:
            log.debug(' '.join(qsub_command))
            retcode = subprocess.check_call(qsub_command, stdout=qid_fh)
            log.debug('retcode: {}'.format(retcode))

        # get the qsub job id
        with open(qid_file, 'r') as qid_fh:
            qid = qid_fh.read()

        # link the log file (last run)
        link_base = os.path.join(os.path.dirname(pbs_file), job_name)
        self._renew_link(link_base + '.log', log_file)
        self._renew_link(link_base + '.pbs', pbs_file)
        self._renew_link(link_base + '.qid', qid_file)

        return qid

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
