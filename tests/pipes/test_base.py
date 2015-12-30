
import os.path
import builtins
import sys
sys.path.append("remsci/")

import unittest
from unittest.mock import patch, Mock, MagicMock, mock_open

from libpipe.pipes.base import BasePipe


import logging
from remsci.lib.utility import customLogging
customLogging.config()
log = logging.getLogger(__name__)


class TestBasePipe(unittest.TestCase):

    def setUp(self):

        # avoid testing command objects
        patcher = patch('libpipe.cmds.base.BaseCmd')
        self.mock_cmd = patcher.start()
        self.addCleanup(patcher.stop)

        # prevent creation and deletion of symbolic links during testing
        patcher = patch.object(BasePipe, '_renew_link')
        self.mock_renew_link = patcher.start()
        self.addCleanup(patcher.stop)

        # prevent actually calling any commands
        patcher = patch('subprocess.check_call')
        self.mock_call = patcher.start()
        self.addCleanup(patcher.stop)

    def setup_mock_read(self):
        patcher = patch.object(builtins, 'open',
                               mock_open(read_data='Hello world'),
                               create=True)
        m = patcher.start()
        self.addCleanup(patcher.stop)
        m.return_value.__iter__ = lambda s: s
        m.return_value.__next__ = lambda s: s.readline()
        return m

    def setup_mock_write(self):
        mm = MagicMock()
        mo = mock_open(mm)
        patcher = patch.object(builtins, 'open',
                               mo,
                               create=True)
        m = patcher.start()
        self.addCleanup(patcher.stop)
        return m

    def test_do_run_returns_false_if_cmd_output_exists(self):
        self.mock_cmd._has_output = Mock(return_value=True)
        bp = BasePipe()
        self.assertFalse(bp._do_run(self.mock_cmd))

    def test_do_run_returns_true_if_cmd_output_does_not_exist(self):
        self.mock_cmd = Mock()
        with patch.object(BasePipe, '_has_output', return_value=False):
            bp = BasePipe()
            self.assertTrue(bp._do_run(self.mock_cmd), 'Will not run cmd')

    def test_do_run_returns_true_if_force_is_true_even_if_no_output(self):
        self.mock_cmd._has_output = Mock(return_value=True)
        bp = BasePipe(force=True)
        self.assertTrue(bp._do_run(self.mock_cmd))

    def test_force_true_replaces_do_run_with_lambda_function(self):

        bp = BasePipe(force=True)
        lf = lambda: True
        self.assertIsInstance(bp._do_run, type(lf))

    def test_force_true_only_replaces_local_do_run(self):
        self.mock_cmd._has_output = Mock(return_value=True)
        bp_force1 = BasePipe(force=True)
        bp = BasePipe(force=False)
        bp_force2 = BasePipe(force=True)

        lf = lambda: True
        self.assertNotIsInstance(bp._do_run, type(lf))
        self.assertIsInstance(bp_force2._do_run, type(lf))
        self.assertIsInstance(bp_force1._do_run, type(lf))

    def test_pbs_template_loaded_on_write_pbs_call(self):
        self.setup_mock_read()

        bp = BasePipe()
        bp._write_pbs('pbs_file')
        if not bp.pbs_template:
            self.fail('PBS template not set during write_pbs')

    def test_run_pbs_calls_qsub(self):
        self.setup_mock_write()

        bp = BasePipe()
        bp._run_pbs('job_name', 'pbs_file')

        # use the context manager to get the intended mock object
        exp_arg_list = ['qsub', '-N', 'job_name', '-o', 'log_file', 'pbs_file']
        with open('file', 'w') as fh:
            self.mock_call.assert_called_once_with(
                exp_arg_list, stdout=fh, stderr=fh)

    def test_run_calls_pbs_script_directly_if_resource_manager_not_found(self):
        self.setup_mock_write()
        self.mock_call.side_effect = [FileNotFoundError, 0]

        bp = BasePipe()
        bp.run(job_name='job_name', pbs_file='pbs_file')

        with open('file', 'w') as fh:
            self.mock_call.assert_called_with(
                ['pbs_file'], stdout=fh, stderr=fh)

if __name__ == '__main__':
    unittest.main()
