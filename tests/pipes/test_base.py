
import os.path
import builtins
import sys
sys.path.append("remsci/")

import unittest
from unittest.mock import patch, Mock, mock_open

from libpipe.pipes.base import BasePipe


import logging
from remsci.lib.utility import customLogging
customLogging.config()
log = logging.getLogger(__name__)


class TestBasePipe(unittest.TestCase):

    def setUp(self):
        patcher = patch('libpipe.cmds.base.BaseCmd')
        self.mock_cmd = patcher.start()
        self.addCleanup(patcher.stop)

    def setup_mock_write(self):
        patcher = patch.object(builtins, 'open',
                               mock_open(read_data='Hello world'),
                               create=True)
        m = patcher.start()
        self.addCleanup(patcher.stop)
        m.return_value.__iter__ = lambda s: s
        m.return_value.__next__ = lambda s: s.readline()
        return m

    def test_do_run_returns_false_if_cmd_output_exists(self):
        self.mock_cmd._has_output = Mock(return_value=True)
        bp = BasePipe()
        self.assertFalse(bp._do_run(self.mock_cmd))

    def test_do_run_returns_true_if_cmd_output_does_not_exist(self):
        self.mock_cmd._has_output = Mock(return_value=False)
        bp = BasePipe()
        self.assertTrue(bp._do_run(self.mock_cmd))

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
        self.setup_mock_write()

        bp = BasePipe()
        bp._write_pbs('pbs_file')
        if not bp.pbs_template:
            self.fail('PBS template not set during write_pbs')

if __name__ == '__main__':
    unittest.main()
