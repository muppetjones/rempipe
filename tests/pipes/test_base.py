
import os.path
import builtins
import sys
sys.path.append("remsci/")

import unittest
from unittest.mock import patch, Mock, MagicMock, mock_open, call

from remsci.lib.utility import path

from libpipe.cmds.base import BaseCmd
from libpipe.pipes.base import BasePipe

import logging
from remsci.lib.utility import customLogging
customLogging.config()
log = logging.getLogger(__name__)


class TestBasePipe_setup(unittest.TestCase):

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

        # prevent removal of any files
        patcher = patch('os.remove')
        self.mock_remove = patcher.start()
        self.addCleanup(patcher.stop)

        # prevent actually calling any commands
        patcher = patch('os.stat')
        self.mock_osstat = patcher.start()
        self.addCleanup(patcher.stop)

        # prevent actually calling any commands
        patcher = patch('os.chmod')
        self.mock_oschmod = patcher.start()
        self.addCleanup(patcher.stop)

    def setup_mock_read(self):
        patcher = patch.object(builtins, 'open',
                               new=mock_open(read_data='Hello world'),
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


class TestBasePipe(TestBasePipe_setup):

    #
    #   Initialization tests
    #

    def test_init_sets_default_template_path(self):
        bp = BasePipe()
        self.assertIsNotNone(bp.pbs_template_path)

    def test_init_sets_given_template_path(self):
        templ_path = '~/template.pbs'
        bp = BasePipe(template_path=templ_path)
        self.assertEqual(bp.pbs_template_path, templ_path)

    def test_init_with_force_true_replaces_do_run_with_lambda_function(self):

        bp = BasePipe(force=True)
        lf = lambda: True
        self.assertIsInstance(bp._do_run, type(lf))

    def test_init_with_force_true_replaces_do_run_in_instance(self):
        self.mock_cmd._has_output = Mock(return_value=True)
        bp_force1 = BasePipe(force=True)
        bp = BasePipe(force=False)
        bp_force2 = BasePipe(force=True)

        lf = lambda: True
        self.assertNotIsInstance(bp._do_run, type(lf))
        self.assertIsInstance(bp_force2._do_run, type(lf))
        self.assertIsInstance(bp_force1._do_run, type(lf))

    def test_init_saves_timestamp_if_given(self):
        kw = {'timestamp': '151012-162900'}
        cmd = BasePipe(**kw)
        self.assertEqual(cmd.timestamp, kw['timestamp'])

    def test_init_saves_a_timestamp_if_not_given(self):
        cmd = BasePipe()
        self.assertRegex(cmd.timestamp, '\d{6}-\d{6}')

    def test_init_defaults_job_name_to_None(self):
        cmd = BasePipe()
        self.assertIsNone(cmd.job_name)

    def test_init_saves_job_name(self):
        job_name = 'foo'
        cmd = BasePipe(job_name=job_name)
        self.assertEqual(cmd.job_name, job_name)

    #
    #   Do / Force run checks
    #

    def test_do_run_calls_has_output_cmd(self):
        self.mock_cmd._has_output = Mock(return_value=True)
        bp = BasePipe()
        bp._do_run(self.mock_cmd)

        self.mock_cmd._has_output.assert_called_once_with()

    def test_do_run_returns_false_if_cmd_output_exists(self):
        self.mock_cmd._has_output = Mock(return_value=True)
        bp = BasePipe()
        self.assertFalse(bp._do_run(self.mock_cmd))

    def test_do_run_returns_true_if_cmd_output_does_not_exist(self):
        self.mock_cmd._has_output = Mock(return_value=False)
        bp = BasePipe()
        self.assertTrue(bp._do_run(self.mock_cmd), 'Will not run cmd')

    def test_do_run_returns_true_if_force_is_true_even_if_no_output(self):
        self.mock_cmd._has_output = Mock(return_value=True)
        bp = BasePipe(force=True)
        self.assertTrue(bp._do_run(self.mock_cmd))

    #
    #   Add commands
    #

    def test_add_extends_command_list(self):
        bp = BasePipe()
        a = self.mock_cmd()
        b = self.mock_cmd()

        bp.add(a)
        bp.add(b)

        self.assertEqual(bp.cmds, [a, b])

    def test_add_supports_chaining(self):
        bp = BasePipe()
        a = self.mock_cmd()
        b = self.mock_cmd()

        bp.add(a).add(b)

        self.assertEqual(bp.cmds, [a, b])

    def test_add_accepts_multiple_commands(self):
        bp = BasePipe()
        a = self.mock_cmd()
        b = self.mock_cmd()

        bp.add(a, b)

        self.assertEqual(bp.cmds, [a, b])

    def test_add_links_commands(self):
        bp = BasePipe()
        a = self.mock_cmd()
        b = self.mock_cmd()

        bp.add(a, b)

        a.link.assert_called_once_with(b)

    def test_add_sets_timestamp_on_each_cmd(self):

        bp = BasePipe()
        a = self.mock_cmd()
        b = self.mock_cmd()

        bp.add(a, b)

        self.assertEqual(bp.timestamp, a.timestamp)
        self.assertEqual(bp.timestamp, b.timestamp)


class TestBasePipe_CmdInterface(TestBasePipe_setup):

    '''Test BasePipe Cmd Interface for usage within another pipe'''

    class CmdPipe(BasePipe):

        NAME = 'cmdpipe'

    def setUp(self):
        super().setUp()

        a = self.mock_cmd()
        b = self.mock_cmd()
        cmd = BasePipe()
        cmd.add(a, b)

        self.cmd = cmd

    def sample_pipe(self):
        a = self.mock_cmd()
        b = self.mock_cmd()
        pipe = BasePipe()
        pipe.add(a, b)
        return pipe

    def test_pipe_may_be_used_as_command(self):
        a = self.mock_cmd()
        b = self.mock_cmd()
        subpipe = BasePipe()
        subpipe.add(a, b)

        pipe = BasePipe()
        pipe.add(a, subpipe)  # should not raise

    def test_subpipe_should_return_expected_output(self):
        a = Mock(output=lambda: list('abc'))
        b = Mock(output=lambda: list('def'))
        c = Mock(output=lambda: list('beg'))

        subpipe = BasePipe()
        subpipe.add(a, b)

        pipe = BasePipe()
        pipe.add(c, subpipe)  # should not raise

        self.assertEqual(pipe.output(), subpipe.output())

    def test_subpipe_input_set_as_expected(self):
        a = Mock(output=lambda: list('abc'))
        b = Mock(output=lambda: list('def'))
        c = Mock(output=lambda: list('beg'))

        def ugh(x):
            x.input = c.output
        c.link = Mock(side_effect=ugh)

        subpipe = BasePipe()
        subpipe.add(a, b)

        pipe = BasePipe()
        pipe.add(c, subpipe)  # should not raise

        c.link.assert_called_once_with(subpipe)
        self.assertEqual(subpipe.input(), c.output())

    #
    #   Output tests
    #

    def test_output_returns_output_list_from_last_command(self):
        cmd = self.mock_cmd()
        cmd.output = lambda: list('abc')

        pipe = self.sample_pipe()
        pipe.add(cmd)

        self.assertEqual(pipe.output(), cmd.output())

    def test_output_returns_unique_combined_list_from_all_if_passed_true(self):

        pipe = BasePipe()

        a = Mock(output=lambda: list('abc'))
        b = Mock(output=lambda: list('def'))
        c = Mock(output=lambda: list('beg'))

        pipe.add(a, b, c)

        self.assertEqual(pipe.output(from_all=True), list('abcdefg'))

    #
    #   Link tests
    #

    def test_link_sets_dest_input_to_src_output(self):
        a = self.sample_pipe()
        b = self.sample_pipe()
        a.link(b)

        self.assertEqual(a.output, b.input)

    def test_link_returns_dest_object(self):
        a = self.sample_pipe()
        b = self.sample_pipe()
        c = a.link(b)

        self.assertNotEqual(a, b)
        self.assertNotEqual(a, c)
        self.assertEqual(b, c)

    def test_link_chaining(self):
        a = self.sample_pipe()
        b = self.sample_pipe()
        c = self.sample_pipe()
        d = a.link(b).link(c)

        self.assertEqual(b.output, c.input)
        self.assertEqual(d, c)

    #
    #   Command tests
    #

    def test_cmd_passes_all_execptions_from_self_cmds(self):
        pipe = self.sample_pipe()
        pipe.cmds[-1].cmd.side_effect = TabError('catch me')

        with self.assertRaises(TabError):
            pipe.cmd()

    def test_cmd_returns_string_of_all_cmds(self):
        self.fail()


class TestBasePipe_Write(TestBasePipe_setup):

    def setUp(self):
        super().setUp()
        self.mock_write = self.setup_mock_write()

    def test_write_script_sets_pbs_file_using_job_name_and_timestamp(self):

        bp = BasePipe(job_name='foo')
        bp.write_script()
        self.assertEqual(bp.pbs_file, 'foo_{}.pbs'.format(bp.timestamp))

    def test_write_script_sets_pbs_file_directory(self):

        bp = BasePipe(job_name='foo')
        bp.write_script(directory='~')
        self.assertEqual(
            bp.pbs_file, os.path.join(
                path.protect('~'), 'foo_{}.pbs'.format(bp.timestamp)
            )
        )

    def test_write_loads_pbs_template_on_initial_call(self):
        self.setup_mock_read()
        bp = BasePipe(job_name='foo')
        bp.write_script('pbs_file')

        self.assertIsNotNone(bp.pbs_template)
        self.assertEqual(bp.pbs_template, 'Hello world')

    def test_write_does_not_load_pbs_template_on_subsequent_calls(self):
        m_open = self.setup_mock_read()
        bp = BasePipe(job_name='foo')
        with patch.object(bp, '_BasePipe__write_pbs'):
            bp.write_script('pbs_file')
            bp.write_script('pbs_file')

        self.assertEqual(m_open.call_count, 1, 'PBS template init > 1x')

    def test_write_attempts_to_update_permissions_on_pbs_script(self):
        bp = BasePipe(job_name='foo')
        bp.write_script('pbs_file')

        self.mock_oschmod.assert_called_once_with(
            bp.pbs_file, self.mock_osstat().st_mode.__or__())


class TestBasePipe_Run(TestBasePipe_setup):

    def setUp(self):
        super().setUp()
        self.mock_write = self.setup_mock_write()

    def create_pipe_to_run(self):
        bp = BasePipe(job_name='job_name')
        bp.pbs_file = 'script.pbs'
        patcher = patch.object(os.path, 'isfile', return_value=True)
        self.mock_isfile = patcher.start()
        self.addCleanup(patcher.stop)
        return bp

    def test_run_raises_AttributeError_if_pbs_file_not_set(self):

        bp = BasePipe()
        with self.assertRaisesRegex(AttributeError, r'file.*not set'):
            bp.run()

    def test_run_raises_OSError_if_pbs_file_not_found(self):
        bp = BasePipe()
        bp.pbs_file = 'script.pbs'
        with patch.object(os.path, 'isfile', return_value=False):
            with self.assertRaises(OSError):
                bp.run()

    def test_run_pbs_calls_qsub(self):
        bp = self.create_pipe_to_run()
        bp.run()

        # use the context manager to get the intended mock object
        exp_arg_list = [
            'qsub', '-N', 'job_name', '-o', 'script.log', 'script.pbs']
        with open('file', 'w') as fh:
            self.mock_call.assert_called_once_with(
                exp_arg_list, stdout=fh, stderr=fh)

    def test_run_calls_pbs_script_directly_if_resource_manager_not_found(self):
        self.mock_call.side_effect = [FileNotFoundError, 0]

        bp = self.create_pipe_to_run()
        bp.run()

        with open('file', 'w') as fh:
            self.mock_call.assert_called_with(
                ['script.pbs'], stdout=fh, stderr=fh)

    def test_run_renews_symbolic_links_based_on_job_name(self):
        bp = self.create_pipe_to_run()
        bp.run()

        expected_calls = [
            call(bp.job_name + os.path.splitext(o)[1], o)
            for o in bp.output
        ] + [call(bp.job_name + '.pbs', bp.pbs_file)]

        self.mock_renew_link.assert_has_calls(expected_calls, any_order=True)

if __name__ == '__main__':
    unittest.main()
