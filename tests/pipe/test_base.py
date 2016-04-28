'''Test the PipeBase class'''

import io
import os.path
import random
import re
import subprocess
from unittest import mock

from tests.base import LibpipeTestCase  # includes read and write mock

from libpipe.cmd.interface import CmdInterface
from libpipe.pipe.base import PipeBase

import logging
log = logging.getLogger(__name__)


class PipeBaseTestCase(LibpipeTestCase):

    def setUp(self):

        # avoid testing the CmdBase class
        patcher = mock.patch(
            'libpipe.cmd.base.CmdBase',
            new=self.create_mock_cmd()
        )
        self.mock_cmd = patcher.start()
        self.addCleanup(patcher.stop)

    def create_mock_cmd(self, i=0):
        '''Use to create individual mock instances (vs. self.mock_cmd)'''
        return mock.MagicMock(
            link=mock.Mock(),
            cmd=mock.Mock(return_value='cmd --foo={}'.format(i)),
            output=mock.Mock(return_value=[])
        )

    def get_n_cmds(self, n):
        return [self.create_mock_cmd(i) for i in range(n)]

    def mock_protect(self):
        patcher = mock.patch(
            'libpipe.util.path.protect', side_effect=lambda x: x)
        self.addCleanup(patcher.stop)
        return patcher.start()


class TestPipeBase(PipeBaseTestCase):

    def test_CmdBase_is_mocked(self):
        pass
        # pipe.CmdBase()  # should not raise

    def test_init_raises_TypeError_if_unexpected_kwarg_and_no_setup(self):
        '''Test we don't get arguments we don't expect

        NOTE: We may want to allow this in later implementations,
            e.g., to pass arguments through to commands
        '''
        u = list('rcdefblame')
        kwargs = {'shot_through_the_heart': u[1:-5]}
        with self.assertRaises(TypeError):
            PipeBase(**kwargs)

    def test_init_saves_timestamp_if_given(self):
        kwargs = {'timestamp': '151012-162900'}
        pipe = PipeBase(**kwargs)
        self.assertEqual(pipe.timestamp, kwargs['timestamp'])

    def test_init_saves_a_timestamp_if_not_given(self):
        pipe = PipeBase()
        self.assertRegex(pipe.timestamp, '\d{6}-\d{6}')

    def test_init_sets_jobname_to_class_name_if_not_given(self):
        pipe = PipeBase()
        self.assertEqual(pipe.job_name, pipe.__class__.__name__.lower())

    def test_init_sets_jobname_if_given(self):
        pipe = PipeBase(job_name="foo_bar")
        self.assertEqual(pipe.job_name, 'foo_bar')

    def test_init_creates_a_dummy_cmd_for_linking(self):
        '''Test that init creates a cmd obj to use for linking'''

        with mock.patch('libpipe.pipe.base.CmdDummy') as mock_dummy:
            PipeBase()
        mock_dummy.assert_called_once_with()

    def test_init_passes_input_to_dummy(self):
        '''Test that dummy output matches the pipe input list'''

        pipe = PipeBase(input=list('abc'))
        self.assertEqual(pipe.dummy.output(), list('abc'))

    def test_init_adds_cmds_from_list(self):
        '''Test commands added and linked directly through init'''
        cmds = self.get_n_cmds(3)
        pipe = PipeBase(input=list('abc'), cmds=cmds)
        self.assertEqual(pipe.cmds, cmds)
        self.assertEqual(pipe.cmds[0].input(), list('abc'))

    def test_input_returns_list_from_dummy_cmd(self):
        pipe = PipeBase(input=list('abc'))
        _input = pipe.input()
        self.assertEqual(_input, list('abc'))
        self.assertEqual(pipe.dummy.input(), _input)

    def test_add_stores_cmds_in_list_in_order_added(self):
        '''Test basic add cmd (implicit test for multiple adding methods)'''

        cmds = self.get_n_cmds(5)
        pipe = PipeBase()

        # tests chaining
        # tests adding commands individually, as *args, and as a list
        pipe.add(cmds[0]).add(cmds[1], cmds[2]).add(cmds[3:])

        self.assertEqual(pipe.cmds, cmds)

    def test_add_sets_timestamp_on_each_cmd(self):
        cmds = self.get_n_cmds(3)
        pipe = PipeBase()
        pipe.add(cmds)

        for cmd in pipe.cmds:
            self.assertEqual(pipe.timestamp, cmd.timestamp)

    def test_add_links_all_cmds_when_called(self):
        '''Test that commands are re-linked on successive calls'''
        n = 3
        cmds = self.get_n_cmds(n)

        pipe = PipeBase()
        for i in range(n):
            pipe.add(cmds[i])

        max_n_call = n - 1
        for i in range(n):
            self.assertEqual(cmds[i].link.call_count, max_n_call - i)

    def test_add_links_dummy_to_first_command(self):
        '''Test that dummy output matchs first command input

        NOTE: Due to decorating and mocking, can't compare functions directly
        '''

        _input = list('abc')
        cmd = self.mock_cmd()
        pipe = PipeBase(input=_input)

        pipe.add(cmd)

        self.assertEqual(cmd.input(), _input)
        self.assertEqual(pipe.dummy.output(), cmd.input())

    def test_setup_raises_NotImplemented_error(self):
        pipe = PipeBase()
        with self.assertRaises(NotImplementedError):
            pipe._setup()

    def test_setup_called_during_init(self):
        with mock.patch.object(PipeBase, '_setup') as mock_setup:
            PipeBase()
        mock_setup.assert_called_once_with()

    def test_setup_not_called_if_cmds_passed_to_init(self):
        cmds = self.get_n_cmds(3)
        with mock.patch.object(PipeBase, '_setup') as mock_setup:
            PipeBase(cmds=cmds)
        self.assertEqual(mock_setup.call_count, 0, '_setup was called!')

    def test_setup_gets_unused_kwargs(self):
        kwargs = {'input': ['not_passed'], 'give_to_setup': 'passed'}
        with mock.patch.object(PipeBase, '_setup') as mock_setup:
            PipeBase(**kwargs)
        mock_setup.assert_called_once_with(give_to_setup='passed')


class TestPipeBase_CmdInterface(PipeBaseTestCase):

    '''Test that PipeBase implements the CmdInterface

    NOTE: Perhaps move to implementation tests?
    '''

    def test_inherits_from_cmd_interface(self):
        pipe = PipeBase()
        self.assertIsInstance(pipe, CmdInterface)

    #
    #   Cmd
    #

    def test_cmd_raises_ValueError_if_pipe_is_empty(self):
        pipe = PipeBase()
        with self.assertRaises(ValueError):
            pipe.cmd()

    def test_cmd_returns_a_string(self):
        pipe = PipeBase()
        pipe.add(self.get_n_cmds(3))
        cmd_str = pipe.cmd()
        self.assertIsInstance(cmd_str, str)

    def test_cmd_returns_string_with_one_cmd_per_line(self):
        pipe = PipeBase()
        cmds = self.get_n_cmds(3)
        pipe.add(cmds)

        expected = '\n'.join(cmd.cmd() for cmd in cmds)
        cmd_str = pipe.cmd()
        self.assertEqual(cmd_str, expected)

    def test_cmd_uses_cmd_sep_parameter_to_join_cmds(self):
        pipe = PipeBase()
        cmds = self.get_n_cmds(3)
        pipe.add(cmds)

        cmd_str = pipe.cmd(cmd_sep='@')
        self.assertEqual(cmd_str.count('@'), len(cmds) - 1)

    def test_cmd_does_not_block_exceptions_from_cmd(self):

        # this list is not complete, but includes the standard we would
        # expect plus a couple we wouldn't
        exception_list = [
            AttributeError, ValueError, KeyError, OSError, TabError,
        ]

        n = 3
        for err in exception_list:
            which_cmd = random.randrange(0, n)
            pipe = PipeBase()
            cmds = self.get_n_cmds(n)
            pipe.add(cmds)

            cmds[which_cmd].cmd.side_effect = err

            with self.subTest(err=err):
                with self.assertRaises(err):
                    pipe.cmd()

    #
    #   Link
    #

    def test_link_raises_ValueError_if_pipe_is_empty(self):
        cmd = self.mock_cmd()
        pipe = PipeBase()
        with self.assertRaises(ValueError):
            pipe.link(cmd)

    def test_link_returns_given_obj(self):
        '''Test link enables chaining by returning linked object'''
        cmd_add, cmd_link = self.get_n_cmds(2)
        pipe = PipeBase()
        pipe.add(cmd_add)
        link_out = pipe.link(cmd_link)
        self.assertEqual(link_out, cmd_link)

    def test_link_sets_linked_obj_input_to_last_cmd_output(self):
        cmd_add, cmd_link = self.get_n_cmds(2)
        cmd_add.output = lambda: list('abc')
        pipe = PipeBase()
        pipe.add(cmd_add)
        pipe.link(cmd_link)
        self.assertEqual(cmd_link.input(), cmd_add.output())

    def test_pipe_can_be_linked(self):
        pipe1 = PipeBase(input=list('xyz'), fall_through=True)
        pipe2 = PipeBase()
        pipe1.add(*self.get_n_cmds(2))
        pipe2.add(*self.get_n_cmds(2))
        link_ret = pipe1.link(pipe2)

        self.assertEqual(link_ret, pipe2)
        self.assertEqual(pipe2.input(), list('xyz'))

    #
    #   Output
    #

    def test_output_raises_ValueError_if_pipe_is_empty(self):
        pipe = PipeBase()
        with self.assertRaises(ValueError):
            pipe.output()

    def test_output_returns_output_from_last_command(self):
        cmd = self.mock_cmd()
        cmd.output = lambda: list('abc')

        pipe = PipeBase()
        pipe.add(cmd)

        self.assertEqual(pipe.output(), cmd.output())

    def test_fall_through_passes_all_input_to_output(self):
        _input = list('abc')
        pipe = PipeBase(input=_input, fall_through=True)
        pipe.add(*self.get_n_cmds(2))
        self.assertEqual(pipe.output(), _input)


class TestPipeBase_run(PipeBaseTestCase):

    '''Test ability to run commands

    TODO(sjbush): Run commands separately
    TODO(sjbush): Submit job if available
    '''

    def setUp(self):
        super().setUp()

        # prevent actually calling any commands
        patcher = mock.patch('subprocess.check_call')
        self.mock_call = patcher.start()
        self.addCleanup(patcher.stop)

    def test_run_calls_each_cmd_run_if_no_script_file(self):
        '''Test each cmd is run sep. if script_file not set'''
        pipe = PipeBase(cmds=self.get_n_cmds(3))
        pipe.run()

        for cmd in pipe.cmds:
            cmd.run.assert_called_once_with()
        self.assertEqual(self.mock_call.call_count, 0)  # cmd is mocked

    def test_run_calls_subprocess_once_if_script_file_set(self):
        '''Test script_file is executed if set (implicit logfile check)'''

        cmds = self.get_n_cmds(3)
        pipe = PipeBase(cmds=self.get_n_cmds(3))
        pipe.script_file = 'script.pbs'
        pipe.run()

        for cmd in cmds:
            self.assertEqual(cmd.call_count, 0)
        self.mock_call.assert_called_once_with(
            pipe.script_file, stdout=pipe.log_file, stderr=pipe.log_file)

    def test_run_script_raises_CalledProcessError_if_problem(self):
        self.mock_call.side_effect = subprocess.CalledProcessError(1, 'AH!')

        pipe = PipeBase(cmds=self.get_n_cmds(3))
        pipe.script_file = 'script.pbs'
        with self.assertRaises(subprocess.CalledProcessError):
            pipe.run()

    def test_run_script_does_not_call_subprocess_directly(self):
        # assumes default mode of local
        pipe = PipeBase(cmds=self.get_n_cmds(3))
        pipe.script_file = 'script.pbs'
        with mock.patch.object(pipe, '_run_local'):
            pipe.run()
        self.assertEqual(self.mock_call.call_count, 0)

    def test_run_script_raises_ValueError_for_unknown_mode(self):
        pipe = PipeBase(cmds=self.get_n_cmds(3))
        pipe.script_file = 'script.pbs'
        with self.assertRaises(ValueError):
            pipe.run(mode='Dorian')

    def test_run_script_local_called_by_default(self):
        '''Test that run mode is 'local' by default'''
        pipe = PipeBase(cmds=self.get_n_cmds(3))
        pipe.script_file = 'script.pbs'
        with mock.patch.object(pipe, '_run_local') as mock_run_local:
            pipe.run()
        mock_run_local.assert_called_once_with()

    def test_run_script_calls_run_pbs_if_mode_is_pbs(self):
        '''Test that _run_pbs is used when mode=pbs'''

        pipe = PipeBase(cmds=self.get_n_cmds(3))
        pipe.script_file = 'script.pbs'
        with mock.patch.object(pipe, '_run_local') as mock_run_local, \
                mock.patch.object(pipe, '_run_pbs') as mock_run_pbs:
            pipe.run(mode='pbs')
        mock_run_pbs.assert_called_once_with()
        self.assertEqual(mock_run_local.call_count, 0)


class TestPipeBase_write(PipeBaseTestCase):

    '''Test ability to write a script'''

    def setUp(self):
        super().setUp()
        self.mock_write = self.setup_mock_write()

        # from __update_pbs_permissions
        patcher = mock.patch('os.chmod')
        self.mock_oschmod = patcher.start()
        self.addCleanup(patcher.stop)

        patcher = mock.patch('os.stat')
        self.mock_osstat = patcher.start()
        self.addCleanup(patcher.stop)

    def mock_splitext(self, _file):
        '''Mock out splitext--cannot easily get filename from mock handle'''
        split_tuple = os.path.splitext(_file)
        patcher = mock.patch('os.path.splitext')
        mock_obj = patcher.start()
        self.addCleanup(patcher.stop)
        mock_obj.return_value = list(split_tuple)
        return mock_obj

    def test_write_raises_ValueError_if_pipe_is_empty(self):
        pipe = PipeBase()
        with self.assertRaises(ValueError):
            pipe.write('script.pbs')

    #
    #   PBS template setup
    #

    def test_default_pbs_template_actually_exists(self):
        pipe = PipeBase()
        self.assertTrue(
            os.path.exists(pipe.pbs_template_path),
            'Default template ({}) does not exist'.format(
                pipe.pbs_template_path)
        )

    def test_init_sets_default_template_path(self):
        pipe = PipeBase()
        self.assertIsNotNone(pipe.pbs_template_path)

    def test_init_sets_given_template_path(self):
        template_path = '~/template.pbs'
        pipe = PipeBase(template_path=template_path)
        self.assertEqual(pipe.pbs_template_path, template_path)

    def test_pbs_template_not_set_before_call_to_write(self):
        pipe = PipeBase()
        self.assertIsNone(pipe.pbs_template)

    def test_pbs_template_not_loaded_if_non_pbs_extension_found(self):
        '''Test non-pbs file formats do not load the pbs template'''
        _file = 'script.csh'
        self.mock_splitext(_file)
        self.setup_mock_read('Hello world')

        pipe = PipeBase()
        with mock.patch.object(pipe, '_load_pbs_template') as mock_load, \
                mock.patch.object(pipe, '_write_cmds'):
            pipe.write(_file)
        self.assertIsNone(pipe.pbs_template)
        self.assertEqual(mock_load.call_count, 0)

    def test_pbs_template_loaded_once_on_initial_call(self):
        mock_read = self.setup_mock_read('Hello world')
        pipe = PipeBase()

        for i in range(5):
            pipe._load_pbs_template()

        mock_read.assert_called_once_with(pipe.pbs_template_path, 'r')
        self.assertEqual(pipe.pbs_template, 'Hello world')

    #
    #   write
    #

    def test_write_opens_given_file(self):
        '''Assert `write` attempts to open the file'''
        _file = 'script.pbs'
        self.mock_splitext(_file)
        pipe = PipeBase(cmds=self.get_n_cmds(3))

        pipe.write(_file)
        self.mock_write.assert_any_call(_file, 'w')

    def test_write_opens_self_named_file_if_not_given(self):
        self.mock_protect()
        odir = 'foo/bar'
        self.mock_splitext('meh.pbs')
        pipe = PipeBase(cmds=self.get_n_cmds(3), odir=odir)
        _file = os.path.join(odir, '{}__{}.pbs'.format(
            pipe.job_name, pipe.timestamp))

        pipe.write()
        self.mock_write.assert_any_call(_file, 'w')

    def test_write_raises_ValueError_if_no_file_or_output_dir_given(self):
        self.mock_protect()
        pipe = PipeBase(cmds=self.get_n_cmds(3))

        with self.assertRaises(ValueError):
            pipe.write()

    def test_pipe_writes_all_commands_to_handle(self):
        pipe = PipeBase(cmds=self.get_n_cmds(3))

        with io.StringIO() as sh:
            pipe.write(sh)  # should not raise!
            sh.seek(0)
            result = sh.read()

        # The script should/could insert extra new lines or comments
        no_comments = re.sub(r'\#.*\n', '', result)
        clean_result = no_comments.replace('\n\n', '\n').lstrip().rstrip()
        self.assertEqual(clean_result, pipe.cmd().rstrip())

    def test_pbs_template_written_if_pbs_or_shell_file_given(self):
        '''Tests pbs template loaded when needed'''
        pipe = PipeBase(cmds=self.get_n_cmds(3))
        pipe.pbs_template = 'Thou cream faced loon'  # prevent template loading
        mock_splitext = self.mock_splitext('script.pbs')

        for extn in ['.pbs', '.sh']:
            with self.subTest(extn=extn):
                mock_splitext.return_value[1] = extn
                _file = 'script' + extn
                pipe.write(_file)
                self.mock_write().write.assert_any_call(
                    pipe.pbs_template + "\n\n")  # don't forget the new line!

    def test_write_attempts_to_update_permissions_if_filename_given(self):
        _file = 'script.pbs'
        self.mock_splitext(_file)
        pipe = PipeBase(cmds=self.get_n_cmds(3))
        pipe.pbs_template = 'You fustilarian!'  # prevent template loading

        with mock.patch.object(pipe, '_write_cmds'):  # prevent write
            pipe.write(_file)

        self.mock_oschmod.assert_called_once_with(
            self.mock_write().name, self.mock_osstat().st_mode.__or__())

    def test_write_does_not_attempt_permissions_if_non_file_handle(self):
        pipe = PipeBase(cmds=self.get_n_cmds(3))
        pipe.pbs_template = 'You fustilarian!'  # prevent template loading

        with io.StringIO() as sh:
            pipe.write(sh)

        self.assertEqual(self.mock_oschmod.call_count, 0)

    def test_write_stores_filename_if_file(self):
        '''Test that the filename written to is stored by pipe'''
        _file = 'script.pbs'
        self.mock_splitext(_file)
        pipe = PipeBase(cmds=self.get_n_cmds(3))
        pipe.write(_file)
        self.assertEqual(pipe.script_file, self.mock_write().name)

    def test_write_includes_cmd_str_that_generated_the_file(self):
        pipe = PipeBase(cmds=self.get_n_cmds(3))

        fake_argv = ['ABC', '--easy as', '123', '*' * 70]
        with io.StringIO() as sh, mock.patch('sys.argv', fake_argv):
            pipe.write(sh)  # should not raise!
            sh.seek(0)
            result = sh.read()

        # replace consecutive spaces, backslashes, and crunches
        # -- the latter two in case the line is broken apart
        result = re.sub(r'[\\\#\s\n]+', ' ', result)
        self.assertIn(' '.join(fake_argv), result)


# ENDFILE
