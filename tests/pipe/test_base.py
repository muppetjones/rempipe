'''Test the PipeBase class'''

import io
import random
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


class TestPipeBase(PipeBaseTestCase):

    def test_CmdBase_is_mocked(self):
        pass
        # pipe.CmdBase()  # should not raise

    def test_init_saves_timestamp_if_given(self):
        kwargs = {'timestamp': '151012-162900'}
        pipe = PipeBase(**kwargs)
        self.assertEqual(pipe.timestamp, kwargs['timestamp'])

    def test_init_saves_a_timestamp_if_not_given(self):
        pipe = PipeBase()
        self.assertRegex(pipe.timestamp, '\d{6}-\d{6}')

    def test_init_creates_a_dummy_cmd_for_linking(self):
        '''Test that init creates a cmd obj to use for linking'''

        with mock.patch('libpipe.pipe.base.CmdDummy') as mock_dummy:
            PipeBase()
        mock_dummy.assert_called_once_with()

    def test_init_passes_input_to_dummy(self):
        '''Test that dummy output matches the pipe input list'''

        pipe = PipeBase(input=list('abc'))
        self.assertEqual(pipe.dummy.output(), list('abc'))

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


class TestPipeBase_write(PipeBaseTestCase):

    '''Test ability to write a script'''

    def setUp(self):
        super().setUp()
        self.mock_write = self.setup_mock_write()

    def test_write_opens_given_file(self):
        cmds = self.get_n_cmds(3)
        pipe = PipeBase(cmds)
        ofile = 'fake_file.txt'

        pipe.write(ofile)
        self.mock_write.assert_called_once_with(ofile, 'w')

    def test_write_accepts_file_handle(self):
        cmds = self.get_n_cmds(3)
        pipe = PipeBase(cmds)
        pipe.add(cmds)

        with io.StringIO() as sh:
            pipe.write(sh)  # should not raise!
            sh.seek(0)
            result = sh.read()

        self.assertEqual(result, pipe.cmd())


# ENDFILE
