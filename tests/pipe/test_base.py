'''Test the PipeBase class'''

import unittest
from unittest import mock

from libpipe.cmd.interface import CmdInterface
from libpipe.pipe.base import PipeBase

import logging
log = logging.getLogger(__name__)


class PipeBaseTestCase(unittest.TestCase):

    def setUp(self):

        # avoid testing the CmdBase class
        patcher = mock.patch(
            'libpipe.cmd.base.CmdBase',
            new=self.create_mock_cmd()
        )
        self.mock_cmd = patcher.start()
        self.addCleanup(patcher.stop)

    def create_mock_cmd(self):
        '''Use to create individual mock instances (vs. self.mock_cmd)'''
        return mock.MagicMock(
            link=mock.Mock(),
            cmd=mock.Mock(return_value='cmd --foo'),
            output=lambda: [],
        )

    def get_n_cmds(self, n):
        return [self.create_mock_cmd() for i in range(n)]


class TestPipeBase(PipeBaseTestCase):

    def test_CmdBase_is_mocked(self):
        pass
        # pipe.CmdBase()  # should not raise

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

        cmds = self.get_n_cmds(3)

        pipe = PipeBase()
        pipe.add(cmds[0], cmds[1]).add(cmds[-1])

        self.assertEqual(pipe.cmds, cmds)

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
