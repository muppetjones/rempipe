'''Test the PipeBase class'''

import unittest
from unittest import mock

from libpipe.cmd.base import CmdBase
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
            wrap=None,
        )

    def get_n_cmds(self, n):
        return [self.create_mock_cmd() for i in range(n)]

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


class TestPipeBase_CmdInterface(unittest.TestCase):

    '''Test that PipeBase implements the CmdInterface

    NOTE: Perhaps move to implementation tests?
    '''
