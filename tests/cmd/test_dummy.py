'''Test the CmdDummy class'''

import unittest
from unittest import mock

from libpipe.cmd.base import CmdBase
from libpipe.cmd.dummy import CmdDummy


import logging
log = logging.getLogger(__name__)


class TestCmdDummy(unittest.TestCase):

    def test_init_does_not_need_parameters(self):
        CmdDummy()  # should not raise

    def test_init_does_NOT_call_super_init(self):
        self.assertNotEqual(CmdDummy.__init__, CmdBase.__init__)
        with mock.patch.object(CmdBase, '__init__') as mock_init:
            CmdDummy()

        self.assertEqual(mock_init.call_count, 0)

    def test_output_contains_all_pargs(self):
        '''Test we pass pos args through (implicit 'input' check)'''
        args = list('abc')
        cmd = CmdDummy(*args)
        output = cmd.output()
        self.assertEqual(output, args)

    def test_output_contains_all_kwargs_sorted_values(self):
        '''Test we pass kw args through (implicit 'input' check)'''
        kwargs = {'a': 'z', 'foo': 'bar', 'abc': '123', 'hello': 'world'}
        expected = sorted(list(kwargs.values()))

        # Dict order is arbitrary, but this should be consistent--test a ton!
        # -- we don't care which test fails, so no subTest
        for i in range(100):
            cmd = CmdDummy(**kwargs)
            output = cmd.output()
            self.assertEqual(output, expected)

    def test_output_lists_args_before_kwargs(self):
        args = list('abc')
        kwargs = {'a': 'z', 'foo': 'bar', 'abc': '123', 'hello': 'world'}
        cmd = CmdDummy(*args, **kwargs)
        output = cmd.output()
        expected = args + sorted(list(kwargs.values()))
        self.assertEqual(output, expected)

    def test_cmd_returns_empty_string(self):
        cmd = CmdDummy()
        self.assertEqual(cmd.cmd(), '')

    def test_link_previous_overrides_dummy_input(self):
        '''Ensure we can still use dummy in a pipeline'''
        args = list('abc')

        mock_cmd = mock.MagicMock(
            CmdBase, link=CmdBase.link, output=lambda: list('xyz'))
        cmd = CmdDummy(*args)
        mock_cmd.link(mock_cmd, cmd)  # link expects 'self'

        self.assertNotEqual(cmd.output(), args)
        self.assertEqual(cmd.output(), mock_cmd.output())

    def test_link_next_correctly_sets_next_input(self):
        '''Ensure dummy will link to subsequent command as expected'''
        args = list('abc')
        mock_cmd = mock.MagicMock(CmdBase)()  # create instance
        cmd = CmdDummy(*args)
        cmd.link(mock_cmd)

        self.assertEqual(mock_cmd.input(), args)
