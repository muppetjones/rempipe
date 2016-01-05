
import os.path
import builtins
import sys
sys.path.append("remsci/")

import unittest
from unittest.mock import patch, Mock, MagicMock, mock_open

from libpipe.cmds.base import BaseCmd
from libpipe.pipes.base import BasePipe

import logging
from remsci.lib.utility import customLogging
customLogging.config()
log = logging.getLogger(__name__)


# class TestCountReadsPipe(unittest.TestCase):
#
#     def setUp(self):
#
#         # prevent creation and deletion of symbolic links during testing
#         patcher = patch.object(BasePipe, '_renew_link')
#         self.mock_renew_link = patcher.start()
#         self.addCleanup(patcher.stop)
#
#         # prevent actually calling any commands
#         patcher = patch('subprocess.check_call')
#         self.mock_call = patcher.start()
#         self.addCleanup(patcher.stop)
#
#         class Cmd(BaseCmd):
#
#             INVOKE_STR = 'test cmd'
#             REQ_TYPE = [
#                 [('-foo', ), ('.txt', )]
#             ]
#
#             def __init__(self, *args, foo='bar', **kwargs):
#                 kwargs['-foo'] = foo
#                 super().__init__(*args, **kwargs)
#
#             def output(self):
#                 return [self.kwargs['-foo']]
#
#         self.CMD = Cmd
#
#     def test_add_links_output_of_first_command_to_input_of_next(self):
#         bp = BasePipe()
#
#         a = self.CMD()
#         b = self.CMD()
#         bp.add(a, b)
#
#         self.assertEqual(a.output, b.input)
#
#     def test_add_double_check_link_integrity(self):
#         bp = BasePipe()
#
#         a = self.CMD(foo='meh')
#         b = self.CMD()
#         bp.add(a, b)
#
#         self.assertEqual(a.output(), ['meh'])
#         self.assertEqual(b.input(), ['meh'])
#
#         a.kwargs['-foo'] = 'hehe'
#         self.assertEqual(b.input(), ['hehe'])
#
#     def test_link_integrity_maintained_through_cmd_output(self):
#         bp = BasePipe()
#
#         a = self.CMD(foo='meh')
#         b = self.CMD()
#         bp.add(a, b)
#
#         a.kwargs['-foo'] = 'hehe.txt'
#         self.assertEqual(b.cmd(readable=False), '{} {} {}'.format(
#             self.CMD.INVOKE_STR, '-foo', 'hehe.txt'))
