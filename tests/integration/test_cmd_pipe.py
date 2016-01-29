
import sys
sys.path.append("remsci/")

import unittest
from unittest.mock import patch  # , Mock, MagicMock, mock_open

from libpipe.cmds.base import BaseCmd, CmdAttributes
from libpipe.pipes.base import BasePipe

import logging
from remsci.lib.utility import customLogging
customLogging.config()
log = logging.getLogger(__name__)


sample_attributes = {
    'name': 'test command',
    'synopsis': 'Use to test the CmdAttributes class',
    'description': 'More info about the CmdAttributes class',
    'invoke_str': 'test_command',

    'arguments': [
        (None, 'INPUT', 'Something we need'),
        ('-f', 'FILE', 'Better to be explicit'),
        ('-o', 'FILE', 'Output file'),
        ('-n', int, 'A number'),
        ('--foo', 'FILE', 'verbose arg'),
        ('-v', None, 'A random flag'),
        ('-x', None, 'A random flag'),
        ('-foo', 'FILE', 'fill'),
    ],
    'defaults': {
        '-n': 5,
    },

    'req_args': 1,
    'req_kwargs': ['-f', ],
    'req_type': [
        [(0, ), ('.txt', '.csv', ), ],
    ],
}


class TestBasePipe_IntegrationWith_BaseCmd(unittest.TestCase):

    def setUp(self):

        # prevent creation and deletion of symbolic links during testing
        patcher = patch.object(BasePipe, '_renew_link')
        self.mock_renew_link = patcher.start()
        self.addCleanup(patcher.stop)

        # prevent actually calling any commands
        patcher = patch('subprocess.check_call')
        self.mock_call = patcher.start()
        self.addCleanup(patcher.stop)

        ca = CmdAttributes(**sample_attributes)
        ca.req_args = 0
        ca.req_kwargs = []
        ca.req_type = []
        ca.defaults = {}

        class ModSample(BaseCmd):
            attr = ca

            def output(self):
                return ['file.txt']

        self.CMD = ModSample

    def test_add_links_output_of_first_command_to_input_of_next(self):
        bp = BasePipe()

        a = self.CMD()
        b = self.CMD()
        bp.add(a, b)

        self.assertEqual(a.output, b.input)

    def test_changes_in_link_maintined_in_cmd_output(self):
        self.CMD.attr.req_type = sample_attributes['req_type']
        bp = BasePipe()

        a = self.CMD()
        b = self.CMD()
        a.output = lambda: ['hehe.txt']
        bp.add(a, b)

        # NOTE: The type requirement is a first positional *.txt
        a.kwargs['--foo'] = 'hehe.txt'
        cmd_str = b.cmd(verbose=False)
        self.assertEqual(cmd_str, '{} {}'.format(
            self.CMD.attr.invoke_str, 'hehe.txt'))
