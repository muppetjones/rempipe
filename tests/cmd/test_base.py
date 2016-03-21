
import unittest

from libpipe.cmd.attr import CmdAttributes
from libpipe.cmd.base import CmdInterface, CmdBase


sample_attributes = {
    'name': 'test command',
    'synopsis': 'Use to test the CmdAttributes class',
    'description': 'More info about the CmdAttributes class',
    'invoke': 'test_command',

    'args': [
        (0, 'INPUT', 'Something we need'),
        ('-f', 'FILE', 'Better to be explicit'),
        ('-o', 'FILE', 'Output file'),
        ('-n', int, 'A number'),
        ('--foo', 'FILE', 'verbose arg'),
        ('-v', None, 'A random flag'),
        ('-x', None, 'A random flag'),
    ],
    'defaults': {
        '-n': 5,
    },

    'req_args': 1,
    'req_kwargs': ['-f', ],
    'req_types': [
        [(0, ), ('.txt', '.csv', ), ],
    ],
}

#-----------------------------------------------------------------------------
#   Base test case


class BaseTestCase(unittest.TestCase):

    '''Setup a dummy base command for testing'''

    def setUp(self):
        ca = CmdAttributes(**sample_attributes)
        ca.req_args = 0
        ca.req_kwargs = []
        ca.req_type = []
        ca.defaults = {}

        class IndirectBase(CmdBase):
            attr = ca

            def output(self):
                return ['file.txt']
        self.CMD = IndirectBase
        self.ATTR = ca

        # prevent error logs from occuring during testing
        # -- any log.error calls **should** be expected!
        patcher = patch.object(libpipe.cmds.base.log, 'error')
        patcher.start()
        self.addCleanup(patcher.stop)


#-----------------------------------------------------------------------------
#   Direct Tests


class TestCmdInterface(unittest.TestCase):

    def test_cannot_be_instatiated_directly(self):
        with self.assertRaises(TypeError):
            CmdInterface()


class TestCmdBase(unittest.TestCase):

    def setUp(self):
        self.cmd_attr = CmdAttributes(**sample_attributes)

    def test_cannot_be_instatiated_directly(self):
        with self.assertRaises(TypeError):
            CmdBase()

#-----------------------------------------------------------------------------
#   Indirect Tests


class TestCmdBase_init(BaseTestCase):

    '''All tests related to command initialization'''

    def test_init_sets_defaults(self):
        self.CMD.attr.defaults = {'-n': 4800, }
        cmd = self.CMD()

        self.assertEqual(cmd.kwargs['-n'], self.CMD.attr.defaults['-n'])
