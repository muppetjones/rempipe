import unittest
import libpipe
from unittest.mock import patch


from libpipe.cmds.base import BaseCmd


import logging
from remsci.lib.utility import customLogging
customLogging.config()
log = logging.getLogger(__name__)


class CmdSample(BaseCmd):

    NAME = 'tmp'
    INVOKE_STR = 'tmp'
    ARGUMENTS = [
        (None, 'INPUT', 'Something we need'),
        ('-f', 'FILE', 'Better to be explicit'),
        ('-o', 'FILE', 'Output file'),
        ('-n', 'INT', 'A number'),
        ('v', None, 'A random flag'),
    ]

    DEFAULTS = {
        '-n': 5,
    }

    REQ_KWARGS = []  # ['-f']
    REQ_ARGS = 0


class TestBaseCmds(unittest.TestCase):

    def setUp(self):
        class ModSample(CmdSample):
            pass
        self.CMD = ModSample

        # prevent error logs from occuring during testing
        patcher = patch.object(libpipe.cmds.base.log, 'error')
        patcher.start()
        self.addCleanup(patcher.stop)

    def sample(self):
        kwargs = {'-f': 'req_kwarg', '-n': 'a'}
        return self.CMD(**kwargs)

    def test_init_raises_AttributeError_if_keyword_args_not_expanded(self):
        '''Test for common mistake of passing keyword dict directly'''
        kw = {'-f': 0, '-n': 2}

        with self.assertRaises(AttributeError):
            self.CMD(kw)

    def test_cmd_raises_ValueError_if_required_kwarg_not_given(self):
        self.CMD.REQ_KWARGS = ['-f']
        cmd = self.CMD()  # should not raise

        with self.assertRaises(ValueError):
            cmd.cmd()

    def test_cmd_raises_ValueError_if_missing_req_number_args(self):
        self.CMD.REQ_ARGS = 2
        kwargs = {'-f': 'req_kwarg'}
        cmd = self.CMD('only_one_cmd', **kwargs)  # shold not raise

        with self.assertRaises(ValueError):
            cmd.cmd()

    def test_cmd_raises_ValueError_if_missing_req_AND_args(self):
        self.CMD.REQ_KWARGS = [('-f', '-o', '-n')]
        self.CMD.DEFAULTS = {}
        kw = {'-f': 0, '-n': 2}

        cmd = self.CMD(**kw)  # should not raise
        with self.assertRaises(ValueError):
            cmd.cmd()

    def test_cmd_raises_ValueError_if_missing_all_req_XOR_args(self):
        self.CMD.REQ_KWARGS = [['-f', '-o', '-n']]
        self.CMD.DEFAULTS = {}

        cmd = self.CMD(**{'-h': 1})
        with self.assertRaises(ValueError):
            cmd.cmd()

    def test_cmd_raises_ValueError_if_more_than_one_req_XOR_arg_given(self):
        self.CMD.REQ_KWARGS = [['-f', '-o', '-n']]
        self.CMD.DEFAULTS = {}
        kw = {'-f': 0, '-n': 2}

        cmd = self.CMD(**kw)
        with self.assertRaises(ValueError):
            cmd.cmd()

    def test_no_exceptions_with_expected_arguments(self):
        self.CMD.REQ_ARGS = 2
        args = ['first_arg', 'second_arg']
        kwargs = {'-f': 'req_kwarg'}

        cmd = self.CMD(*args, **kwargs)  # no error raised
        cmd.cmd()  # no error raised

    def test_defaults_set_on_init(self):
        cmd = self.CMD()

        self.assertEqual(cmd.kwargs['-n'], self.CMD.DEFAULTS['-n'])

    def test_defaults_overridden_if_args_given(self):
        kwargs = {'-f': 'req_kwarg', '-n': 8}
        cmd = self.CMD(**kwargs)

        self.assertNotEqual(cmd.kwargs['-n'], self.CMD.DEFAULTS['-n'])
        self.assertEqual(cmd.kwargs['-n'], 8)

    def test_hyphens_added_to_kwargs_if_omitted_during_init(self):

        kwargs = {'f': 'req_kwarg', 'n': 8, 'a': 'foo'}
        cmd = self.CMD(**kwargs)

        kwargs = {'-' + k: v for k, v in kwargs.items()}
        self.assertDictEqual(kwargs, cmd.kwargs)

    def test_double_hyphens_added_to_str_kwargs_if_omitted_during_init(self):

        kwargs = {'f': 'req_kwarg', 'foo': 'bar'}
        cmd = self.CMD(**kwargs)

        kwargs = {k: v for k, v in self.CMD.DEFAULTS.items()}
        kwargs.update({'-f': 'req_kwarg', '--foo': 'bar'})
        self.assertEqual(kwargs, cmd.kwargs)

    def test_trubase_returns_path_wo_dir_or_extn(self):

        basename = BaseCmd._trubase('~/test/path/with/file.name.txt')
        self.assertEqual(basename, 'file.name')

    def test_defaults_unchanged_after_init(self):

        defaults = {}
        defaults.update(self.CMD.DEFAULTS)

        kwargs = {'f': 'req_kwarg', 'n': 'a'}
        cmd = self.CMD(**kwargs)

        # Ensure we're deep copying defaults when we set kwargs
        # 1) Check against expected (set above)
        # 2) Check our local copy worked (very basic control)
        # 3) Check object defaults not changed
        # 4) Check that the kwargs are not equal to the defaults
        self.assertEqual(self.CMD.DEFAULTS['-n'], 5)
        self.assertEqual(self.CMD.DEFAULTS['-n'], defaults['-n'])
        self.assertEqual(self.CMD.DEFAULTS['-n'], cmd.DEFAULTS['-n'])
        self.assertNotEqual(self.CMD.DEFAULTS['-n'], cmd.kwargs['-n'])

    def test_cmd_returns_expected_cmd_string(self):

        kwargs = {'-f': 'req_kwarg', '-n': 'a'}
        cmd = self.sample()

        expected_kwargs = ' '.join([
            '{} {}'.format(k, v)
            for k, v in sorted(kwargs.items())
        ])
        expected_args = None
        expected_flags = None
        expected_cmd = ' '.join(filter(None, [
            self.CMD.INVOKE_STR,
            expected_flags,
            expected_kwargs,
            expected_args,
        ]))

        self.assertEqual(
            cmd.cmd(readable=False).rstrip(), expected_cmd.rstrip())

    def test_help_text_contains_name(self):
        cmd = self.sample()
        self.assertIn(self.CMD.NAME, cmd.help())

    def test_help_contains_expected_args_text(self):
        cmd = self.sample()
        expected_help = 'Something we need'

        self.assertIn(
            expected_help, cmd.help(),
            'Help does not contain expected args text'
        )

    def test_help_contains_expected_kwargs_text(self):
        cmd = self.sample()
        expected_help = 'Better to be explicit'

        self.assertIn(
            expected_help, cmd.help(),
            'Help does not contain expected kwargs text'
        )

    def test_help_contains_expected_flags_text(self):
        cmd = self.sample()
        expected_help = 'A random flag'

        self.assertIn(
            expected_help, cmd.help(),
            'Help does not contain expected flags text'
        )

    def test_link_sets_dest_input_to_src_output(self):
        a = self.sample()
        b = self.sample()
        a.link(b)

        self.assertEqual(a.output, b.input)

    def test_link_returns_dest_object(self):
        a = self.sample()
        b = self.sample()
        c = a.link(b)

        self.assertNotEqual(a, b)
        self.assertNotEqual(a, c)
        self.assertEqual(b, c)

    def test_link_chaining(self):
        a = self.sample()
        b = self.sample()
        c = self.sample()
        d = a.link(b).link(c)

        self.assertEqual(b.output, c.input)
        self.assertEqual(d, c)
