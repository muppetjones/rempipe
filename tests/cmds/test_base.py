import unittest

from textwrap import dedent
from libpipe.cmds.base import BaseCmd


import logging
from remsci.lib.utility import customLogging
customLogging.config()
log = logging.getLogger(__name__)


class TmpCmd(BaseCmd):

    NAME = 'tmp'
    INVOKE_STR = 'tmp'

    attributes = {
        '-f': 'FILE\tInput file',
        '-o': 'FILE\tOutput file',
        '-n': 'INT\tA number',
    }
    DEFAULTS = {
        '-n': 5,
    }

    REQ_KWARGS = ['-f']
    REQ_ARGS = 0


class TestBaseCmds(unittest.TestCase):

    def setUp(self):
        pass

    def test_raises_ValueError_if_required_kwarg_not_given(self):
        with self.assertRaises(ValueError):
            self.cmd = TmpCmd()

    def test_raises_ValueError_if_missing_needed_number_args(self):
        TmpCmd.REQ_ARGS = 2
        kwargs = {'-f': 'req_kwarg'}
        with self.assertRaises(ValueError):
            self.cmd = TmpCmd('only_one_cmd', **kwargs)

    def test_no_exceptions_with_expected_arguments(self):
        TmpCmd.REQ_ARGS = 2
        args = ['first_arg', 'second_arg']
        kwargs = {'-f': 'req_kwarg'}

        self.cmd = TmpCmd(*args, **kwargs)  # no error raised

    def test_defaults_set_on_init(self):
        kwargs = {'-f': 'req_kwarg'}
        cmd = TmpCmd(**kwargs)

        self.assertEqual(cmd.kwargs['-n'], TmpCmd.DEFAULTS['-n'])

    def test_defaults_overridden_if_args_given(self):
        kwargs = {'-f': 'req_kwarg', '-n': 8}
        cmd = TmpCmd(**kwargs)

        self.assertNotEqual(cmd.kwargs['-n'], TmpCmd.DEFAULTS['-n'])
        self.assertEqual(cmd.kwargs['-n'], 8)

    def test_hyphens_added_to_kwargs_if_omitted_during_init(self):

        kwargs = {'f': 'req_kwarg', 'n': 8, 'a': 'foo'}
        cmd = TmpCmd(**kwargs)

        kwargs = {'-' + k: v for k, v in kwargs.items()}
        self.assertDictEqual(kwargs, cmd.kwargs)

    def test_double_hyphens_added_to_str_kwargs_if_omitted_during_init(self):

        kwargs = {'f': 'req_kwarg', 'foo': 'bar'}
        cmd = TmpCmd(**kwargs)

        kwargs = {k: v for k, v in TmpCmd.DEFAULTS.items()}
        kwargs.update({'-f': 'req_kwarg', '--foo': 'bar'})
        self.assertEqual(kwargs, cmd.kwargs)

    def test_trubase_returns_path_wo_dir_or_extn(self):

        basename = BaseCmd._trubase('~/test/path/with/file.name.txt')
        self.assertEqual(basename, 'file.name')

    def test_defaults_unchanged_after_init(self):

        defaults = {}
        defaults.update(TmpCmd.DEFAULTS)

        kwargs = {'f': 'req_kwarg', 'n': 'a'}
        cmd = TmpCmd(**kwargs)

        # Ensure we're deep copying defaults when we set kwargs
        # 1) Check against expected (set above)
        # 2) Check our local copy worked (very basic control)
        # 3) Check object defaults not changed
        # 4) Check that the kwargs are not equal to the defaults
        self.assertEqual(TmpCmd.DEFAULTS['-n'], 5)
        self.assertEqual(TmpCmd.DEFAULTS['-n'], defaults['-n'])
        self.assertEqual(TmpCmd.DEFAULTS['-n'], cmd.DEFAULTS['-n'])
        self.assertNotEqual(TmpCmd.DEFAULTS['-n'], cmd.kwargs['-n'])

    def test_cmd_returns_expected_cmd_string(self):

        kwargs = {'-f': 'req_kwarg', '-n': 'a'}
        cmd = TmpCmd(**kwargs)

        expected_kwargs = ' '.join([
            '{} {}'.format(k, v)
            for k, v in sorted(kwargs.items())
        ])
        expected_args = ' '.join(TmpCmd.ARGS)
        expected_flags = ' '.join(TmpCmd.FLAGS)
        expected_cmd = ' '.join(filter(None, [
            TmpCmd.INVOKE_STR,
            expected_flags,
            expected_kwargs,
            expected_args,
        ]))

        self.assertEqual(
            cmd.cmd(readable=False).rstrip(), expected_cmd.rstrip())

    def test_help_returns_expected_name_synopsis_description(self):

        kwargs = {'-f': 'req_kwarg', '-n': 'a'}
        cmd = TmpCmd(**kwargs)

        expected_help = dedent('''
            NAME
            \t{}

            SYNOPSIS
            \t{}

            DESCRIPTION
            \t{}
        ''').format(
            TmpCmd.NAME,
            TmpCmd.HELP_DICT['synopsis'],
            TmpCmd.HELP_DICT['description'],
        ).lstrip()

        self.assertTrue(
            cmd.help().startswith(expected_help),
            'Help does not start with expected text'
        )

    def test_help_contains_expected_args_text(self):
        kwargs = {'-f': 'req_kwarg', '-n': 'a'}
        cmd = TmpCmd(**kwargs)

        expected_help = '\tFILE\tExample input one'

        self.assertTrue(
            expected_help in cmd.help(),
            'Help does not contain expected args text'
        )

    def test_help_contains_expected_kwargs_text(self):
        kwargs = {'-f': 'req_kwarg', '-n': 'a'}
        cmd = TmpCmd(**kwargs)

        expected_help = '-f|--file\tFILE\tExample keyword argument'

        self.assertTrue(
            expected_help in cmd.help(),
            'Help does not contain expected args text'
        )

    def test_help_contains_expected_flags_text(self):
        kwargs = {'-f': 'req_kwarg', '-n': 'a'}
        cmd = TmpCmd(**kwargs)

        expected_help = '-v\tExample flag argument'

        self.assertTrue(
            expected_help in cmd.help(),
            'Help does not contain expected args text'
        )
