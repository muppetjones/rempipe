import unittest

from libpipe.cmds.base import BaseCmd


import logging
from remsci.lib.utility import customLogging
customLogging.config()
log = logging.getLogger(__name__)


class TmpCmd(BaseCmd):

    bin_name = 'tmp'

    attributes = {
        '-f': 'FILE\tInput file',
        '-o': 'FILE\tOutput file',
        '-n': 'INT\tA number',
    }
    defaults = {
        '-n': 5,
    }

    required_kwargs = ['-f']
    required_args = 0


class TestBaseCmds(unittest.TestCase):

    def setUp(self):
        pass

    def test_raises_ValueError_if_required_kwarg_not_given(self):
        with self.assertRaises(ValueError):
            self.cmd = TmpCmd()

    def test_raises_ValueError_if_missing_needed_number_args(self):
        TmpCmd.required_args = 2
        kwargs = {'-f': 'req_kwarg'}
        with self.assertRaises(ValueError):
            self.cmd = TmpCmd('only_one_cmd', **kwargs)

    def test_no_exceptions_with_expected_arguments(self):
        TmpCmd.required_args = 2
        args = ['first_arg', 'second_arg']
        kwargs = {'-f': 'req_kwarg'}

        self.cmd = TmpCmd(*args, **kwargs)  # no error raised

    def test_defaults_set_on_init(self):
        kwargs = {'-f': 'req_kwarg'}
        cmd = TmpCmd(**kwargs)

        self.assertEqual(cmd.kwargs['-n'], TmpCmd.defaults['-n'])

    def test_defaults_overridden_if_args_given(self):
        kwargs = {'-f': 'req_kwarg', '-n': 8}
        cmd = TmpCmd(**kwargs)

        self.assertNotEqual(cmd.kwargs['-n'], TmpCmd.defaults['-n'])
        self.assertEqual(cmd.kwargs['-n'], 8)

    def test_hyphens_added_to_kwargs_if_omitted_during_init(self):

        kwargs = {'f': 'req_kwarg', 'n': 8, 'a': 'foo'}
        cmd = TmpCmd(**kwargs)

        kwargs = {'-' + k: v for k, v in kwargs.items()}
        self.assertDictEqual(kwargs, cmd.kwargs)

    def test_double_hyphens_added_to_str_kwargs_if_omitted_during_init(self):

        kwargs = {'f': 'req_kwarg', 'foo': 'bar'}
        cmd = TmpCmd(**kwargs)

        kwargs = {k: v for k, v in TmpCmd.defaults.items()}
        kwargs.update({'-f': 'req_kwarg', '--foo': 'bar'})
        self.assertEqual(kwargs, cmd.kwargs)
