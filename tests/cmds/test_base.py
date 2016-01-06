import unittest
import libpipe
from unittest.mock import patch, Mock


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

    #
    #   Initialization tests
    #

    def test_init_sets_defaults(self):
        cmd = self.CMD()

        self.assertEqual(cmd.kwargs['-n'], self.CMD.DEFAULTS['-n'])

        def test_init_defaults_overridden_if_args_given(self):
            kwargs = {'-f': 'req_kwarg', '-n': 8}
            cmd = self.CMD(**kwargs)

            self.assertNotEqual(cmd.kwargs['-n'], self.CMD.DEFAULTS['-n'])
            self.assertEqual(cmd.kwargs['-n'], 8)

    def test_init_saves_timestamp_if_given(self):
        kw = {'timestamp': '151012-162900'}
        cmd = self.CMD(**kw)
        self.assertEqual(cmd.timestamp, kw['timestamp'])

    def test_init_saves_a_timestamp_if_not_given(self):
        cmd = self.CMD()
        self.assertRegexpMatches(cmd.timestamp, '\d{6}-\d{6}')

    def test_init_raises_AttributeError_if_keyword_args_not_expanded(self):
        '''Test for common mistake of passing keyword dict directly'''
        kw = {'-f': 0, '-n': 2}

        with self.assertRaises(AttributeError):
            self.CMD(kw)

    def test_init_raises_AttributeError_if_unrecognized_kwargs_given(self):
        self.fail()

    def test_init_raises_IndexError_if_too_many_args_given(self):
        self.fail()

    def test_init_raises_ArgumentError_if_unrecognized_flag_given(self):
        self.fail()

    def test_init_sets_unrecognized_arguments_if_open_eq_True(self):
        self.fail()

    def test_init_adds_hyphens_to_kwargs_if_omitted_during(self):

        kwargs = {'f': 'req_kwarg', 'n': 8, 'a': 'foo'}
        cmd = self.CMD(**kwargs)

        kwargs = {'-' + k: v for k, v in kwargs.items()}
        self.assertDictEqual(kwargs, cmd.kwargs)

    def test_init_adds_double_hyphens_to_str_kwargs_if_omitted(self):

        kwargs = {'f': 'req_kwarg', 'foo': 'bar'}
        cmd = self.CMD(**kwargs)

        kwargs = {k: v for k, v in self.CMD.DEFAULTS.items()}
        kwargs.update({'-f': 'req_kwarg', '--foo': 'bar'})
        self.assertEqual(kwargs, cmd.kwargs)

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

    #
    #   Misc
    #

    def test_trubase_returns_path_wo_dir_or_extn(self):

        basename = BaseCmd._trubase('~/test/path/with/file.name.txt')
        self.assertEqual(basename, 'file.name')

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

    #
    #   Help tests
    #

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

    #
    #   Link tests
    #

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

    def test_cmd_matches_linked_input_without_error(self):
        self.CMD.REQ_KWARGS = ['-f']
        self.CMD.REQ_TYPE = [
            [('-f', ), ('.txt', )],
        ]

        a = self.CMD()
        a.output = lambda: ['file.txt']
        b = self.CMD()
        a.link(b)

        b.cmd()  # should not raise
        self.assertEqual(b.kwargs['-f'], a.output()[0])

    def test_cmd_raises_TypeError_if_req_not_fulfiled_by_link(self):
        self.CMD.REQ_KWARGS = ['-f']
        self.CMD.REQ_TYPE = [
            [('-f', ), ('.txt', )],
        ]

        a = self.CMD()
        b = self.CMD()
        a.link(b)

        with self.assertRaises(TypeError):
            b.cmd()

    def test_cmd_raises_TypeError_if_unable_to_use_linked_input(self):
        self.CMD.REQ_KWARGS = ['-f']
        self.CMD.REQ_TYPE = [
            [('-f', ), ('.txt', )],
        ]

        a = self.CMD()
        a.output = lambda: ['file.csv']
        b = self.CMD()
        a.link(b)

        with self.assertRaises(TypeError):
            b.cmd()

    def test_cmd_raises_TypeError_if_linked_input_unused_w_exact_match(self):
        self.CMD.REQ_KWARGS = ['-1']
        self.CMD.REQ_TYPE = [
            [('-1', '-2'), ('.txt', ), False],
        ]

        a = self.CMD()
        a.output = lambda: ['file.txt']
        b = self.CMD()
        a.link(b)

        with self.assertRaises(TypeError):
            b.cmd()

    def test_cmd_sets_linked_input_to_correct_flag_if_not_match_any(self):
        self.CMD.REQ_TYPE = [
            [('-1', '-2'), ('.txt', ), False],
            [('-U', ), ('.txt', ), False],
        ]

        a = self.CMD()
        a.output = lambda: ['file.txt']
        b = self.CMD()
        a.link(b)

        b.cmd()  # should not raise!
        self.assertEqual(b.kwargs['-U'], a.output()[0])

    def test_cmd_raises_TypeError_if_input_not_callable(self):
        self.CMD.REQ_KWARGS = ['-1']
        self.CMD.REQ_TYPE = [
            [('-1', '-2'), ('.txt', ), False],
        ]

        a = self.CMD()
        a.output = ['file.txt']
        b = self.CMD()
        a.link(b)

        with self.assertRaises(TypeError):
            b.cmd()

    def test_cmd_sets_linked_input_to_pos_arg_if_req(self):
        self.CMD.REQ_ARGS = 1
        self.CMD.DEFAULTS = []
        self.CMD.REQ_TYPE = [
            [('-1', '-2'), ('.txt', ), False],
        ]

        a = self.CMD()
        b = self.CMD()
        a.output = lambda: ['file.txt', 'file.1', 'file2.txt']
        a.link(b)

        b.cmd()  # should not raise
        self.assertEqual(b.kwargs, {'-1': 'file.txt', '-2': 'file2.txt'})
        self.assertIn('file.1', b.args)

    def test_cmd_sets_linked_input_to_pos_arg_if_given(self):
        self.CMD.REQ_ARGS = 0
        self.CMD.DEFAULTS = []
        self.CMD.REQ_TYPE = [
            [(0, ), ('.txt', )],
        ]

        a = self.CMD()
        b = self.CMD()
        a.output = lambda: ['file.txt', 'file.1', 'file2.txt']
        a.link(b)

        b.cmd()  # should not raise
        self.assertEqual(b.args, ['file.txt'])

    #
    #   Required argument tests
    #

    def test_cmd_raises_IndexError_if_missing_req_number_args(self):
        self.CMD.REQ_ARGS = 2
        args = ['only_one_arg']
        cmd = self.CMD(*args)  # shold not raise

        with self.assertRaises(IndexError):
            cmd.cmd()

    def test_no_exceptions_with_expected_arguments(self):
        self.CMD.REQ_ARGS = 2
        args = ['first_arg', 'second_arg']

        cmd = self.CMD(*args)  # no error raised
        cmd.cmd()  # no error raised

    def test_cmd_raises_AttributeError_if_required_kwarg_not_given(self):
        self.CMD.REQ_KWARGS = ['-f']
        cmd = self.CMD()  # should not raise

        with self.assertRaises(AttributeError):
            cmd.cmd()

    def test_cmd_reqAND_does_not_raise_if_none_of_group_are_given(self):
        self.CMD.REQ_KWARGS = [('-f', '-o', '-n')]  # req all three
        self.CMD.DEFAULTS = {}
        kw = {'-x': 'something different'}

        cmd = self.CMD(**kw)  # should not raise
        cmd.cmd()  # should NOT raise!!

    def test_cmd_raises_AttributeError_if_missing_req_AND_args(self):
        self.CMD.REQ_KWARGS = [('-f', '-o', '-n')]  # req all three
        self.CMD.DEFAULTS = {}
        kw = {'-f': 0, '-n': 2}

        cmd = self.CMD(**kw)  # should not raise
        with self.assertRaises(AttributeError):
            cmd.cmd()

    def test_cmd_raises_AttributeError_if_missing_all_req_XOR_args(self):
        self.CMD.REQ_KWARGS = [['-f', '-o', '-n']]
        self.CMD.DEFAULTS = {}

        cmd = self.CMD(**{'-h': 1})
        with self.assertRaises(AttributeError):
            cmd.cmd()

    def test_cmd_raises_AttributeError_if_gt_1_req_XOR_arg_given(self):
        self.CMD.REQ_KWARGS = [['-f', '-o', '-n']]
        self.CMD.DEFAULTS = {}
        kw = {'-f': 0, '-n': 2}

        cmd = self.CMD(**kw)
        with self.assertRaises(AttributeError):
            cmd.cmd()

    def test_cmd_checks_for_additional_requirements_defined_by_child(self):
        m = Mock()
        self.CMD._additional_requirements = m
        cmd = self.CMD()
        cmd.cmd()

        m.assert_called_once_with()

    #
    #   File type tests
    #

    def test_cmd_type_checks_positional_arguments(self):
        self.CMD.REQ_TYPE = [
            [(0, 1), ('.txt', )],
        ]
        a = self.sample()
        a.args = ['req_args.txt']
        # should not raise; implicit check against missing [1]
        a.cmd()

    def test_cmd_does_not_raise_TypeError_if_checked_kwarg_not_given(self):
        self.CMD.REQ_TYPE = [
            [('-x', ), ('.txt', )],
        ]
        a = self.sample()
        a.cmd()  # should not raise

    def test_cmd_raises_TypeError_if_wrong_filetype_given(self):
        self.CMD.REQ_TYPE = [
            [('-f', ), ('.txt', )],
        ]
        a = self.sample()
        with self.assertRaises(TypeError):
            a.cmd()

    def test_cmd_raises_TypeError_for_bad_pos_arg_type(self):
        self.CMD.REQ_TYPE = [
            [(0, 1), ('.txt', )],
        ]
        a = self.sample()
        a.args = ['req_args.txt', 'req_args.csv']
        with self.assertRaises(TypeError):
            a.cmd()

    def test_cmd_does_not_raise_if_expected_filetype_given(self):
        self.CMD.REQ_TYPE = [
            [('-f', ), ('.txt', )],
        ]
        kwargs = {'-f': 'req_kwarg.txt'}
        a = self.sample()
        a.kwargs = kwargs
        a.cmd()  # should not raise

    #
    #   Execution
    #

    def test_filter_type_returns_files_with_expected_type(self):
        args = ['seq.1.fq', 'seq.2.fq', 'seq.txt']
        extn = ['.fq']

        self.assertEqual(self.CMD._filter_by_type(args, extn), args[:-1])

    def test_filter_type_returns_empty_if_not_found(self):
        args = ['seq.1.fq', 'seq.2.fq', 'seq.txt']
        extn = ['.csv']

        self.assertFalse(
            self.CMD._filter_by_type(args, extn),
            'Filter returned non-empty list'
        )
