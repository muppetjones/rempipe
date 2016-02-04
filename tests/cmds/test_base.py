import unittest
import libpipe
from unittest.mock import patch, Mock


from libpipe.cmds.base import BaseCmd, CmdAttributes


import logging
from logging import RootLogger
from remsci.lib.utility import customLogging
customLogging.config()
log = logging.getLogger(__name__)


# class CmdSample(BaseCmd):
#
#     NAME = 'tmp'
#     INVOKE_STR = 'tmp'
#     ARGUMENTS = [
#         (None, 'INPUT', 'Something we need'),
#         ('-f', 'FILE', 'Better to be explicit'),
#         ('-o', 'FILE', 'Output file'),
#         ('-n', 'INT', 'A number'),
#         ('--foo', 'FILE', 'verbose arg'),
#         ('-v', None, 'A random flag'),
#         ('-x', None, 'A random flag'),
#     ]
#
#     attr.defaults = {
#         '-n': 5,
#     }
#
#     attr.req_kwargs = []  # ['-f']
#     attr.req_args = 0
#
#     def output(self):
#         return ['file.txt']


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


class TestCmdAttributes(unittest.TestCase):

    def setUp(self):
        self.kwargs = dict.copy(sample_attributes)

    def test_init_sets_invoke_str_to_name_if_not_given(self):
        del self.kwargs['invoke_str']
        ca = CmdAttributes(**self.kwargs)
        self.assertEqual(ca.invoke_str, ca.name)

    def test_init_deep_copies_all_attributes(self):
        ca = CmdAttributes(**self.kwargs)
        for k, v in self.kwargs.items():
            with self.subTest(attr=k):
                self.kwargs[k] = 'Something completely different'
                self.assertNotEqual(getattr(ca, k), self.kwargs[k])

    def test_duplicate_creates_new_object(self):
        ca = CmdAttributes(**self.kwargs)
        self.assertNotEqual(ca, ca.duplicate())

    def test_duplicate_deep_copies_object(self):
        ca = CmdAttributes(**self.kwargs)
        dup = ca.duplicate()

        for k, v in ca.__dict__.items():
            with self.subTest(attr=k):
                self.assertEqual(getattr(ca, k), getattr(dup, k))
                setattr(ca, k, 'not equal')
                self.assertNotEqual(getattr(ca, k), getattr(dup, k))


class TestBase(unittest.TestCase):

    def setUp(self):
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
        self.ATTR = ca

        # prevent error logs from occuring during testing
        patcher = patch.object(libpipe.cmds.base.log, 'error')
        patcher.start()
        self.addCleanup(patcher.stop)

    def sample(self, *args, **kwargs):
        dkwargs = {'-f': 'somefile.txt', '-n': 'a'}
        dkwargs.update(kwargs)
        return self.CMD(*args, **dkwargs)

    def linked_samples(self):
        a = self.CMD()
        a.output = lambda: ['file.txt']
        b = self.CMD()
        a.link(b)
        return a, b


class TestBaseCmd_init(TestBase):

    '''BaseCmd initialization tests'''

    def test_BaseCmd_cannot_be_initialized(self):
        with self.assertRaises(TypeError):
            BaseCmd()

    def test_init_raises_NotImplementedError_if_ARGUMENTS_not_set(self):
        '''ID10T error'''

        class Id10tCmd(BaseCmd):

            def output(self):
                pass

        with self.assertRaisesRegex(NotImplementedError, 'attr.arguments'):
            Id10tCmd()

    def test_init_sets_defaults(self):
        self.CMD.attr.defaults = {'-n': 4800, }
        cmd = self.CMD()

        self.assertEqual(cmd.kwargs['-n'], self.CMD.attr.defaults['-n'])

    def test_init_defaults_overridden_if_args_given(self):
        kwargs = {'-f': 'req_kwarg', '-n': 8}
        self.CMD.attr.defaults['-n'] = 25
        cmd = self.CMD(**kwargs)

        self.assertNotEqual(cmd.kwargs['-n'], self.CMD.attr.defaults['-n'])
        self.assertEqual(cmd.kwargs['-n'], 8)

    def test_init_saves_timestamp_if_given(self):
        kw = {'timestamp': '151012-162900'}
        cmd = self.CMD(**kw)
        self.assertEqual(cmd.timestamp, kw['timestamp'])

    def test_init_saves_a_timestamp_if_not_given(self):
        cmd = self.CMD()
        self.assertRegex(cmd.timestamp, '\d{6}-\d{6}')

    def test_init_raises_AttributeError_if_keyword_args_not_expanded(self):
        '''Test for common mistake of passing keyword dict directly'''
        kw = {'-f': 0, '-n': 2}

        with self.assertRaises(IndexError):
            self.CMD(kw)

    def test_init_adds_hyphens_to_kwargs_if_omitted_during(self):

        kwargs = {'f': 'req_kwarg', 'n': 8}
        cmd = self.CMD(**kwargs)

        kwargs = {'-' + k: v for k, v in kwargs.items()}
        self.assertDictEqual(kwargs, cmd.kwargs)

    def test_init_adds_double_hyphens_to_str_kwargs_if_omitted(self):
        kwargs = {'f': 'req_kwarg', 'foo': 'bar'}
        cmd = self.CMD(**kwargs)

        kwargs = {k: v for k, v in self.CMD.attr.defaults.items()}
        kwargs.update({'-f': 'req_kwarg', '--foo': 'bar'})
        self.assertEqual(kwargs, cmd.kwargs)

    def test_defaults_unchanged_after_init(self):
        defaults = {'-n': 805}
        self.CMD.attr.defaults = dict.copy(defaults)

        kwargs = {'f': 'req_kwarg', 'n': 'a'}
        cmd = self.CMD(**kwargs)

        # Ensure we're deep copying defaults when we set kwargs
        # 1) Check against expected (set above)
        # 2) Check object defaults not changed
        # 3) Check that the kwargs are not equal to the defaults
        self.assertEqual(self.CMD.attr.defaults['-n'], defaults['-n'])
        self.assertEqual(self.CMD.attr.defaults['-n'], cmd.attr.defaults['-n'])
        self.assertNotEqual(self.CMD.attr.defaults['-n'], cmd.kwargs['-n'])

    def test_init_adds_args_w_hyphens_as_flags(self):

        args = ['test.txt', '-v', '-x', 'test.xls']
        cmd = self.CMD(*args)

        expected_flags = [arg for arg in args if arg.startswith('-')]
        expected_args = [arg for arg in args if arg not in expected_flags]

        self.assertEqual(cmd.flags, expected_flags)
        self.assertEqual(cmd.args, expected_args)


class TestBaseCmd_requirements(TestBase):

    '''Requirement tests'''

    @unittest.skip('Need to upgrade positional arg handling first')
    def test_init_raises_IndexError_if_too_many_args_given(self):
        with self.assertRaises(IndexError):
            self.CMD('file.txt')

    def test_init_raises_ArgumentError_if_unrecognized_flag_given(self):
        with self.assertRaises(AttributeError):
            self.CMD('--bar')

    def test_init_sets_unrecognized_kwarg_if_strict_eq_False(self):
        kw = {'--unknown_flag': 'haha'}
        cmd = self.CMD(strict=False, **kw)  # should not raise
        self.assertIn('--unknown_flag', cmd.kwargs)

    def test_init_raises_AttributeError_if_unrecognized_kwargs_given(self):
        kw = {'-p': 'unknown kwarg'}
        with self.assertRaises(AttributeError):
            self.CMD(**kw)

    #
    #   Required argument tests
    #

    def test_cmd_raises_IndexError_if_missing_req_number_args(self):
        self.CMD.attr.req_args = 2
        args = ['only_one_arg']
        cmd = self.CMD(*args)  # shold not raise

        with self.assertRaises(IndexError):
            cmd.cmd()

    def test_no_exceptions_with_expected_arguments(self):
        self.CMD.attr.req_args = 2
        self.CMD.attr.req_type = [[(0, ), ('.txt', )]]
        self.CMD.attr.req_kwargs = []
        args = ['first_arg.txt', 'second_arg']

        cmd = self.CMD(*args)  # no error raised
        cmd.cmd()  # no error raised

    def test_cmd_raises_AttributeError_if_required_kwarg_not_given(self):
        self.CMD.attr.req_kwargs = ['-f']
        cmd = self.CMD()  # should not raise

        with self.assertRaises(AttributeError):
            cmd.cmd()

    def test_cmd_reqAND_does_not_raise_if_non_members_are_given(self):
        self.CMD.attr.req_kwargs = [('-f', '-o', '-n')]  # req all three
        self.CMD.attr.defaults = {}
        kw = {'--foo': 'something different'}

        cmd = self.CMD(**kw)  # should not raise
        cmd.cmd()  # should NOT raise!!

    def test_cmd_raises_AttributeError_if_missing_reqAND_args(self):
        self.CMD.attr.req_kwargs = [('-f', '-o', '-n')]  # req all three
        self.CMD.attr.defaults = {}
        kw = {'-f': 0, '-n': 2}

        cmd = self.CMD(**kw)  # should not raise
        with self.assertRaises(AttributeError):
            cmd.cmd()

    def test_cmd_raises_AttributeError_if_missing_all_reqXOR_args(self):
        self.CMD.attr.req_kwargs = [['-f', '-o', '-n']]
        self.CMD.attr.defaults = {}

        cmd = self.CMD(**{'--foo': 1})
        with self.assertRaises(AttributeError):
            cmd.cmd()

    def test_cmd_raises_AttributeError_if_gt_1_req_XOR_arg_given(self):
        self.CMD.attr.req_kwargs = [['-f', '-o', '-n']]
        self.CMD.attr.defaults = {}
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
        self.CMD.attr.req_type = [
            [(0, 1), ('.txt', )],
        ]
        a = self.sample()
        a.args = ['req_args.txt']
        # should not raise; implicit check against missing [1]
        a.cmd()

    def test_cmd_does_not_raise_TypeError_if_checked_kwarg_not_given(self):
        self.CMD.attr.req_type = [
            [('-x', ), ('.txt', )],
        ]
        a = self.sample()
        a.cmd()  # should not raise

    def test_cmd_raises_TypeError_if_wrong_filetype_given(self):
        self.CMD.attr.req_type = [
            [('-f', ), ('.42', )],
        ]
        a = self.sample()
        with self.assertRaises(TypeError):
            a.cmd()

    def test_cmd_raises_TypeError_for_bad_pos_arg_type(self):
        self.CMD.attr.req_type = [
            [(0, 1), ('.txt', )],
        ]
        a = self.sample()
        a.args = ['req_args.txt', 'req_args.csv']
        with self.assertRaises(TypeError):
            a.cmd()

    def test_cmd_does_not_raise_if_expected_filetype_given(self):
        self.CMD.attr.req_type = [
            [('-f', ), ('.txt', )],
        ]
        kwargs = {'-f': 'req_kwarg.txt'}
        a = self.sample()
        a.kwargs = kwargs
        a.cmd()  # should not raise


class TestBaseCmd_misc(TestBase):

    '''Miscellaneous tests; help, etc.'''

    def test_has_output_returns_True_if_output_file_exists(self):
        with patch('os.path.isfile', return_value=True):
            cmd = self.sample()
            self.assertTrue(cmd._has_output())

    def test_has_output_returns_False_if_output_file_does_not_exist(self):
        with patch('os.path.isfile', return_value=False):
            cmd = self.sample()
            self.assertFalse(cmd._has_output())

    def test_trubase_returns_path_wo_dir_or_extn(self):

        basename = BaseCmd._trubase('~/test/path/with/file.name.txt')
        self.assertEqual(basename, 'file.name')

    def test_cmd_returns_expected_cmd_string(self):

        kwargs = {'-f': 'req_kwarg', '-n': 'a'}
        cmd = self.sample(**kwargs)

        expected_kwargs = ' '.join([
            '{} {}'.format(k, v)
            for k, v in sorted(kwargs.items())
        ])
        expected_args = None
        expected_flags = None
        expected_cmd = ' '.join(filter(None, [
            self.CMD.attr.invoke_str,
            expected_flags,
            expected_kwargs,
            expected_args,
        ]))

        self.assertEqual(
            cmd.cmd(verbose=False).rstrip(), expected_cmd.rstrip())

    #
    #   Help tests
    #

    def test_help_text_contains_name(self):
        cmd = self.sample()
        self.assertIn(self.CMD.attr.name, cmd.help())

    def test_help_contains_expected_args_text(self):
        cmd = self.sample()
        help_str = cmd.help()

        for _, _, expected_txt in cmd.attr.arguments:
            with self.subTest(text=expected_txt):
                self.assertIn(
                    expected_txt, help_str,
                )

    def test_fall_through_command_includes_ALL_input_in_output(self):
        class Cmd2(self.CMD):

            def output(self):
                return ['ha!']

        a = self.CMD()
        out = ['something_unexpected', 42]
        a.output = lambda: out
        b = Cmd2(fall_through=True)
        a.link(b)

        self.assertEqual(out + ['ha!'], b.output())

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


class TestBaseCmd_link(TestBase):

    '''Link tests'''

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
        self.CMD.attr.req_type = [
            [('-f', ), ('.txt', )],
        ]
        a = self.sample()
        b = self.sample()
        c = self.sample()
        d = a.link(b).link(c)

        self.assertEqual(b.output, c.input)
        self.assertEqual(d, c)

    def test_linked_input_is_matched_to_args_without_error(self):
        self.CMD.attr.req_kwargs = ['-f']
        self.CMD.attr.req_type = [
            [('-f', ), ('.txt', )],
        ]

        a = self.CMD()
        a.output = lambda: ['file.txt']
        b = self.CMD()
        a.link(b)

        b._match_input_with_args()  # should not raise
        self.assertEqual(b.kwargs['-f'], a.output()[0])

    def test_cmd_raises_TypeError_if_req_not_fulfiled_by_link(self):
        self.CMD.attr.req_kwargs = ['-f']
        self.CMD.attr.req_type = [
            [('-f', ), ('.csv', )],
        ]

        a = self.CMD()
        b = self.CMD()
        a.link(b)

        with self.assertRaises(TypeError):
            b.cmd()

    def test_cmd_raises_TypeError_if_unable_to_use_linked_input(self):

        a = self.CMD()
        a.output = lambda: ['file.txt']
        b = self.CMD()
        a.link(b)

        with self.assertRaises(TypeError):
            b.cmd()

    def test_cmd_warns_about_unused_linked_input_if_not_strict(self):
        logger = logging.getLogger('libpipe.cmds.base')
        with patch.object(logger, 'warning') as mock_warn:
            a, b = self.linked_samples()
            b.strict = False
            b.cmd()  # should not raise
        self.assertTrue(mock_warn.called)

    def test_cmd_raises_TypeError_if_linked_input_unused_w_exact_match(self):
        self.CMD.attr.req_kwargs = ['-1']
        self.CMD.attr.req_type = [
            [('-1', '-2'), ('.txt', ), True],
        ]

        a = self.CMD()
        b = self.CMD()
        a.link(b)

        with self.assertRaises(TypeError):
            b.cmd()

    def test_linked_input_matches_correct_args_if_match_only_set(self):
        self.CMD.attr.req_type = [
            [('-1', '-2'), ('.txt', ), True],
            [('-U', ), ('.txt', ), True],
        ]

        a = self.CMD()
        b = self.CMD()
        c = self.CMD()

        a.output = lambda: ['file1.txt']
        b.output = lambda: ['file1.txt', 'file2.txt']
        a.link(b).link(c)

        expected = a.output()[0]
        b.cmd()  # should not raise!
        self.assertEqual(b.kwargs['-U'], expected)

        c.cmd()  # should not raise!
        self.assertEqual(c.kwargs['-1'], expected)

    def test_cmd_sets_unused_linked_input_to_req_args_if_not_strict(self):
        self.CMD.attr.req_args = 1
        self.CMD.attr.defaults = {}
        self.CMD.attr.req_type = [
            [('-1', '-2'), ('.txt', ), False],
        ]

        a = self.CMD()
        b = self.CMD(strict=False)
        a.output = lambda: ['file.txt', 'file.1', 'file2.txt']
        a.link(b)

        b.cmd()  # should not raise
        self.assertEqual(b.kwargs, {'-1': 'file.txt', '-2': 'file2.txt'})
        self.assertIn('file.1', b.args)

    def test_cmd_sets_linked_input_to_pos_arg_if_req_type_set(self):
        self.CMD.attr.req_args = 0
        self.CMD.attr.defaults = {}
        self.CMD.attr.req_type = [
            [(0, ), ('.txt', )],
        ]

        a = self.CMD()
        b = self.CMD()
        a.output = lambda: ['file.txt', 'file.1', 'file2.txt']
        a.link(b)

        b.cmd()  # should not raise
        self.assertEqual(b.args, ['file.txt'])


class TestBaseCmd_cmd(TestBase):

    '''Cmd tests'''

    def test_cmd_raises_TypeError_if_input_not_callable(self):
        a = self.CMD()
        a.input = ['file.txt']

        with self.assertRaises(TypeError):
            a.cmd()

    def test_cmd_escapes_all_unsafe_char(self):
        kwargs = {'-f': '&& rm * &&'}
        cmd = self.CMD(**kwargs)

        c = cmd.cmd(verbose=False, strict=False)
        expected = '{} -f \&\& rm \* \&\&'.format(self.ATTR.invoke_str)
        self.assertEqual(c, expected)

    def test_cmd_strict_removes_all_unsafe_char(self):
        kwargs = {'-f': '&& rm * &&'}
        cmd = self.CMD(**kwargs)

        c = cmd.cmd(verbose=False, strict=True)
        expected = '{} -f rm'.format(cmd.attr.invoke_str)
        self.assertEqual(c, expected)

    def test_cmd_retains_redirect_after_str_protect(self):
        cmd = self.CMD()
        cmd.redirect = ('>', 'out.txt')

        c = cmd.cmd(verbose=False, strict=True)
        expected = '{} > out.txt'.format(self.ATTR.invoke_str)
        self.assertEqual(c, expected)

    def test_cmd_is_wrapped_in_bash_assignment(self):
        a = self.CMD(wrap='foo')
        cmd_str = a.cmd()

        self.assertTrue(cmd_str.startswith('foo="$('))
        self.assertTrue(cmd_str.endswith(')"'))

    def test_cmd_puts_priority_args_before_kwargs(self):
        self.CMD.attr.n_priority_args = 2
        args = list('abc')
        kwargs = {'-x': '1', }
        cmd = self.CMD(*args, **kwargs)

        cmd_str = cmd.cmd(verbose=False)
        expected_str = '{} a b -x 1 c'.format(cmd.invoke_str)
        self.assertEqual(cmd_str, expected_str)

    def test_cmd_priority_args_neg_one_puts_all_args_first(self):
        self.CMD.attr.n_priority_args = -1
        args = list('abc')
        kwargs = {'-x': '1', }
        cmd = self.CMD(*args, **kwargs)

        cmd_str = cmd.cmd(verbose=False)
        expected_str = '{} a b c -x 1'.format(cmd.invoke_str)
        self.assertEqual(cmd_str, expected_str)
