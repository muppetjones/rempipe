'''Test the CmdBase class

Test CmdBase directly (no instantiation) and indirectly (fake child).
Due to the size of the class, tests are split across several test cases.
Indirect tests are in alphabetical order, with the exception of init tests.

    * init
    * classmethods
    * cmd
    * link
    * requirements
'''

import subprocess
import unittest
from unittest import mock

import libpipe
from libpipe.cmd.attr import CmdAttributes
from libpipe.cmd.base import CmdBase

from tests.base import LibpipeTestCase  # includes read and write mock

import logging
log = logging.getLogger(__name__)

# NOTE: flake8 complains about use of '_' for unused variables
#       e.g., _ = Obj()
#       So...just remove.


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


class BaseTestCase(LibpipeTestCase):

    '''Setup a dummy base command for testing'''

    def setUp(self):
        ca = CmdAttributes(**sample_attributes)
        ca.req_args = 0
        ca.req_kwargs = []
        ca.req_types = []
        ca.defaults = {}

        class IndirectBase(CmdBase):
            attr = ca

            def output(self):
                return ['file.txt']
        self.CMD = IndirectBase
        self.ATTR = ca

        # prevent error logs from occuring during testing
        # -- any log.error calls **should** be expected!
        patcher = mock.patch.object(libpipe.cmd.base.log, 'error')
        patcher.start()
        self.addCleanup(patcher.stop)

    def mock_log_warning(self):
        # prevent error logs from occuring during testing
        # -- any log.error calls **should** be expected!
        patcher = mock.patch.object(libpipe.cmd.base.log, 'warning')
        patcher.start()
        self.addCleanup(patcher.stop)
        return


#-----------------------------------------------------------------------------
#   Direct Tests

class TestCmdBase(unittest.TestCase):

    def setUp(self):
        self.cmd_attr = CmdAttributes(**sample_attributes)

    def test_cannot_be_instatiated_directly(self):
        with self.assertRaises(TypeError):
            CmdBase()

    def test_children_must_define_output_and_attr(self):
        class FakeCmd(CmdBase):
            attr = mock.Mock(defaults={}, args=[])
            output = lambda x: None
        FakeCmd()  # should not raise!

    def test_children_raise_AttributeError_if_attr_not_set(self):
        class FakeCmd(CmdBase):
            output = lambda x: None
        with self.assertRaises(AttributeError):
            FakeCmd()  # should not raise!

    def test_name_property_returns_attr_name(self):
        FakeCmd = type('FakeCmd', (CmdBase,), dict(
            attr=self.cmd_attr, output=lambda: []))
        cmd = FakeCmd()
        self.assertEqual(cmd.name, cmd.attr.name)

    def test_invoke_property_returns_attr_invoke(self):
        FakeCmd = type('FakeCmd', (CmdBase,), dict(
            attr=self.cmd_attr, output=lambda: []))
        cmd = FakeCmd()
        self.assertEqual(cmd.invoke, cmd.attr.invoke)


#-----------------------------------------------------------------------------
#   Indirect Tests


class TestCmdBase_init(BaseTestCase):

    '''All tests related to command initialization'''

    def test_init_sets_defaults(self):
        '''Test that defaults defined by CmdAttributes are set properly'''
        self.CMD.attr.defaults = {'-n': 4800, }
        cmd = self.CMD()

        self.assertEqual(cmd.kwargs['-n'], self.CMD.attr.defaults['-n'])

    def test_init_defaults_overridden_if_kwargs_given(self):
        '''Test that given kwargs will override defaults

        NOTE: positional args and flags should NOT have defaults
        '''

        kwargs = {'-n': 8}
        self.CMD.attr.defaults['-n'] = 25
        cmd = self.CMD(**kwargs)

        self.assertEqual(cmd.kwargs['-n'], 8)  # set as expected
        self.assertEqual(self.CMD.attr.defaults['-n'], 25)  # unchanged

    def test_init_saves_timestamp_if_given(self):
        kw = {'timestamp': '151012-162900'}
        cmd = self.CMD(**kw)
        self.assertEqual(cmd.timestamp, kw['timestamp'])

    def test_init_saves_a_timestamp_if_not_given(self):
        cmd = self.CMD()
        self.assertRegex(cmd.timestamp, '\d{6}-\d{6}')

    def test_init_sets_redirect_to_None(self):
        cmd = self.CMD()
        self.assertIsNone(cmd.redirect)

    def test_init_sets_given_args_as_positional_args(self):
        args = list('abc')
        cmd = self.CMD(*args)
        self.assertEqual(cmd.args, args)

    def test_init_sets_hyphenated_input_args_as_flags(self):
        '''Test that hypenated *args are set as flags instead'''

        args = ['test.txt', '-v', '-x', 'test.xls']
        cmd = self.CMD(*args)

        expected_flags = [arg for arg in args if arg.startswith('-')]
        expected_args = [arg for arg in args if arg not in expected_flags]

        self.assertEqual(cmd.flags, expected_flags)
        self.assertEqual(cmd.args, expected_args)

    def test_init_sets_flags_given_as_kwarg_list(self):
        '''Test that flags given as a kwargs list are set after args flags'''

        args = ['test.txt', '-v', 'test.xls']
        cmd = self.CMD(*args, flags=['-x'])

        expected_flags = [arg for arg in args if arg.startswith('-')] + ['-x']
        self.assertEqual(cmd.flags, expected_flags)

    def test_init_ignores_unknown_flag_if_strict(self):
        '''Test that unknown kwargs are ignored by default'''

        cmd = self.CMD('--unknown')
        self.assertNotIn('--unknown', cmd.flags)

    def test_init_ignores_unknown_kwarg_if_strict(self):
        '''Test that unknown kwargs are ignored by default'''
        self.mock_log_warning()

        kwargs = {'--unknown': 'bad_kwarg'}
        cmd = self.CMD(**kwargs)
        self.assertNotIn('--unknown', cmd.kwargs)

    def test_init_sets_unknown_flag_if_not_strict(self):
        '''Test that unknown flags are set if not strict'''

        cmd = self.CMD('--unknown', strict=False)
        self.assertIn('--unknown', cmd.flags)

    def test_init_sets_unknown_kwarg_if_not_strict(self):
        '''Test that unknown kwargs are ignored by default'''

        kwargs = {'--unknown': 'bad_kwarg'}
        cmd = self.CMD(strict=False, **kwargs)
        self.assertIn('--unknown', cmd.kwargs)

    def test_init_raises_ValueError_if_unknown_kwarg_and_complain(self):
        '''Test init fails if unknown kwarg and complain=True'''
        self.mock_log_warning()

        kwargs = {'--unknown': 'bad_kwarg'}
        with self.assertRaises(ValueError):
            self.CMD(complain=True, **kwargs)

    def test_init_sets_odir_if_given(self):
        cmd = self.CMD(odir='la/la/la')
        self.assertEqual(cmd.output_dir, 'la/la/la')


class TestCmdBase_classmethods(BaseTestCase):

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


class TestCmdBase_cmd(BaseTestCase):

    '''Test cmd() output

    TODO: Protect against unsafe characters.
    TODO: Allow BASH variables in cmd string (removed if protected).
    TODO: Define 'priority_args' to allow some args to come before kwargs.

    NOTE: The above are required for some programs, e.g., velvet.
        These features are already implemented in rempipe, but are not
        included here.
    '''

    def test_magic_method_str_calls_non_verbose_cmd(self):
        cmd = self.CMD()

        with mock.patch.object(
                cmd, 'cmd', return_value="I've lost the bleeps") as mock_cmd:
            str(cmd)

        mock_cmd.assert_called_once_with(verbose=False)

    def test_cmd_calls_match_method(self):
        '''Test that input is matched to args from via cmd'''

        cmd = self.CMD()
        with mock.patch.object(cmd, '_match_input_with_args') as mock_match:
            cmd.cmd()
        mock_match.assert_called_once_with()

    def test_cmd_calls_update_output_directory_method(self):

        cmd = self.CMD()
        with mock.patch.object(cmd, '_update_output_directory') as mock_match:
            cmd.cmd()
        mock_match.assert_called_once_with()

    def test_cmd_calls_each_user_defined_customization(self):
        '''Test cmd calls each available user defined methods for cmd prep'''

        avail_user_custom = ['_pre_req', '_post_req', '_pre_cmd']

        cmd = self.CMD()
        mocked = {}
        for cust in avail_user_custom:
            mocked[cust] = mock.Mock()
            setattr(cmd, cust, mocked[cust])

        cmd.cmd()
        for cust in avail_user_custom:
            with self.subTest(custom_method=cust):
                mocked[cust].assert_called_once_with()

    def test_cmd_verbose_adds_bash_safe_line_breaks(self):
        '''Test that verbose adds '\\n between each arg'''

        kwargs = {'-f': 'req_kwarg', '-n': 'a'}
        cmd = self.CMD(**kwargs)
        cmd_str = cmd.cmd()

        line_break = '\\\n'
        n_breaks = cmd_str.count(line_break)

        # +1 line break for break btwn invoke and first arg
        self.assertEqual(n_breaks, len(cmd.kwargs))

    def test_cmd_returns_expected_kwargs_in_cmd_string(self):
        '''Test that 'cmd' pieces together given kwargs in an expected manner

        The command should be "[invoke] [kwargs] [args] > [redirect]".
        '''

        args = ['seq1.fq', 'seq2.fq', '-x']
        kwargs = {'-f': 'req_kwarg', '-n': 'a'}
        cmd = self.CMD(*args, **kwargs)

        expected_kwargs = ' '.join([
            '{} {}'.format(k, v)
            for k, v in sorted(kwargs.items())
        ])
        expected_args = ' '.join(['seq1.fq', 'seq2.fq'])
        expected_flags = ' '.join(['-x'])
        expected_cmd = ' '.join(filter(None, [
            self.CMD.attr.invoke,
            expected_flags,
            expected_kwargs,
            expected_args,
        ]))

        self.assertEqual(
            cmd.cmd(verbose=False).rstrip(), expected_cmd.rstrip())

    def test_cmd_uses_flag_sep_variable(self):

        self.CMD.attr.flag_sep = '='
        kwargs = {'-x': 'bar', }
        cmd = self.CMD(**kwargs)

        cmd_str = cmd.cmd(verbose=False)
        expected_str = '{} -x=bar'.format(cmd.attr.invoke)
        self.assertEqual(cmd_str, expected_str)

    def test_cmd_does_not_duplicate_given_args_after_each_run(self):
        '''Test that match is not re-adding args'''

        args = ['a_text_file.txt', ]
        kwargs = {'-f': 'more_text.txt', }
        cmd = self.CMD(*args, **kwargs)

        for i in range(100):
            cmd.cmd()
            self.assertEqual(cmd.args, args)
            self.assertEqual(cmd.kwargs, cmd.kwargs)

    def test_cmd_does_not_duplicate_input_args_after_each_run(self):
        '''Test that match is not re-adding args'''

        self.CMD.attr.req_kwargs = ['-f']
        self.CMD.attr.req_types = [
            [(0, '-f', ), ('.txt', )],
        ]

        cmd = self.CMD()
        cmd.input = lambda: ['fake', 'file0.txt', 'file-f.txt']
        args = ['file0.txt']
        kwargs = {'-f': 'file-f.txt'}

        for i in range(100):
            cmd.cmd()
            self.assertEqual(cmd.args, args)
            self.assertEqual(cmd.kwargs, kwargs)

    def test_cmd_str_includes_pargs(self):
        '''Test that a single arg is included (from bug)'''

        # This is a known bug.
        self.CMD.attr.req_args = 1
        self.CMD.attr.req_kwargs = {}
        self.CMD.attr.req_types = [
            [(0, 1), ('.bam', )]
        ]
        cmd = self.CMD()
        cmd.input = lambda: ['data/seq.s.bam', 'data/seq2.s.bam']

        cmd_str = cmd.cmd()
        self.assertIn('data/seq.s.bam', cmd_str)
        self.assertIn('data/seq2.s.bam', cmd_str)

    def test_cmd_updates_args_in_attr_output_args_with_given_odir(self):
        '''Test the dir path is updated on listed args'''

        self.CMD.attr.output_args = [1, '--blue']
        cmd = self.CMD(odir='brave/new')
        cmd.args = ['hello', 'path/to/world']
        cmd.kwargs = {'--blue': 'suede/shoes', '-o': 'suede/sandals'}

        cmd.cmd()
        expected_args = ['hello', 'brave/new/world']
        expected_kwargs = {'--blue': 'brave/new/shoes', '-o': 'suede/sandals'}

        self.assertEqual(cmd.args, expected_args)
        self.assertEqual(cmd.kwargs, expected_kwargs)

    def test_cmd_updates_redirect_tuple_with_odir(self):
        '''Redirect should be handled by the child, but child may not know'''

        cmd = self.CMD(odir='brave/new')
        cmd.redirect = ('>', 'hello/world')
        expected_redirect = ('>', 'brave/new/world')

        cmd.cmd()
        self.assertEqual(cmd.redirect, expected_redirect)

    def test_cmd_does_not_update_redirect_str_with_odir(self):
        '''Redirect value not modified if str'''

        cmd = self.CMD(odir='brave/new')
        cmd.redirect = '> hello/world'
        expected_redirect = '> hello/world'

        cmd.cmd()
        self.assertEqual(cmd.redirect, expected_redirect)

    def test_update_output_dir_does_not_raise_if_output_arg_not_set(self):
        self.CMD.attr.output_args = [1, '--blue']
        cmd = self.CMD(odir='brave/new')
        cmd._update_output_directory()  # should not raise

    def test_update_output_directory_does_not_raise_if_no_odir_given(self):

        self.CMD.attr.output_args = [1, '--blue']
        cmd = self.CMD()
        expected_args = ['hello', 'path/to/world']
        expected_kwargs = {'--blue': 'suede/shoes', '-o': 'suede/sandals'}

        cmd.args = expected_args
        cmd.kwargs = expected_kwargs

        cmd._update_output_directory()  # should not raise
        self.assertEqual(cmd.args, expected_args)
        self.assertEqual(cmd.kwargs, expected_kwargs)


class TestBaseCmd_link(BaseTestCase):

    '''Test linking commands together

    NOTE: This group tests '_match_input_with_args', too.
    '''

    def test_link_sets_dest_input_to_src_output(self):
        a = self.CMD()
        b = self.CMD()
        a.link(b)

        self.assertEqual(a.output, b.input)

    def test_link_chaining(self):
        a = self.CMD()
        b = self.CMD()
        c = a.link(b)

        self.assertNotEqual(a, b)
        self.assertNotEqual(a, c)
        self.assertEqual(b, c)

    def test_match_sets_args_with_basic_types(self):
        self.CMD.attr.req_types = [
            [(0, ), (int, )],
        ]

        cmd = self.CMD()
        cmd.input = lambda: ['hello', 'world', 42, 0.1]

        cmd._match_input_with_args()
        self.assertEqual(cmd.args[0], 42)

    def test_match_sets_args_with_file_types(self):
        self.CMD.attr.req_types = [
            [('-f', ), ('.foo', )],
        ]

        cmd = self.CMD()
        cmd.input = lambda: ['hello', 'world.foo', 42, 0.1]

        cmd._match_input_with_args()
        self.assertEqual(cmd.kwargs['-f'], 'world.foo')

    def test_linked_input_is_matched_by_order_given_in_req_types(self):
        self.CMD.attr.req_kwargs = ['-f']
        self.CMD.attr.req_types = [
            [(0, '-f', ), ('.txt', )],
        ]

        cmd = self.CMD()
        cmd.input = lambda: ['fake', 'file-f.txt', 'file0.txt']

        cmd._match_input_with_args()  # should not raise
        self.assertEqual(cmd.kwargs['-f'], cmd.input()[2])
        self.assertEqual(cmd.args[0], cmd.input()[1])

    def test_match_raises_AttributeError_if_not_linked(self):
        cmd = self.CMD()
        with self.assertRaises(AttributeError):
            cmd._match_input_with_args()  # ignored by cmd!

    def test_match_raises_TypeError_if_input_not_callable(self):
        cmd = self.CMD()
        cmd.input = list('abc')
        with self.assertRaises(TypeError):
            cmd.cmd()  # raised by match, caught and raised by cmd

    def test_match_warns_if_input_not_used_and_strict(self):
        logger = logging.getLogger('libpipe.cmd.base')
        with mock.patch.object(logger, 'warning') as mock_warn:
            cmd = self.CMD()
            cmd.input = lambda: ['file.txt']
            cmd._match_input_with_args()
        self.assertTrue(mock_warn.called)

    def test_match_raises_ValueError_if_input_not_used_and_complaining(self):
        self.mock_log_warning()  # already tested for warning
        cmd = self.CMD(complain=True)
        cmd.input = lambda: ['file.txt']
        with self.assertRaises(ValueError):
            cmd.cmd()  # raised by match, caught and raised by cmd

    def test_match_exact_sets_args_with_perfect_match_1(self):
        '''Test for exact match single file to '-U' arg'''
        self.CMD.attr.req_types = [
            [('-U', ), ('.txt', ), True],
            [('-1', '-2'), ('.txt', ), True],
        ]

        cmd = self.CMD(complain=True)
        cmd.input = lambda: ['fileU.txt', ]
        cmd._match_input_with_args()
        self.assertEqual(cmd.kwargs['-U'], 'fileU.txt')
        self.assertNotIn('-1', cmd.kwargs)
        self.assertNotIn('-2', cmd.kwargs)

    def test_match_exact_sets_args_with_perfect_match_2(self):
        '''Test for exact match two files to '-1' and '-2' args'''
        self.CMD.attr.req_types = [
            [('-U', ), ('.txt', ), True],
            [('-1', '-2'), ('.txt', ), True],
        ]

        cmd = self.CMD(complain=True)
        cmd.input = lambda: ['file1.txt', 'file2.txt']
        cmd._match_input_with_args()
        self.assertEqual(cmd.kwargs['-1'], 'file1.txt')
        self.assertEqual(cmd.kwargs['-2'], 'file2.txt')
        self.assertNotIn('-U', cmd.kwargs)


class TestCmdBase_requirements(BaseTestCase):

    #
    #   Positional requirements
    #

    def test_cmd_raises_IndexError_if_missing_req_number_args(self):
        self.CMD.attr.req_args = 2
        args = ['only_one_arg']
        cmd = self.CMD(*args)  # shold not raise

        with self.assertRaises(IndexError):
            cmd.cmd()

    #
    #   Keyword requirements
    #

    def test_cmd_raises_KeyError_if_required_kwarg_not_given(self):
        self.CMD.attr.req_kwargs = ['-f']
        cmd = self.CMD()  # should not raise

        with self.assertRaises(KeyError):
            cmd.cmd()

    def test_cmd_does_not_raise_if_all_AND_req_members_given(self):
        self.CMD.attr.req_kwargs = [('-f', '-o', '-n'), ]
        self.CMD.attr.defaults = {}
        kwargs = {'-f': 0, '-n': 2, '-o': 42}

        cmd = self.CMD(**kwargs)  # should not raise
        cmd.cmd()  # should NOT raise!!

    def test_cmd_does_not_raise_if_no_AND_req_members_given(self):
        self.CMD.attr.req_kwargs = [('-f', '-o', '-n'), ]
        self.CMD.attr.defaults = {}
        kwargs = {'--foo': 'something different', }

        cmd = self.CMD(**kwargs)  # should not raise
        cmd.cmd()  # should NOT raise!!

    def test_cmd_raises_KeyError_if_do_not_have_all_AND_kwargs(self):
        self.CMD.attr.req_kwargs = [('-f', '-o', '-n')]  # req all three
        self.CMD.attr.defaults = {}
        kwargs = {'-f': 0, '-n': 2}

        cmd = self.CMD(**kwargs)  # should not raise
        with self.assertRaises(KeyError):
            cmd.cmd()

    def test_cmd_does_not_raise_if_exactly_one_XOR_arg_given(self):
        self.CMD.attr.req_kwargs = [['-f', '-o', '-n'], ]
        self.CMD.attr.defaults = {}
        kwargs = {'-f': 0, }

        cmd = self.CMD(**kwargs)
        cmd.cmd()  # should not raise

    def test_cmd_raises_KeyError_if_missing_all_XOR_args(self):
        self.CMD.attr.req_kwargs = [['-f', '-o', '-n'], ]
        self.CMD.attr.defaults = {}

        cmd = self.CMD(**{'--foo': 1})
        with self.assertRaises(KeyError):
            cmd.cmd()

    def test_cmd_raises_KeyError_if_gt_1_req_XOR_arg_given(self):
        self.CMD.attr.req_kwargs = [['-f', '-o', '-n'], ]
        self.CMD.attr.defaults = {}
        kw = {'-f': 0, '-n': 2}

        cmd = self.CMD(**kw)
        with self.assertRaises(KeyError):
            cmd.cmd()

    #
    #   Type requirements
    #

    def test_cmd_type_checks_pos_and_kw_arguments(self):
        self.CMD.attr.req_types = [
            [(0, 1, '-f'), ('.txt', )],
        ]
        args = ['a_text_file.txt', ]
        kwargs = {'-f': 'more_text.txt', }
        cmd = self.CMD(*args, **kwargs)

        cmd.cmd()  # should not raise; implicit check against missing [1]

    def test_cmd_raises_ValueError_if_wrong_pos_filetype_given(self):
        '''Test for file type checking (implicitly check non-basic type)'''

        self.CMD.attr.req_types = [
            [(1, ), ('.42', )],
        ]

        args = ['not_correct_arg.42', 'unexpected.txt', 'no_swap']
        cmd = self.CMD(*args)

        with self.assertRaises(ValueError):
            cmd.cmd()

    def test_cmd_raises_ValueError_if_wrong_kw_filetype_given(self):
        '''Test for file type checking (implicitly check non-basic type)'''

        self.CMD.attr.req_types = [
            [('-f', ), ('.42', )],
        ]

        kwargs = {'-f': 'unexpected.txt'}
        cmd = self.CMD(**kwargs)
        with self.assertRaises(ValueError):
            cmd.cmd()

    def test_cmd_raises_ValueError_if_arg_has_wrong_type(self):
        for t, v in [(int, 0.5), (float, 'bar')]:
            with self.subTest(type=t):
                self.CMD.attr.req_types = [
                    [('-f', ), (t, )],
                ]

                kwargs = {'-f': v}
                cmd = self.CMD(**kwargs)

                with self.assertRaises(ValueError):
                    cmd.cmd()

    def test_cmd_checks_multiple_basic_types(self):
        self.CMD.attr.req_types = [
            [('-f', ), (int, float)],
        ]

        kwargs = {'-f': 0.01}
        cmd = self.CMD(**kwargs)

        cmd.cmd()  # should NOT raise--int will fail, but float should not

    def test_cmd_checks_mixed_types(self):
        '''[YAGNI?] Test that args can have basic and file type req'''
        self.CMD.attr.req_types = [
            [('-f', '-x', ), (int, '.txt')],
        ]

        kwargs = {'-f': 3, '-x': 'hi.txt'}
        cmd = self.CMD(**kwargs)

        cmd.cmd()  # should NOT raise

    def test_cmd_tries_reversing_args_on_req_file_fail_IFF_len_eq_2(self):
        '''Test that a simple swap (2 args only) won't break things'''
        self.CMD.attr.req_types = [
            [(0, ), ('.txt', )],
        ]
        args = ['not_file', 'file.txt']
        cmd = self.CMD(*args)
        cmd.cmd()  # should not raise
        self.assertNotEqual(cmd.args, args)  # reversed!


class TestCmdBase_run(BaseTestCase):

    def setUp(self):
        super().setUp()

        # prevent actually calling any commands
        patcher = mock.patch('subprocess.check_call')
        self.mock_call = patcher.start()
        self.addCleanup(patcher.stop)

    #
    #   Locally
    #

    def test_runs_locally_by_default(self):
        cmd = self.CMD()
        with mock.patch.object(cmd, '_run_local') as m:
            cmd.run()
        m.assert_called_once_with()

    def test_local_run_calls_subprocess_with_cmd_str(self):
        cmd = self.CMD()
        cmd.run()

        cmd_str = cmd.cmd()
        self.mock_call.assert_called_once_with(cmd_str)

    def test_local_opens_and_writes_to_logfile_if_given(self):
        mock_open = self.setup_mock_write()
        cmd = self.CMD()

        log_file = 'cmd.log'
        cmd_str = cmd.cmd()
        cmd.run(log_file=log_file)

        mock_open.assert_called_once_with(log_file, 'w')
        self.mock_call.assert_called_once_with(
            cmd_str, stdout=mock_open(), stderr=mock_open())

    def test_local_raises_CalledProcessError_if_problem(self):
        self.mock_call.side_effect = subprocess.CalledProcessError(1, 'AH!')

        cmd = self.CMD()
        with self.assertRaises(subprocess.CalledProcessError):
            cmd.run()

    #
    #   Resource management queue
    #

    def test_calls_run_queue_if_queue_given(self):
        cmd = self.CMD()
        with mock.patch.object(cmd, '_run_queue') as m:
            cmd.run(queue='qsub')
        m.assert_called_once_with(queue='qsub')

    def test_run_queue_raises_ValueError_if_unknown_queue(self):
        cmd = self.CMD()
        with self.assertRaises(ValueError):
            cmd.run(queue='qsub')
