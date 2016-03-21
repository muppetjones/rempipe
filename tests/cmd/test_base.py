
import unittest
from unittest import mock

import libpipe
from libpipe.cmd.attr import CmdAttributes
from libpipe.cmd.base import CmdInterface, CmdBase

# NOTE: flake8 complains about use of '_' for unused variables
#       e.g., _ = Obj()


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
        patcher = mock.patch.object(libpipe.cmd.base.log, 'error')
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

        kwargs = {'--unknown': 'bad_kwarg'}
        with self.assertRaises(ValueError):
            self.CMD(complain=True, **kwargs)


class TestCmdBase_cmd(BaseTestCase):

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

    def test_cmd_does_not_raise_if_no_AND_req_members_given(self):
        self.CMD.attr.req_kwargs = [('-f', '-o', '-n'), ]
        self.CMD.attr.defaults = {}
        kw = {'--foo': 'something different', }

        cmd = self.CMD(**kw)  # should not raise
        cmd.cmd()  # should NOT raise!!

    def test_cmd_raises_KeyError_if_do_not_have_all_AND_kwargs(self):
        self.CMD.attr.req_kwargs = [('-f', '-o', '-n')]  # req all three
        self.CMD.attr.defaults = {}
        kw = {'-f': 0, '-n': 2}

        cmd = self.CMD(**kw)  # should not raise
        with self.assertRaises(KeyError):
            cmd.cmd()

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
        self.CMD.attr.req_type = [
            [(0, 1, '-f'), ('.txt', )],
        ]
        args = ['a_text_file.txt', ]
        kwargs = {'-f': 'more_text.txt', }
        cmd = self.CMD(*args, **kwargs)

        cmd.cmd()  # should not raise; implicit check against missing [1]

    def test_cmd_raises_ValueError_if_wrong_kw_filetype_given(self):
        self.CMD.attr.req_type = [
            [('-f', ), ('.42', )],
        ]

        kwargs = {'-f': 'unexpected.txt'}
        cmd = self.CMD(**kwargs)
        with self.assertRaises(ValueError):
            cmd.cmd()
