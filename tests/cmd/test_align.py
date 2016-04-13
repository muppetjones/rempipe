import os.path
import unittest
from unittest import mock

import libpipe
from libpipe.cmd import align
from libpipe.cmd.dummy import CmdDummy
from libpipe.type import index

from tests.type import test_index

import logging
log = logging.getLogger(__name__)


DEFAULT_INPUT = ['data/sample.1.fq', 'data/sample.2.fq', 'genome/hisat_index']
TEST_CMD = align.Hisat2Cmd


class TestHisatIndex(test_index.IndexTypeTestCase):

    def test_creates_index_obj(self):
        self.setup_mock_check_extns(align.Hisat2Index)
        obj = align.Hisat2Index('genome/index')
        self.assertIsInstance(obj, index.IndexType)

    def test_expects_8_ht2_files(self):
        cls = align.Hisat2Index
        self.assertEqual(cls.extns, ['.ht2'])
        self.assertEqual(cls.counts, [8])


class TestHisat2Cmd(unittest.TestCase):

    def setUp(self):

        # log test name (for debugging)
        # id_split = self.id().split('.')
        # log.debug('-' * 50 + '\n\t' + '.'.join(id_split[-2:]))

        # override walk_file method
        # -- used when checking for index
        patcher = mock.patch('libpipe.util.path.walk_file')
        self.mock_walk = patcher.start()
        self.mock_walk.return_value = [
            'idx.{}.ht2'.format(i) for i in range(8)]
        self.addCleanup(patcher.stop)

    def get_cmd(self, _input=None):
        '''Initialize the basic command'''
        if not _input:
            _input = DEFAULT_INPUT
        self.dummy = CmdDummy(*_input)
        cmd = TEST_CMD()
        return self.dummy.link(cmd)

    def test_command_can_be_initialized(self):
        TEST_CMD()  # should not fail

    #
    #   Test Output
    #

    def test_output_includes_exactly_one_samfile(self):

        cmd = self.get_cmd()
        cmd.cmd()  # must call to process args
        sam_files = [o for o in cmd.output() if o.lower().endswith('.sam')]
        self.assertEqual(len(sam_files), 1)

    def test_output_includes_two_unal_fq_for_pe_input(self):
        cmd = self.get_cmd()
        cmd.cmd()  # must call to process args

        prefix, extn = os.path.splitext(cmd.kwargs['-S'])
        expected = [
            '{}_unal.{}{}'.format(prefix, i + 1, '.fastq') for i in range(2)
        ]
        unal_files = [o for o in cmd.output() if 'unal' in o]
        self.assertEqual(unal_files, expected)

    def test_output_includes_one_unal_fq_for_se_input(self):
        cmd = self.get_cmd(_input=DEFAULT_INPUT[1:])
        cmd.cmd()  # must call to process args

        prefix, extn = os.path.splitext(cmd.kwargs['-S'])
        expected = ['{}_unal{}'.format(prefix, '.fastq')]
        unal_files = [o for o in cmd.output() if 'unal' in o]
        self.assertEqual(unal_files, expected)

    def test_output_includes_given_index(self):
        cmd = self.get_cmd()
        cmd.cmd()
        self.assertIn(cmd.kwargs['-x'], cmd.output())

    #
    #   Test match
    #

    def test_match_sets_kwarg_1_and_2_if_two_input_fastq(self):
        cmd = self.get_cmd()

        # NOTE: _match_input_with_args does NOT ensure only -1 or -U are set
        cmd._match_input_with_args()
        self.assertEqual(cmd.kwargs['-1'], DEFAULT_INPUT[0])
        self.assertEqual(cmd.kwargs['-2'], DEFAULT_INPUT[1])

    def test_match_sets_kwarg_x_to_index_type(self):
        '''Ensure custom matching to Hisat2Index works

        NOTE: Unused input checking now extraneous
        '''

        cmd = self.get_cmd()
        cmd._match_input_with_args()  # walk_file_patched
        self.assertIn('-x', cmd.kwargs)
        self.assertIsInstance(cmd.kwargs['-x'], align.Hisat2Index)
        self.assertEqual(cmd.kwargs['-x'], DEFAULT_INPUT[2])

    def test_check_req_raises_ValueError_if_bad_index(self):
        '''Test for exception if index not as expected (8 *.ht2 files)

        NOTE: if not *.ht2, won't be returned by walk_file (empty list)
        '''

        self.mock_walk.return_value = [
            'idx.{}.ht2'.format(i) for i in range(5)]
        cmd = align.Hisat2Cmd()
        cmd.kwargs.update({
            '-x': 'genomes/index_name',
            '-U': 'something.fq',
        })
        with self.assertRaises(ValueError):
            cmd.cmd()

    #
    #   Test _prepcmd
    #

    def test_pre_cmd_sets_S_based_on_input_if_S_not_given_paired_end(self):
        cmd = self.get_cmd()
        cmd._match_input_with_args()
        cmd._pre_cmd()

        self.assertIn('-S', cmd.kwargs)
        self.assertEqual(cmd.kwargs['-S'], 'data/sample.sam')

    def test_pre_cmd_sets_S_based_on_input_if_S_not_given_single_end(self):
        cmd = self.get_cmd(_input=DEFAULT_INPUT[1:])
        cmd._match_input_with_args()
        cmd._pre_cmd()

        self.assertIn('-S', cmd.kwargs)
        self.assertEqual(cmd.kwargs['-S'], 'data/sample.2.sam')

    def test_pre_cmd_removes_trailing_underscore_from_S_common_prefix(self):
        _input = ['seq_1.fq', 'seq_2.fq']
        cmd = self.get_cmd(_input=_input)
        cmd._match_input_with_args()
        cmd._pre_cmd()

        self.assertIn('-S', cmd.kwargs)
        self.assertEqual(cmd.kwargs['-S'], 'seq.sam')

    def test_pre_cmd_sets_redirect_to_log_file(self):
        cmd = self.get_cmd()
        cmd._match_input_with_args()
        cmd._pre_cmd()

        self.assertIsNotNone(cmd.redirect)
        found = cmd.redirect[-1]
        expected = 'data/hisat_sample.log'

        # NOTE: The log filename is created from the output file, so the
        #       genome name will only be in the file if the user puts it
        #       or if the user does not provide an output file.
        self.assertEqual(
            found, expected,
            'Redirect not set to expected log file ({})'.format(cmd.redirect),
        )

    def test_prepcmd_sets_unal_based_on_given_samfile_name_se(self):
        cmd = self.get_cmd(_input=DEFAULT_INPUT[1:])
        cmd._match_input_with_args()
        cmd._pre_cmd()

        expected_file = os.path.splitext(cmd.kwargs['-S'])[0] + '_unal.fastq'

        self.assertIn('--un', cmd.kwargs)
        self.assertEqual(cmd.kwargs['--un'], expected_file)

    def test_prepcmd_sets_unal_based_on_given_samfile_name_pe(self):
        cmd = self.get_cmd()
        cmd._match_input_with_args()
        cmd._pre_cmd()

        expected_file = os.path.splitext(cmd.kwargs['-S'])[0] + '_unal.fastq'

        self.assertIn('--un-conc', cmd.kwargs)
        self.assertEqual(cmd.kwargs['--un-conc'], expected_file)

    def test_pre_cmd_sets_redirect_for_stdout_and_stderr_to_tee(self):
        cmd = self.get_cmd()
        cmd._match_input_with_args()
        cmd._pre_cmd()

        self.assertIsNotNone(cmd.redirect)
        found = cmd.redirect[:-1]
        expected = ('2>&1', '|', 'tee -a')

        self.assertEqual(
            found, expected,
            'Redirect not set properly: {}'.format(cmd.redirect),
        )

    #
    #   Cmd
    #

    def test_cmd_raises_KeyError_if_1_and_U_set(self):
        cmd = self.get_cmd()
        # cmd._match_input_with_args()
        cmd.kwargs['-U'] = 'data/seqU.fq'

        with mock.patch.object(libpipe.cmd.base.log, 'error'):
            with self.assertRaises(KeyError):
                cmd.cmd()
