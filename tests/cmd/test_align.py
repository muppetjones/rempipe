import os.path
import unittest
from unittest import mock

import libpipe
from libpipe.cmd.align import Hisat2Cmd
from libpipe.cmd.dummy import CmdDummy

import logging
log = logging.getLogger(__name__)


DEFAULT_INPUT = ['data/sample.1.fq', 'data/sample.2.fq', 'genome/hisat_index']
TEST_CMD = Hisat2Cmd


class TestHistatCmd(unittest.TestCase):

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

    #
    #   Test match
    #

    def test_match_sets_kwarg_1_and_2_if_two_input_fastq(self):
        cmd = self.get_cmd()

        # NOTE: _match_input_with_args does NOT ensure only -1 or -U are set
        cmd._match_input_with_args()
        self.assertEqual(cmd.kwargs['-1'], DEFAULT_INPUT[0])
        self.assertEqual(cmd.kwargs['-2'], DEFAULT_INPUT[1])

    def test_match_checks_each_unused_input_as_possible_bowtie_index(self):
        '''Test unused input for bowtie index name (check for files)'''

        with mock.patch('libpipe.cmd.base.log.warning'):
            with mock.patch.object(
                    Hisat2Cmd,
                    '_check_for_index_files',
                    side_effect=ValueError
            ) as mock_check:
                cmd = self.get_cmd(_input=list('abc'))
                cmd._match_input_with_args()
            self.assertEqual(mock_check.call_count, 3)

    def test_match_sets_unused_input_as_bowtie_index_if_files_found(self):
        '''Test matching of -X (implicit check index success)'''

        cmd = self.get_cmd()
        cmd._match_input_with_args()  # walk_file_patched
        self.assertIn('-x', cmd.kwargs)
        self.assertEqual(cmd.kwargs['-x'], DEFAULT_INPUT[2])

    def test_check_index_raises_ValueError_if_bad_index(self):
        '''Test for exception if index not as expected (8 *.ht2 files)

        NOTE: if not *.ht2, won't be returned by walk_file (empty list)
        '''

        self.mock_walk.return_value = [
            'idx.{}.ht2'.format(i) for i in range(5)]
        with self.assertRaises(ValueError):
            Hisat2Cmd._check_for_index_files('genomes/index_name')

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

    #
    #   Cmd
    #

    # def test_cmd_raises_TypeError_if_1_and_U_set(self):
    #     hc = self.sample_cmd()
    #     hc.kwargs['-1'] = 'seq1.fq'
    #
    #     with self.assertRaises(AttributeError):
    #         hc.cmd()
