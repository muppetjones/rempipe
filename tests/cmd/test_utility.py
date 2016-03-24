import os.path
import unittest
from unittest import mock

import libpipe
from libpipe.cmd.dummy import CmdDummy
from libpipe.cmd.utility import SamtoolsSortCmd

import logging
log = logging.getLogger(__name__)


class UtilityTestCase(unittest.TestCase):

    def get_cmd(self, _input=None):
        '''Initialize the basic command'''
        if not _input:
            _input = self.default_input
        self.dummy = CmdDummy(*_input)
        cmd = self.CMD()
        return self.dummy.link(cmd)


class TestSamtoolsSortCmd(UtilityTestCase):

    '''Test samtools sort

    NOTE: The first couple of tests ensure that CmdBase is working as
        expected.
    '''

    def setUp(self):
        self.default_input = [
            'data/sample.sam', 'data/sample_unal.fq', 'genome/hisat_index']
        self.CMD = SamtoolsSortCmd

    def test_cmd_raises_IndexError_if_not_given_bam_or_sam_parg(self):
        cmd = self.get_cmd(_input=['not/a/bam', ])

        with mock.patch.object(libpipe.cmd.base.log, 'warning'):
            with self.assertRaises(IndexError):
                cmd.cmd()

    def test_match_sets_sam_or_bam_to_first_parg(self):
        extns = ['.bam', '.sam']

        for extn in extns:
            _input = ['data/sample' + extn]
            with self.subTest(extn=extn):
                cmd = self.get_cmd(_input=_input)
                cmd._match_input_with_args()
                self.assertEqual(cmd.args, _input)

    def test_cmd_uses_given_output(self):
        '''Test that output file is not overwritten if given'''

        cmd = self.get_cmd()

        expected = 'data/sorted.sam'
        cmd.kwargs['-o'] = expected
        cmd.cmd()

        self.assertEqual(cmd.kwargs['-o'], expected)

    def test_pre_cmd_ensures_o_value_has_bam_extension_by_default(self):
        '''Test that cmd defaults output name to <basename>.s.bam'''
        cmd = self.get_cmd()
        cmd._match_input_with_args()
        cmd._pre_cmd()

        self.assertIn('-o', cmd.kwargs)
        self.assertTrue(
            cmd.kwargs['-o'].endswith('.s.bam'),
            'Output not set with BAM extension',
        )

    def test_pre_cmd_sets_T_option_if_o_is_set(self):
        '''Ensure that -T is set with .tmp extn if -o '''
        cmd = self.get_cmd()
        cmd._match_input_with_args()
        cmd._pre_cmd()

        self.assertIn('-o', cmd.kwargs)
        self.assertIn('-T', cmd.kwargs)

        self.assertTrue(
            cmd.kwargs['-T'].endswith('.tmp'),
            'Temporary flag not set with .tmp extension',
        )

    def test_output_contains_only_o_flag(self):
        cmd = self.get_cmd()
        cmd.cmd()

        expected = [self.default_input[0].replace('.sam', '.s.bam')]
        self.assertEqual(cmd.output(), expected)
