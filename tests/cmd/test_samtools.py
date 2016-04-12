import unittest
from unittest import mock

import libpipe
from libpipe.cmd import align
from libpipe.cmd import dummy
from libpipe.cmd import samtools
from libpipe.type import index

import logging
log = logging.getLogger(__name__)


class SamtoolsTestCase(unittest.TestCase):

    # TODO(sjbush): This is duplicated in several tests--refactor

    def get_cmd(self, *args, _input=None, **kwargs):
        '''Initialize the basic command'''
        if not _input:
            _input = self.default_input
        self.dummy = dummy.CmdDummy(*_input)
        cmd = self.CMD(*args, **kwargs)
        return self.dummy.link(cmd)


class TestSamtoolsIndexCmd(SamtoolsTestCase):

    '''Test samtools index

    NOTE: The first couple of tests ensure that CmdBase is working as
        expected.
    '''

    def setUp(self):
        self.default_input = ['data/sample.s.bam', ]
        self.CMD = samtools.SamtoolsIndexCmd

    def test_cmd_raises_IndexError_if_not_given_bam_parg(self):
        cmd = self.get_cmd(_input=['not/a/bam', ])

        with mock.patch.object(libpipe.cmd.base.log, 'warning'):
            with self.assertRaises(IndexError):
                cmd.cmd()

    def test_match_sets_sam_or_bam_to_first_parg(self):
        cmd = self.get_cmd()
        cmd._match_input_with_args()
        self.assertEqual(cmd.args, self.default_input)

    def test_b_flag_set_by_default(self):
        '''Test that a BAI-format index is created by default'''
        cmd = self.get_cmd()
        self.assertIn('-b', cmd.flags)

    def test_fmt_kwarg_does_not_overwrite_fmt_passed_by_arg(self):
        cmd = self.get_cmd('-c')
        self.assertIn('-c', cmd.flags)
        self.assertNotIn('-b', cmd.flags)

    def test_output_contains_only_given_bam_file(self):
        '''Test that the output includes the given bam file, but not index'''

        cmd = self.get_cmd()
        cmd.cmd()

        self.assertEqual(cmd.output(), cmd.args)

    def test_output_includes_IndexType_from_input_if_any(self):
        with mock.patch.object(align.Hisat2Index, '_check_extns'):
            idx = align.Hisat2Index('path/to/index')
        cmd = self.get_cmd(_input=[idx] + self.default_input)
        cmd.cmd()
        self.assertIn(idx, cmd.output())


class TestSamtoolsSortCmd(SamtoolsTestCase):

    '''Test samtools sort

    NOTE: The first couple of tests ensure that CmdBase is working as
        expected.
    '''

    def setUp(self):
        # mimic expected hisat2 output
        self.default_input = [
            'data/sample.sam', 'data/sample_unal.fq', 'genome/hisat_index']
        self.CMD = samtools.SamtoolsSortCmd

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

    def test_output_contains_o_flag(self):
        cmd = self.get_cmd()
        cmd.cmd()

        expected = [self.default_input[0].replace('.sam', '.s.bam')]
        m = mock.MagicMock(side_effect=ValueError)  # no index matching
        with mock.patch.object(align.Hisat2Index, '_check_extns', m):
            self.assertEqual(cmd.output(), expected)

    def test_output_includes_IndexType_from_input_if_any(self):
        # NOTE: The same test from SamtoolsIndexCmd did NOT pass...
        #   Why not?
        cmd = self.get_cmd(_input=self.default_input)
        cmd.cmd()
        with mock.patch.object(libpipe.util.path, 'walk_safe') as m:
            m.return_value = {
                'file': ['hisat_index{}.ht2'.format(i) for i in range(8)]}
            self.assertIn('genome/hisat_index', cmd.output())
