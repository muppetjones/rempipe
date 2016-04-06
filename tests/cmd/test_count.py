import os.path
import unittest
from unittest import mock

from libpipe.cmd import align
from libpipe.cmd import count
from libpipe.cmd import dummy

import logging
log = logging.getLogger(__name__)


class CountTestCase(unittest.TestCase):

    # TODO(sjbush): This is duplicated in several tests--refactor

    def get_cmd(self, *args, _input=None, **kwargs):
        '''Initialize the basic command'''
        if not _input:
            _input = self.default_input
        self.dummy = dummy.CmdDummy(*_input)
        cmd = self.CMD(*args, **kwargs)
        return self.dummy.link(cmd)


class TestHtseqCountCmd(CountTestCase):

    '''Test htseq-count

    NOTE: The first couple of tests ensure that CmdBase is working as
        expected.
    '''

    def setUp(self):
        patcher = mock.patch.object(align.Hisat2Index, '_check_extns')
        patcher.start()
        self.addCleanup(patcher.stop)

        index = align.Hisat2Index('genome/index')
        self.default_input = ['data/sample.s.bam', index]
        self.CMD = count.HtseqCountCmd

    def test_init_sets_order_arg_r_to_pos_by_defaults(self):
        '''Test that the default alignment order is pos'''
        cmd = self.get_cmd()
        self.assertIn('-r', cmd.kwargs)
        self.assertEqual(cmd.kwargs['-r'], 'pos')

    def test_init_sets_strandedness_arg_r_to_no_by_defaults(self):
        '''Test that strandedness is off by default'''
        cmd = self.get_cmd()
        self.assertIn('-s', cmd.kwargs)
        self.assertEqual(cmd.kwargs['-s'], 'no')

    def test_init_accepts_genome_kwarg_as_annotation_file(self):
        '''Test init(genome=val) sets val as gff file'''
        cmd = self.get_cmd(**{'genome': 'genome/hisat_index'})
        self.assertIn('genome/hisat_index', cmd.args)

    def test_pre_req_swaps_args_if_bam_or_sam_not_listed_first(self):
        '''Test that args are swapped iff BAM or SAM not first'''
        cmd = self.get_cmd()

        cmd.args = ['genome/genome.gtf', 'data/sample.bam']
        expected = list(cmd.args)
        expected.reverse()
        cmd._pre_req()

        self.assertEqual(cmd.args, expected)

    def test_pre_req_does_not_reverse_args_if_bam_or_sam_listed_first(self):
        cmd = self.get_cmd()

        cmd.args = ['data/sample.bam', 'genome/genome.gtf']
        expected = list(cmd.args)
        cmd._pre_req()

        self.assertEqual(cmd.args, expected)

    def test_pre_req_automatically_sets_f_flag_with_lc_input_extn(self):

        cmd = self.get_cmd()
        cmd.args = ['data/random.TXT', ]
        cmd._pre_req()

        self.assertIn('-f', cmd.kwargs)
        self.assertEqual(cmd.kwargs['-f'], 'txt')

    def test_pre_req_adds_gff_extn_to_arg_if_missing(self):
        '''Test 2nd arg given gff extn if missing (in case of hisat index)'''

        cmd = self.get_cmd()
        cmd._match_input_with_args()

        cmd._pre_req()
        self.assertEqual(cmd.args[1], self.default_input[1] + '.gtf')

    def test_pre_cmd_sets_redirect_to_align_file_w_count_extn(self):
        '''Test redirect points to a count file'''

        cmd = self.get_cmd()
        cmd._match_input_with_args()
        cmd._pre_req()
        cmd._pre_cmd()

        count_file = os.path.splitext(cmd.args[0])[0] + '.count'
        expected = ('>', count_file)
        self.assertEqual(cmd.redirect, expected)

    def test_output_contains_only_redirect_file(self):

        cmd = self.get_cmd()
        cmd._match_input_with_args()
        cmd.cmd()

        count_file = os.path.splitext(cmd.args[0])[0] + '.count'
        expected = [count_file]
        self.assertEqual(cmd.output(), expected)
