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
