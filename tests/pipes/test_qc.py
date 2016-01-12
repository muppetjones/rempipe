import os.path
import unittest
from unittest.mock import patch, Mock

from libpipe.cmds.base import BaseCmd
from libpipe.pipes.base import BasePipe
from libpipe.pipes.qc import TrimPipe


import logging
log = logging.getLogger(__name__)


class TestTrimPipe(unittest.TestCase):

    def test_job_name_defaults_to_TrimPipe(self):
        cmd = TrimPipe()
        self.assertEqual(cmd.job_name, 'TrimPipe')

    def test_init_calls_setup(self):
        with patch.object(TrimPipe, '_setup') as m:
            TrimPipe()
        m.assert_called_once_with(odir='', genome='', input_list=[])

    def test_output_is_fastq(self):
        pipe = TrimPipe(job_name="ptest", input_list=['ugh.fq'], odir='./foo')
        pipe.cmd()
        self.assertIn(
            os.path.join('./foo', 'ptest') + '-trimmed.fastq',
            pipe.output(),
        )

    def test_two_fq_in_means_two_fq_out(self):
        pipe = TrimPipe(
            job_name="ptest", input_list=['seq1.fq', 'seq2.fq'], odir='./foo')
        pipe.cmd()

        out_fq = [f for f in pipe.output() if f.endswith('.fastq')]
        self.assertEqual(len(out_fq), 2)
