import os.path
import unittest
from unittest.mock import patch

import libpipe
from libpipe.cmds.align import HisatCmd

import logging
log = logging.getLogger(__name__)


class TestHistatCmd(unittest.TestCase):

    def setUp(self):
        # prevent error logs from occuring during testing
        patcher = patch.object(libpipe.cmds.base.log, 'error')
        patcher.start()
        self.addCleanup(patcher.stop)

        # override base cmd method
        patcher = patch.object(libpipe.cmds.base.BaseCmd, '__cmd__')
        patcher.start()
        self.addCleanup(patcher.stop)

    def sample_cmd(self):
        kw = {
            '-U': 'path/seq.fa',
            '-x': 'gpath/gen',
            'timestamp': '000',
            '-S': 'path/al.sam',
        }
        return HisatCmd(**kw)

    def test_prepcmd_sets_redirect_to_log_file(self):
        hc = self.sample_cmd()
        hc.__prepcmd__()

        self.assertTrue(
            hc.redirect.endswith('path/seq_gen_000_hisat.log'),
            'Redirect not set to expected log file',
        )

    def test_prepcmd_sets_redirect_for_stdout_and_stderr_to_tee(self):
        hc = self.sample_cmd()
        hc.__prepcmd__()

        self.assertTrue(
            hc.redirect.startswith('2>&1 | tee -a'),
            'Redirect not set properly: {}'.format(hc.redirect),
        )

    def test_prepcmd_sets_unal_based_on_given_samfile_name(self):
        hc = self.sample_cmd()
        hc.__prepcmd__()

        expected_file = os.path.splitext(hc.kwargs['-S'])[0] + '.unal.fastq'

        self.assertIn('--un', hc.kwargs)
        self.assertEqual(hc.kwargs['--un'], expected_file)

    def test_check_requirements_raises_ValueError_if_only_one_ppe_given(self):
        hc = self.sample_cmd()
        hc.kwargs['-1'] = hc.kwargs['-U']
        del hc.kwargs['-U']
        with self.assertRaises(ValueError):
            hc._check_requirements()
