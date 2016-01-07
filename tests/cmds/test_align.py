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
        patcher = patch.object(libpipe.cmds.base.BaseCmd, '_cmd')
        patcher.start()
        self.addCleanup(patcher.stop)

    def sample_cmd(self):
        kw = {
            '-U': 'upath/seq.fa',
            '-x': 'gpath/gen',
            'timestamp': '000',
            '-S': 'path/al.sam',
        }
        return HisatCmd(**kw)

    #
    #   Test _prepcmd
    #

    def test_prepcmd_sets_S_if_not_given(self):
        hc = self.sample_cmd()
        del hc.kwargs['-S']
        hc._prepcmd()

        self.assertEqual(
            hc.kwargs['-S'],
            'upath/seq_gen.sam',
        )

    def test_prepcmd_sets_redirect_to_log_file(self):
        hc = self.sample_cmd()
        hc._prepcmd()

        self.assertTrue(
            hc.redirect.endswith('path/al_gen_000_hisat.log'),
            'Redirect not set to expected log file ({})'.format(hc.redirect),
        )

    def test_prepcmd_sets_redirect_for_stdout_and_stderr_to_tee(self):
        hc = self.sample_cmd()
        hc._prepcmd()

        self.assertTrue(
            hc.redirect.startswith('2>&1 | tee -a'),
            'Redirect not set properly: {}'.format(hc.redirect),
        )

    def test_prepcmd_sets_unal_based_on_given_samfile_name(self):
        hc = self.sample_cmd()
        hc._prepcmd()

        expected_file = os.path.splitext(hc.kwargs['-S'])[0] + '.unal.fastq'

        self.assertIn('--un', hc.kwargs)
        self.assertEqual(hc.kwargs['--un'], expected_file)

    #
    #   Test cmd
    #

    def test_cmd_raises_AttributeError_if_only_one_ppe_given(self):
        hc = self.sample_cmd()
        hc.kwargs['-1'] = hc.kwargs['-U']
        del hc.kwargs['-U']
        with self.assertRaises(AttributeError):
            hc.cmd()

    def test_addreq_raises_FileNotFoundError_if_n_idx_ne_expected(self):

        with patch('remsci.lib.utility.path.walk_file') as m:
            for i in [0, 100]:
                with self.subTest(n_indx=i):
                    m.return_value = [0] * i
                    hc = self.sample_cmd()
                    with self.assertRaises(FileNotFoundError):
                        hc._additional_requirements()

    #
    #   Test _prepreq
    #

    def test_prepreq_raises_TypeError_if_linked_input_not_used(self):
        with patch.object(
                HisatCmd, 'output', autospec=True, return_value=['seq.txt']):
            ohc = self.sample_cmd()
            ihc = self.sample_cmd()
            ohc.link(ihc)

        with self.assertRaises(TypeError):
            ihc._prepreq()

    def test_prepreq_sets_single_link_input_to_U_kwarg(self):

        with patch.object(HisatCmd, 'output', return_value=['seq.fq']):
            ohc = self.sample_cmd()
            ihc = self.sample_cmd()
            ohc.link(ihc)
        ihc._prepreq()

        self.assertEqual(ihc.kwargs['-U'], 'seq.fq')

    def test_prepreq_sets_double_link_input_to_1_and_2_kwarg(self):

        args = ['seq.1.fq', 'seq.2.fq']
        with patch.object(HisatCmd, 'output', return_value=args):
            ohc = self.sample_cmd()
            ihc = self.sample_cmd()
            ohc.link(ihc)
        ihc._prepreq()

        self.assertEqual(ihc.kwargs['-1'], 'seq.1.fq')
        self.assertEqual(ihc.kwargs['-2'], 'seq.2.fq')

    def test_prepreq_preserves_kwargs_if_no_input_given(self):

        ihc = self.sample_cmd()
        ihc._prepreq()

        self.assertEqual(ihc.kwargs['-U'], 'upath/seq.fa')


if __name__ == '__main__':
    unittest.main()
