from unittest.mock import patch

from remsci.lib.tests.base import BaseTestCase

from libpipe.cmds.count import HtseqCountCmd

import logging
log = logging.getLogger(__name__)


class Test_CountCmd(BaseTestCase):

    def setUp(self):
        self.mock_isfile = self.mock_class('os.path.isfile')
        self.mock_isfile.return_value = True

    def setupCmd(self, **extra_kwargs):
        args = ['test.bam', '~/genome/test.gtf']
        kwargs = {}
        kwargs.update(extra_kwargs)

        return HtseqCountCmd(*args, **kwargs)

    def test_cmd_opens_config_file_in_genome_dir(self):
        m = self.setup_mock_read('-t gene')
        cmd = self.setupCmd()

        cmd.cmd()
        m.assert_called_once_with('~/genome/htseq-count.config', 'r')

    def test_cmd_includes_arguments_from_config_file_in_genome_dir(self):
        self.setup_mock_read('-t gene')
        cmd = self.setupCmd()
        cmd_str = cmd.cmd()

        self.assertIn('-t gene', cmd_str)

    def test_cmd_does_not_complain_if_no_config_file_in_genome_dir(self):
        m = self.setup_mock_read('-t gene')
        m.side_effect = FileNotFoundError()
        cmd = self.setupCmd()

        cmd.cmd()  # should not raise anything!!

    def test_cmd_raises_AttributeError_if_unknown_args_in_config_file(self):
        self.setup_mock_read('--unknown bad_value')
        cmd = self.setupCmd()

        with self.assertRaises(AttributeError):
            cmd.cmd()

    def test_cmd_raises_ValueError_if_unsafe_char_in_config_file(self):
        self.setup_mock_read('&& rm * &&')
        cmd = self.setupCmd()

        with self.assertRaises(ValueError):
            cmd.cmd(verbose=False)

    def test_cmd_raises_TypeError_if_cannot_detect_annotation(self):
        self.mock_isfile.return_value = False  # no file found
        cmd = self.setupCmd()
        cmd.args[1] = '~/genome/test'  # no extension

        with self.assertRaises(TypeError):
            cmd.cmd()

        # is file should be called more than once
        # -- may not be if extension list changed!
        self.assertGreater(self.mock_isfile.call_count, 1)

    def test_cmd_does_not_check_isfile_if_extn_given(self):
        self.mock_isfile.return_value = False  # no file found
        cmd = self.setupCmd()
        cmd.args[1] = '~/genome/test.gtf'  # no extension

        cmd.cmd()  # should not raise

        # isfile not called
        self.assertEqual(self.mock_isfile.call_count, 0)
