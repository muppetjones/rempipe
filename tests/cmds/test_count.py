from unittest.mock import patch

from remsci.lib.tests.base import BaseTestCase

from libpipe.cmds.count import HtseqCountCmd

import logging
log = logging.getLogger(__name__)


class Test_CountPipe(BaseTestCase):

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
