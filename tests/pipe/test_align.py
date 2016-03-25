
import os.path
from unittest import mock

from libpipe.pipe.align import AlignPipe
from tests.base import LibpipeTestCase  # includes read and write mock


import logging
log = logging.getLogger(__name__)


class TestAlignPipe(LibpipeTestCase):

    def setUp(self):

        self.used_cmds = [
            'align.Hisat2Cmd',
            'samtools.SamtoolsSortCmd',
            'samtools.SamtoolsIndexCmd',
            'count.HtseqCountCmd',
        ]
        short_names = [
            (os.path.splitext(cmd)[1])[1:-3].lower()
            for cmd in self.used_cmds
        ]

        # create mocks for [mod.NameCmd] -> {name: patch}
        patchers = [
            mock.patch('libpipe.cmd.{}'.format(cmd))
            for cmd in self.used_cmds
        ]
        self.mock_cmds = {
            short_name: patcher.start()
            for short_name, patcher in zip(short_names, patchers)
        }
        for patcher in patchers:
            self.addCleanup(patcher.stop)

        log.debug(self.mock_cmds)

    def test_init_creates_default_cmd_objects(self):
        '''Ensure hisat2, samtools sort & index, & htseq-count cmds created'''

        AlignPipe()
        for name, mock_cmd in self.mock_cmds.items():
            with self.subTest(cmd=name):
                self.assertEqual(
                    mock_cmd.call_count, 1, '{} not called'.format(name))
