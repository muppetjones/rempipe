'''Test the PipeBase class'''

import unittest
from unittest import mock


class PipeBaseTestCase(unittest.TestCase):

    def setUp(self):

        # avoid testing the CmdBase class
        patcher = mock.patch(
            'libpipe.cmds.CmdBase',
            new=mock.MagicMock(
                link=mock.Mock(),
                cmd=mock.Mock(return_value='cmd --foo'),
                wrap=None,
            ),
        )
        self.mock_cmd = patcher.start()
        self.addCleanup(patcher.stop)
