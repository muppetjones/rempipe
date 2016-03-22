'''Test the PipeBase class'''

import unittest
from unittest import mock


from libpipe.pipe import PipeBase


class PipeBaseTestCase(unittest.TestCase):

    def setUp(self):

        # avoid testing the CmdBase class
        patcher = mock.patch(
            'libpipe.cmd.CmdBase',
            new=mock.MagicMock(
                link=mock.Mock(),
                cmd=mock.Mock(return_value='cmd --foo'),
                wrap=None,
            ),
        )
        self.mock_cmd = patcher.start()
        self.addCleanup(patcher.stop)

    def test_init_creates_a_dummy_cmd_for_linking(self):
        '''Test that init creates a cmd obj to use for linking'''

        with mock.patch('libpipe.cmd.CmdDummy') as mock_dummy:
            pipe = PipeBase()
        mock_dummy.assert_called_once_with()
