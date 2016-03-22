'''Test the CmdDummy class'''

import unittest

from libpipe.cmd.dummy import CmdDummy


class TestCmdDummy(unittest.TestCase):

    def test_init_does_not_need_parameters(self):

        CmdDummy()  # should not raise
