
import unittest


from libpipe.cmd import CmdInterface


class TestCmdInterface(unittest.TestCase):

    def test_cannot_be_instatiated_directly(self):
        with self.assertRaises(TypeError):
            CmdInterface()
