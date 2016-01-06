import os.path
import unittest
from unittest.mock import patch

import libpipe
from libpipe.cmds.utility import FastqcCmd

import logging
log = logging.getLogger(__name__)


class TestFastqcCmd(unittest.TestCase):

    def test_output_list_starts_with_input(self):
        self.fail()
