
import argparse
import unittest
from unittest import mock

from libpipe import argp

import logging
log = logging.getLogger(__name__)


class ArgpTestCase(unittest.TestCase):
    pass


class TestArgpCore(ArgpTestCase):

    def test_parser_always_returns_same_argparse_object(self):
        parser1 = argp.core.parser()
        parser2 = argp.core.parser()
        self.assertIsInstance(parser1, argparse.ArgumentParser)
        self.assertEqual(parser1, parser2)
