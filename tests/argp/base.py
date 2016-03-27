

import argparse
import unittest
from unittest import mock


# from libpipe.util import path

import logging
log = logging.getLogger(__name__)


class ArgpTestCase(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        if cls is ArgpTestCase:
            raise unittest.SkipTest("Skip BaseTest tests, it's a base class")
        super(ArgpTestCase, cls).setUpClass()

    def setUp(self):

        # prevent unwanted usage message
        patcher = mock.patch('argparse.ArgumentParser.print_usage')
        patcher.start()
        self.addCleanup(patcher.stop)

        # prevent unwanted printing of error message (but same effect!)
        patcher = mock.patch(
            'argparse.ArgumentParser.error',
            side_effect=SystemExit('parser error')
        )
        patcher.start()
        self.addCleanup(patcher.stop)
        pass

    def get_args(self, arg_str, **kwargs):
        '''Setup a paser using child-defined function and return args

        Attributes:
            arg_str: A string or dict (long opt only) of args.
                e.g, '-f file.txt --odir foo/bar', or
                {'dir': 'foo/bar', 'project': 'my_project'}
            kwargs: args to pass directly to the parser
        '''

        parser = self.PARSER(**kwargs)
        mock.patch.object(parser, 'print_usage')

        try:
            args = parser.parse_args(arg_str.split())
        except AttributeError:
            arg_str = ' '.join(
                '--{} {}'.format(k, v) for k, v in arg_str.items())
            args = parser.parse_args(arg_str.split())
        return args

    #
    #    Common tests
    #

    def test_creates_new_parser_on_every_call(self):
        parser1 = self.PARSER()
        parser2 = self.PARSER()
        self.assertIsInstance(parser1, argparse.ArgumentParser)
        self.assertNotEqual(parser1, parser2)

    def test_does_not_create_new_parser_if_given(self):
        parser = argparse.ArgumentParser()
        with mock.patch('argparse.ArgumentParser') as mock_parser:
            self.PARSER(parser)
        self.assertEqual(mock_parser.call_count, 0)
