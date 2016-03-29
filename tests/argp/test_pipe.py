# import argparse
# import unittest
import os.path

from unittest import mock

from libpipe import argp
from tests.argp.base import ArgpTestCase

import logging
log = logging.getLogger(__name__)


class TestArgpPipe(ArgpTestCase):

    def setUp(self):
        super().setUp()
        self.PARSER = argp.pipe.pipe_parser  # MUST be in child __init__

        self.default = {
            'project': 'my_project',
            'root': '~/foo/bar',
            # 'summary': 'path/to/summary.txt',
            'data': 'path/to/data',
            # 'genome': 'genome/index_ref',
            # 'filter': 'genome/index_filter',
            # 'debug': '',  # flag only -- test elsewhere
        }

    def test_basic_args_with_variables_set(self):
        self.setup_mock_read('summary')
        with mock.patch('libpipe.util.path.protect') as mock_protect:
            mock_protect.side_effect = lambda x: x
            args = self.get_args(self.default)
        for name, val in self.default.items():
            with self.subTest(arg=name):
                self.assertIn(name, args)
                self.assertEqual(getattr(args, name), val)

    def test_pipe_parser_calls_path_protect_abspath_for_each_arg_dir(self):
        self.setup_mock_read('summary')
        self.default['filter'] = 'genome/index_filter'  # add filter
        dir_args = {
            k: v for k, v in self.default.items()
            if '/' in v
        }

        with mock.patch('libpipe.util.path.protect') as mock_protect:
            self.get_args(self.default)
        expected = [
            mock.call(v) for v in dir_args.values()]
        mock_protect.assert_has_calls(expected, any_order=True)

    def test_pipe_parser_stores_genome_args_in_genome_list(self):
        args = self.get_args('--genome=gen/idx0 --genome gen/idx1')
        filters = [os.path.join(os.getcwd(), 'gen/idx{}'.format(i))
                   for i in range(2)]
        self.assertEqual(args.genome_list, filters)

    def test_pipe_parser_stores_filter_args_in_filter_list(self):
        args = self.get_args('--filter=gen/idx0 --filter gen/idx1')
        filters = [os.path.join(os.getcwd(), 'gen/idx{}'.format(i))
                   for i in range(2)]
        self.assertEqual(args.filter_list, filters)

    def test_pipe_parser_opens_args_summary_file(self):
        mock_open = self.setup_mock_read('A\tA1.fq;A2.fq\nB\tB1.fq;B2.fq\n')
        self.get_args('--summary=data/summary.txt')
        mock_open.assert_called_once_with(
            os.path.abspath('data/summary.txt'), 'r')

    def test_pipe_parser_sets_summary_as_dict_of_lists(self):
        dat = 'A\tA1.fq;A2.fq\nB\tB1.fq\n'
        self.setup_mock_read(dat)
        args = self.get_args('--summary=data/summary.txt')

        expected = {
            'A': ['A1.fq', 'A2.fq'],
            'B': ['B1.fq', ],
        }
        self.assertEqual(args.summary, expected)

    def test_pipe_parser_raises_FileNotFoundError_if_summary_file_DNE(self):
        dat = 'A\tA1.fq;A2.fq\nB\tB1.fq\n'
        self.setup_mock_read(dat, side_effect=FileNotFoundError)
        with self.assertRaises(FileNotFoundError):
            args = self.get_args('--summary=data/summary.txt')
