# import argparse
# import unittest
# import os.path

from unittest import mock

from libpipe import argp
from tests.argp.base import ArgpTestCase

import logging
log = logging.getLogger(__name__)


class TestArgpPipeTestCase(ArgpTestCase):

    def setUp(self):
        super().setUp()
        self.PARSER = argp.pipe.pipe_parser  # MUST be in child __init__

        # prevents addition of pwd while testing due to abspath
        # -- (poor?) design choices were made to avoid this checking...
        patcher = mock.patch(
            'os.path.abspath', side_effect=lambda x: 'abs/' + x)
        patcher.start()
        self.addCleanup(patcher.stop)

        # prevents addition of pwd while testing due to abspath
        # -- (poor?) design choices were made to avoid this checking...
        patcher = mock.patch(
            'libpipe.util.path.protect', side_effect=lambda x: 'protect/' + x)
        self.mock_protect = patcher.start()
        self.addCleanup(patcher.stop)

        self.default = {
            'project': 'my_project',
            'root': '~/foo/bar',
            'summary': 'path/to/summary.txt',
            'data': 'path/to/data',
            'genome': 'genome/index_ref',
            'filter': 'genome/index_filter',
            # 'debug': '',  # flag only -- test elsewhere
        }


class TestArgpPipe(TestArgpPipeTestCase):

    def test_basic_args_with_variables_set(self):
        self.setup_mock_read('summary')
        test_args = ['protect']
        expected = {k: v for k, v in self.default.items() if k in test_args}
        args = self.get_args(self.default)
        for name, val in expected.items():
            with self.subTest(arg=name):
                self.assertIn(name, args)
                self.assertEqual(getattr(args, name), val)

    def test_path_args_with_variables_set(self):
        self.setup_mock_read('summary')
        test_args = ['protect']
        expected = {
            k: self.mock_protect(v)
            for k, v in self.default.items()
            if k in test_args
        }
        args = self.get_args(self.default)
        for name, val in expected.items():
            with self.subTest(arg=name):
                self.assertIn(name, args)
                self.assertEqual(getattr(args, name), val)

    def test_path_list_args_with_variables_set(self):
        self.setup_mock_read('summary')
        test_args = ['filter', 'genome']
        expected = {
            k: [self.mock_protect(v)]
            for k, v in self.default.items()
            if k in test_args
        }
        args = self.get_args(self.default)
        for name, val in expected.items():
            attr = name + '_list'
            with self.subTest(arg=name):
                self.assertIn(name + '_list', args)
                self.assertEqual(getattr(args, attr), val)

    def test_pipe_parser_calls_path_protect_abspath_for_each_arg_dir(self):
        self.setup_mock_read('summary')
        self.default['filter'] = 'genome/index_filter'  # add filter
        dir_args = {
            k: v for k, v in self.default.items()
            if '/' in v
        }

        # with mock.patch('libpipe.util.path.protect') as mock_protect:
        self.get_args(self.default)
        expected = [
            mock.call(v) for v in dir_args.values()]
        self.mock_protect.assert_has_calls(expected, any_order=True)

    def test_pipe_parser_stores_genome_args_in_genome_list(self):
        args = self.get_args('--genome=gen/idx0 --genome gen/idx1')
        filters = [self.mock_protect('gen/idx{}'.format(i)) for i in range(2)]
        self.assertEqual(args.genome_list, filters)

    def test_pipe_parser_stores_multiple_filter_args_in_filter_list(self):
        args = self.get_args('--filter=gen/idx0 --filter gen/idx1')
        filters = [self.mock_protect('gen/idx{}'.format(i)) for i in range(2)]
        self.assertEqual(args.filter_list, filters)


class TestArgpPipe__build_args(TestArgpPipeTestCase):

    '''Test the function argp.pipe.build_args

    `build_args` should be set as a default on the arguments object.

    '''

    def test_pipe_parser_sets_build_args_on_args_obj(self):
        args = self.get_args('')
        self.assertIn('build_args', args)

    def test_build_args_renames_summary_to_summary_file(self):
        self.setup_mock_read('A\tA1.fq;A2.fq\nB\tB1.fq;B2.fq\n')
        args = self.get_args('--summary=data/summary.txt')
        summary = args.summary
        args.build_args(args)

        self.assertIn('summary_file', args)
        self.assertEqual(args.summary_file, summary)

    def test_build_args_opens_args_summary_file(self):
        mock_open = self.setup_mock_read('A\tA1.fq;A2.fq\nB\tB1.fq;B2.fq\n')
        args = self.get_args('--summary=data/summary.txt')
        args.build_args(args)
        mock_open.assert_called_once_with(
            self.mock_protect('data/summary.txt'), 'r')

    def test_build_args_raises_FileNotFoundError_if_summary_file_DNE(self):
        dat = 'A\tA1.fq;A2.fq\nB\tB1.fq\n'
        self.setup_mock_read(dat, side_effect=FileNotFoundError)
        args = self.get_args('--summary=data/summary.txt')
        with self.assertRaises(FileNotFoundError):
            args = argp.pipe.build_args(args)

    def test_build_args_sets_summary_as_dict_of_lists(self):
        dat = 'A\tA1.fq;A2.fq\nB\tB1.fq\n'
        self.setup_mock_read(dat)
        args = self.get_args('--summary=data/summary.txt')
        args = argp.pipe.build_args(args)

        expected = {
            'A': ['A1.fq', 'A2.fq'],
            'B': ['B1.fq', ],
        }
        self.assertEqual(args.summary, expected)
