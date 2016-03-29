
import argparse
import os.path

from unittest import mock

from libpipe import argp
from tests.base import LibpipeTestCase

import pipe as driver

import logging
log = logging.getLogger(__name__)


class TestPipeDriver(LibpipeTestCase):

    def setUp(self):
        super().setUp()
        self.PARSER = argp.pipe.pipe_parser  # MUST be in child __init__

        patcher = mock.patch.object(driver.log, 'info')  # block unwanted msgs
        patcher.start()
        self.addCleanup(patcher.stop)

    def get_args(self, arg_str, **kwargs):
        '''Setup a paser using child-defined function and return args

        Attributes:
            arg_str: A string or dict (long opt only) of args.
                e.g, '-f file.txt --odir foo/bar', or
                {'dir': 'foo/bar', 'project': 'my_project'}
            kwargs: args to pass directly to the parser
        '''

        parser = driver.setup_parser(argparse.ArgumentParser())
        mock.patch.object(parser, 'print_usage')

        try:
            args = parser.parse_args(arg_str.split())
        except AttributeError:
            arg_str = ' '.join(
                '--{} {}'.format(k, v) for k, v in arg_str.items())
            args = parser.parse_args(arg_str.split())
        return args

    def test_main_passes_summary_arg_dict_to_run_pipes(self):
        dat = 'A\tA1.fq;A2.fq\nB\tB1.fq\n'
        self.setup_mock_read(dat)
        args = self.get_args('--summary=data/summary.txt')

        expected = {
            'A': ['A1.fq', 'A2.fq'],
            'B': ['B1.fq', ],
        }

        with mock.patch.object(driver, 'run_pipes') as mock_run:
            driver.main(args)

        mock_run.assert_called_once_with(expected)

    def test_main_passes_dir_dict_to_run_pipes(self):
        dir_files = ['seq{}.fq'.format(i) for i in range(3)]
        m = mock.Mock(return_value=dir_files)
        args = self.get_args('-d foo/bar')
        with mock.patch('libpipe.util.path.walk_file', m):
            args.find_files(args)

        expected = {
            'seq0': ['seq0.fq'],
            'seq1': ['seq1.fq'],
            'seq2': ['seq2.fq'],
        }

        with mock.patch.object(driver, 'run_pipes') as mock_run:
            driver.main(args)
        mock_run.assert_called_once_with(expected)
