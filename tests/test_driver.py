
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

        # prevents addition of pwd while testing due to abspath
        # -- (poor?) design choices were made to avoid this checking...
        patcher = mock.patch('os.path.abspath', side_effect=lambda x: x)
        patcher.start()
        self.addCleanup(patcher.stop)

    def mock_protect(self):
        patcher = mock.patch(
            'libpipe.util.path.protect', side_effect=lambda x: x)
        self.addCleanup(patcher.stop)
        return patcher.start()

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
        # Do not create an argp object and pass it. Doing so
        # tests the ability of argp to create a dict from a
        # summary file. Instead, check that args.summary is
        # passed to run_pipes.
        args = mock.MagicMock(
            summary={k: v for v, k in enumerate(list('abc'))},
            data=None,
            genome_list=None,
            project='project',
            root='root',
        )

        with mock.patch.object(driver, 'run_pipes') as mock_run:
            driver.main(args)

        mock_run.assert_called_once_with(
            args.summary, None, data=None, odir='root/project/samples')

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
        mock_run.assert_called_once_with(expected, None, data=None, odir=None)

    def test_main_passes_out_dir_to_run_pipes(self):
        self.mock_protect()
        args = self.get_args('--project foo --root ~/projects')
        setattr(args, 'summary', {})
        expected = os.path.join('~/projects', 'foo', 'samples')
        with mock.patch.object(driver, 'run_pipes') as mock_run:
            driver.main(args)
        mock_run.assert_called_once_with({}, None, data=None, odir=expected)

    def test_main_calls_AlignPipe_with_genome_in_input(self):
        args = self.get_args('--genome genome/hisat_index')
        setattr(args, 'file_list', ['a.fq'])
        with mock.patch('libpipe.pipe.align.AlignPipe') as mock_pipe:
            driver.main(args)

        expected = args.genome_list + args.file_list
        mock_pipe.assert_called_once_with(input=expected)

    def test_data_dir_is_added_to_file_names_if_given(self):
        args = self.get_args('--data=foo/bar --genome genome/hisat_index')
        setattr(args, 'file_list', ['a.fq'])
        with mock.patch('libpipe.pipe.align.AlignPipe') as mock_pipe:
            driver.main(args)

        expected = args.genome_list + [
            os.path.join('foo/bar', f) for f in args.file_list]
        mock_pipe.assert_called_once_with(input=expected)
