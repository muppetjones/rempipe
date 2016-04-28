
import argparse
import os.path

from unittest import mock

from libpipe import argp
from libpipe.pipe import align
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

    def mock_abspath(self):
        # prevents addition of pwd while testing due to abspath
        patcher = mock.patch('os.path.abspath', side_effect=lambda x: x)
        patcher.start()
        self.addCleanup(patcher.stop)

    def mock_expanduser(self):
        patcher = mock.patch('os.path.expanduser', side_effect=lambda x: x)
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

    #
    #   Test summary arg
    #

    def test_main_passes_summary_arg_dict_to_run_pipes(self):
        # Do not create an argp object and pass it. Doing so
        # tests the ability of argp to create a dict from a
        # summary file. Instead, check that args.summary is
        # passed to run_pipes.

        # Impicitly tests addition of summary file dir to each
        # file in the summary
        args = mock.MagicMock(
            summary_file='ha/ha/ha/summary.txt',
            summary={k: str(v) for v, k in enumerate(list('abc'))},
            data=None,
            genome_list=None,
            project='project',
            root='root',
        )

        with mock.patch.object(driver, 'run_pipes') as mock_run:
            driver.main(args)

        mock_run.assert_called_once_with(
            args.summary, None, data=None, odir='root/project/samples')

    #
    #   Test file and dir args
    #

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
        mock_run.assert_called_once_with(expected, [], data=None, odir=None)

    def test_data_dir_is_added_to_file_names_if_given(self):
        self.mock_abspath()
        args = self.get_args('--data=foo/bar --genome genome/hisat_index')
        setattr(args, 'file_list', ['a.fq'])
        with mock.patch('libpipe.pipe.align.AlignPipe') as mock_pipe:
            driver.main(args)

        expected = args.genome_list + [
            os.path.join('foo/bar', f) for f in args.file_list]
        mock_pipe.assert_called_once_with(
            input=expected, job_name='a', odir=None)

    #
    #   Test output directory setup
    #

    def test_main_sets_odir_via_root_project(self):
        self.mock_protect()
        args = self.get_args('--project foo --root ~/projects')
        setattr(args, 'summary', {})
        expected = os.path.join('~/projects', 'foo', 'samples')
        with mock.patch.object(driver, 'run_pipes') as mock_run:
            driver.main(args)
        mock_run.assert_called_once_with({}, [], data=None, odir=expected)

    def test_main_does_not_set_odir_without_root_and_project(self):
        tests = [
            '--root ~/projects',
            '--project foo',
            '',
        ]

        for test in tests:
            with self.subTest(arg_str=test):
                args = self.get_args(test)
                setattr(args, 'summary', {})
                setattr(args, 'genome_list', [])
                with mock.patch.object(driver, 'run_pipes') as mock_run:
                    driver.main(args)
                mock_run.assert_called_once_with(
                    {}, [], data=None, odir=None)

    def test_run_pipe_tries_to_create_odir(self):
        self.mock_expanduser()
        self.mock_abspath()
        args = self.get_args('--project foo --root ~/projects')
        setattr(args, 'summary', {})
        setattr(args, 'genome_list', [])
        with mock.patch('libpipe.pipe.align.AlignPipe') as mock_pipe, \
                mock.patch('libpipe.util.path.makedirs') as m:
            driver.run_pipes(
                {'name': ['file']}, ['genome'], odir='path/to/odir')
        m.assert_called_once_with('path/to/odir/name')

    def test_run_pipes_passes_out_dir_and_job_name_to_pipe_obj(self):
        with mock.patch('libpipe.pipe.align.AlignPipe') as m:
            driver.run_pipes(
                {'name_key': ['file']}, ['genome'], odir='path/to/odir')
        m.assert_called_once_with(
            input=['genome', 'file'],
            job_name='name_key',
            odir='path/to/odir/name_key',
        )

    #
    #   Test genome
    #

    def test_main_calls_AlignPipe_with_genome_in_input(self):
        args = self.get_args('--genome genome/hisat_index')
        setattr(args, 'file_list', ['a.fq'])
        with mock.patch('libpipe.pipe.align.AlignPipe') as mock_pipe:
            driver.main(args)

        expected = args.genome_list + args.file_list
        mock_pipe.assert_called_once_with(
            input=expected, job_name='a', odir=None)


# __END__
