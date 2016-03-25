

import argparse
import os
import unittest
from unittest import mock

from libpipe import argp
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

        patcher = mock.patch('argparse.ArgumentParser.print_usage')
        patcher.start()
        self.addCleanup(patcher.stop)

        # prevent unwanted printing of error message (but same effect!)
        patcher = mock.patch(
            'argparse.ArgumentParser.error', side_effect=SystemExit)
        patcher.start()
        self.addCleanup(patcher.stop)
        pass

    def get_args(self, arg_str, **kwargs):
        parser = self.PARSER(**kwargs)
        mock.patch.object(parser, 'print_usage')
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


class TestArgpInput(ArgpTestCase):

    def setUp(self):
        super().setUp()
        self.PARSER = argp.io.input_parser

    def test_input_parser_stores_f_and_file_args_in_file_list(self):
        args = self.get_args('-f file1.txt --file file2.txt')
        self.assertEqual(args.file_list, ['file1.txt', 'file2.txt'])

    def test_input_parser_calls_path_protect_for_each_input_file(self):
        with mock.patch('libpipe.util.path.protect') as mock_protect:
            self.get_args('-f file0.txt --file file1.txt')
        expected = [mock.call('file{}.txt'.format(i)) for i in range(2)]
        mock_protect.assert_has_calls(expected)

    def test_input_parser_stores_d_and_dir_args_in_dir_list(self):
        args = self.get_args('-d dir/subdir0 --dir dir/subdir1')
        dirs = [os.path.join(os.getcwd(), 'dir/subdir{}'.format(i))
                for i in range(2)]
        self.assertEqual(args.dir_list, dirs)

    def test_find_files_adds_dir_files_to_file_list(self):
        dir_files = ['f{}.txt'.format(i) for i in range(3)]
        m = mock.Mock(return_value=dir_files)
        args = self.get_args('-f file.txt -d data/foo')
        with mock.patch('libpipe.util.path.walk_file', m):
            args.find_files(args)
        self.assertEqual(args.file_list, ['file.txt'] + dir_files)

    def test_find_files_deletes_dir_list_attribute(self):
        args = self.get_args('-d data/foo')
        self.assertIn('dir_list', args)  # checks namespace
        with mock.patch('libpipe.util.path.walk_file'):
            args.find_files(args)
        self.assertNotIn('dir_list', args)  # checks namespace

    def test_input_parser_stores_extn_in_extn_list_with_leading_period(self):
        args = self.get_args('--extn txt --extn .fq')
        self.assertEqual(args.extn_list, ['.txt', '.fq'])

    def test_input_parser_stores_pattern_in_pattern_list(self):
        # implicit test: removes quotes
        args = self.get_args('--pattern "*.txt" --pattern \'_1.fq\'')
        self.assertEqual(args.pattern_list, ['*.txt', '_1.fq'])

    def test_find_files_filters_dir_files_by_extension(self):
        extns = ['fq', 'txt', 'fastq']
        dir_files = ['f.{}'.format(extns[i]) for i in range(3)]
        m = mock.Mock(return_value={'file': dir_files})
        args = self.get_args('-f file.txt -d data/foo --extn fq')
        with mock.patch('libpipe.util.path.walk_safe', m):
            args.find_files(args)
        self.assertEqual(args.file_list, ['file.txt', 'f.fq'])

    def test_find_files_filters_dir_files_by_pattern(self):
        extns = ['fq', 'txt', 'fastq']
        dir_files = ['f.{}'.format(extns[i]) for i in range(3)]
        m = mock.Mock(return_value={'file': dir_files})
        args = self.get_args('-f file.txt -d dat/foo --pattern="\.f(ast)?q$"')
        with mock.patch('libpipe.util.path.walk_safe', m):
            args.find_files(args)
        self.assertEqual(args.file_list, ['file.txt', 'f.fq', 'f.fastq'])

    def test_raises_SystemExit_if_accept_dir_False_and_dir_opt_given(self):
        opts = ['-d dat/foo', '--extn .foo', '--pattern="foo"']
        for opt in opts:
            with self.subTest(opt=opt):
                with self.assertRaises(SystemExit):
                    self.get_args(opt, accept_dirs=False)

    def test_raises_SystemExit_if_accept_file_False_and_file_opt_given(self):
        opts = ['-f file.txt']
        for opt in opts:
            with self.subTest(opt=opt):
                with self.assertRaises(SystemExit):
                    self.get_args(opt, accept_files=False)


class TestArgpOutput(ArgpTestCase):

    def setUp(self):
        super().setUp()
        self.PARSER = argp.io.output_parser

    def test_sets_ofile_if_o_or_out_given(self):
        opts = ['-o out.txt', '--out out.txt']
        for opt in opts:
            with self.subTest(opt=opt):
                # implicit check for accept_files=True
                args = self.get_args(opt)
                self.assertEqual(args.ofile, 'out.txt')

    def test_sets_odir_if_odir_given(self):
        args = self.get_args('--odir dat/out', accept_dirs=True)
        self.assertEqual(args.odir, 'dat/out')

    def test_raises_SystemExit_if_same_out_opt_given_gt_1x(self):
        opts = ['-o f.txt --out f2.txt', '--odir a/b --odir b/c']
        for opt in opts:
            with self.subTest(opt=opt):
                with self.assertRaises(SystemExit):
                    self.get_args(opt, accept_dirs=True)

    def test_raises_SystemExit_if_accept_dir_False_and_outdir_opt_given(self):
        arg_str = '--outdir dat/foo'
        with self.assertRaises(SystemExit):
            self.get_args(arg_str, accept_dirs=False)

    def test_raises_SystemExit_if_accept_file_False_n_outfile_opt_given(self):
        opts = ['-o file.txt', '--out file2.txt']
        for opt in opts:
            with self.subTest(opt=opt):
                with self.assertRaises(SystemExit):
                    self.get_args(opt, accept_files=False)
