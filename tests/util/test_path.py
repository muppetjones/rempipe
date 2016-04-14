'''
Created on Apr 29, 2015

@author: biiremployee
'''

import os
import unittest
from unittest import mock

from libpipe.util import path


import logging
log = logging.getLogger(__name__)


class TestProtect(unittest.TestCase):

    def test_expands_user(self):
        given = '~/foo/bar/'
        clean = path.protect(given)
        self.assertNotIn('~', clean)

    def test_removes_trailing_slash(self):
        given = '~/foo/bar/'
        clean = path.protect(given)
        self.assertFalse(clean.endswith('/'), 'Still has trailing slash.')

    def test_raises_AssertionError_for_unsafe_char(self):
        given = '~/foo/b?ar/'
        with self.assertRaises(AssertionError):
            path.protect(given)

    def test_converts_to_absolute_path_if_dir_detected(self):
        given = 'foo/bar'
        clean = path.protect(given)
        self.assertTrue(clean.startswith(os.getcwd()), 'Not absolutized.')

    def test_does_not_convert_to_absolute_path_if_no_dir_detected(self):
        given = 'foo_bar'
        clean = path.protect(given)
        self.assertFalse(clean.startswith(os.getcwd()), 'Is absolutized.')

    def test_does_not_convert_to_absolute_path_if_abspath_is_False(self):
        given = 'foo/bar'
        clean = path.protect(given, abspath=False)
        self.assertFalse(clean.startswith(os.getcwd()), 'Is absolutized.')


class TestWalkFile(unittest.TestCase):

    def setUp(self):
        # create os.walk return values based on following directory structure
        # /analyses
        #    /fpkm
        #        /sample0N
        #            /tophat
        #                accepted_hits.bam
        #            sample0N.bam
        #            sample0N.count
        #    /deseq
        #        results.csv
        #    /reads
        #        /samples
        #            sample0N.fq
        #    study_design.txt
        mock_tuples = iter([('/analyses',
                             ['deseq', 'fpkm', 'reads', ],
                             ['study_design.txt', ],
                             ),
                            ('/analyses/deseq',
                             [],
                             ['results.csv', ],
                             ),
                            ('/analyses/fpkm',
                             ['sample01', 'sample02', ],
                             [],
                             ),
                            ('/analyses/fpkm/sample01',
                             ['tophat', ],
                             ['sample01.bam', 'sample01.count', ],
                             ),
                            ('/analyses/fpkm/sample01/tophat',
                             [],
                             ['accepted_hits.bam', ],
                             ),
                            ('/analyses/fpkm/sample02',
                             ['tophat', ],
                             ['sample02.bam', 'sample02.count', ],
                             ),
                            ('/analyses/fpkm/sample02/tophat',
                             [],
                             ['accepted_hits.bam', ],
                             ),
                            ('/analyses/reads',
                             [],
                             ['sample01.fq', 'sample02.fq', ],
                             ),
                            ])

        patch_walk = mock.patch('os.walk')
        self.addCleanup(patch_walk.stop)
        self.mock_walk = patch_walk.start()
        self.mock_walk.return_value = mock_tuples

        # ------------------------------ because we're giving a fake directory,
        # ------------------------------------------ we also need to mock isdir
        patch_isdir = mock.patch('os.path.isdir')
        self.addCleanup(patch_isdir.stop)
        self.mock_isdir = patch_isdir.start()
        self.mock_isdir.return_value = True

        pass

    def tearDown(self):
        pass

    #
    #   walk_safe
    #

    def test_walk_safe_passes_AttributeError_from_protect(self):
        with mock.patch.object(path, 'protect', side_effect=AssertionError):
            with self.assertRaises(AssertionError):
                path.walk_safe('foo/bar')

    def test_walk_safe_level_1_returns_contents_from_given_dir_only(self):
        contents = path.walk_safe("/analyses", level=1)
        expected = {
            'file': [
                '/analyses/study_design.txt',
            ],
            'dir': [
                '/analyses/deseq',
                '/analyses/fpkm',
                '/analyses/reads',
            ]
        }
        self.assertEqual(contents, expected)

    #
    #   walk_dir
    #

    def test_walk_dir_returns_all_directories(self):
        contents = path.walk_dir("/analyses")
        expected = [
            '/analyses/deseq',
            '/analyses/fpkm',
            '/analyses/reads',
            '/analyses/fpkm/sample01',
            '/analyses/fpkm/sample02',
            '/analyses/fpkm/sample01/tophat',
            '/analyses/fpkm/sample02/tophat',
        ]
        self.assertEqual(contents, expected)

    def test_walk_dir_returns_dir_matching_one_of_mult_pattern(self):
        contents = path.walk_dir("/analyses", pattern=['re', '02'])
        expected = [
            '/analyses/reads',
            '/analyses/fpkm/sample02',
            '/analyses/fpkm/sample02/tophat',
        ]
        self.assertEqual(contents, expected)

    #
    #   walk_file
    #

    def test__walk_file__AllFiles(self):
        '''Test walk_file retrieve all files'''
        contents = path.walk_file("/analyses")
        expected = [
            '/analyses/study_design.txt',
            '/analyses/deseq/results.csv',
            '/analyses/fpkm/sample01/sample01.bam',
            '/analyses/fpkm/sample01/sample01.count',
            '/analyses/fpkm/sample01/tophat/accepted_hits.bam',
            '/analyses/fpkm/sample02/sample02.bam',
            '/analyses/fpkm/sample02/sample02.count',
            '/analyses/fpkm/sample02/tophat/accepted_hits.bam',
            '/analyses/reads/sample01.fq',
            '/analyses/reads/sample02.fq',
        ]

        self.assertListEqual(contents, expected)

    #=========================================================================
    # Test walk_file with extension
    #=========================================================================

    def test__walk_file__Extension_string__4files(self):
        '''Test walk_file retrieve files given extension string'''
        contents = path.walk_file("/analyses",
                                  extension="bam")
        expected = [
            '/analyses/fpkm/sample01/sample01.bam',
            '/analyses/fpkm/sample01/tophat/accepted_hits.bam',
            '/analyses/fpkm/sample02/sample02.bam',
            '/analyses/fpkm/sample02/tophat/accepted_hits.bam',
        ]

        self.assertListEqual(contents, expected)

    def test__walk_file__Extension_suffix_string__2files(self):
        '''Test walk_file retrieve files given a suffix (> extension)'''
        contents = path.walk_file("/analyses",
                                  extension="_hits.bam")
        expected = [
            '/analyses/fpkm/sample01/tophat/accepted_hits.bam',
            '/analyses/fpkm/sample02/tophat/accepted_hits.bam',
        ]

        self.assertListEqual(contents, expected)

    def test__walk_file__Extension_list__4files(self):
        '''Test walk_file retrieve all files'''
        contents = path.walk_file("/analyses",
                                  extension=['fq', 'count', 'sam'])
        expected = [
            '/analyses/fpkm/sample01/sample01.count',
            '/analyses/fpkm/sample02/sample02.count',
            '/analyses/reads/sample01.fq',
            '/analyses/reads/sample02.fq',
        ]

        self.assertListEqual(contents, expected)

    def test__walk_file__Extension_list__NoFiles(self):
        '''Test walk_file retrieve all files'''
        contents = path.walk_file("/analyses",
                                  extension=['sam', 'xls'])
        expected = []

        self.assertListEqual(contents, expected)

    #=========================================================================
    # Test walk_file with pattern
    #=========================================================================

    def test__walk_file__Pattern_str__2files(self):
        '''Test walk_file retrieve all files'''
        contents = path.walk_file("/analyses",
                                  pattern='tophat')
        expected = [
            '/analyses/fpkm/sample01/tophat/accepted_hits.bam',
            '/analyses/fpkm/sample02/tophat/accepted_hits.bam',
        ]

        self.assertListEqual(contents, expected)

        pass

    def test__walk_file__Pattern_list__7files(self):
        '''Test walk_file retrieve all files'''
        contents = path.walk_file("/analyses",
                                  pattern=[r'01\/', r'02'])
        expected = [
            '/analyses/fpkm/sample01/sample01.bam',
            '/analyses/fpkm/sample01/sample01.count',
            '/analyses/fpkm/sample01/tophat/accepted_hits.bam',
            '/analyses/fpkm/sample02/sample02.bam',
            '/analyses/fpkm/sample02/sample02.count',
            '/analyses/fpkm/sample02/tophat/accepted_hits.bam',
            '/analyses/reads/sample02.fq',
        ]

        self.assertListEqual(contents, expected)

    def test__walk_file__Pattern_lookbehind__2files(self):
        '''Test walk_file retrieve all files'''
        contents = path.walk_file("/analyses",
                                  pattern='(?<!tophat)\/\w*?\.bam')
        expected = [
            '/analyses/fpkm/sample01/sample01.bam',
            '/analyses/fpkm/sample02/sample02.bam',
        ]

        self.assertListEqual(contents, expected)

    def test__walk_file__Pattern_lookbehind2__2files(self):
        '''Test walk_file retrieve all files'''
        contents = path.walk_file("/analyses",
                                  pattern='(?<!accepted_hits)\.bam')
        expected = [
            '/analyses/fpkm/sample01/sample01.bam',
            '/analyses/fpkm/sample02/sample02.bam',
        ]

        self.assertListEqual(contents, expected)

    def test__walk_file__Pattern_lookahead__2files(self):
        '''Test walk_file retrieve all files'''
        contents = path.walk_file("/analyses",
                                  pattern='sample(?=\d+\.bam)')
        expected = [
            '/analyses/fpkm/sample01/sample01.bam',
            '/analyses/fpkm/sample02/sample02.bam',
        ]

        self.assertListEqual(contents, expected)

    def test__walk_file__Pattern_Extension_2files(self):
        '''Test walk_file retrieve all files'''
        contents = path.walk_file("/analyses",
                                  extension='bam',
                                  pattern='^(?!.*tophat)')
        expected = [
            '/analyses/fpkm/sample01/sample01.bam',
            '/analyses/fpkm/sample02/sample02.bam',
        ]

        self.assertListEqual(contents, expected)

    #=========================================================================
    # Test walking with level
    #    Unable to test level with current mock implementation
    #    because os.walk is a generator function.
    #=========================================================================

    # @unittest.skip('Need to re-mock to handle generators')
    # def test__walk_file__level0_returns_expected_list(self):
    #     '''Test walk_file level 0'''
    #     contents = path.walk_file("/analyses", level=0)
    #     expected = ['/analyses/study_design.txt',
    #                 ]
    #
    #     self.assertEqual(contents, expected)

    # @unittest.skip('Need to re-mock to handle generators')
    # def test__walk_file__level1_returns_expected_list(self):
    #     '''Test walk_file level 1'''
    #     contents = path.walk_file("/analyses", level=1)
    #     expected = ['/analyses/study_design.txt',
    #                 '/analyses/deseq/results.csv',
    #                 '/analyses/reads/sample01.fq',
    #                 '/analyses/reads/sample02.fq',
    #                 ]
    #
    #     self.assertEqual(contents, expected)

    # @unittest.skip('Need to re-mock to handle generators')
    # def test__walk_file__level2_returns_expected_list(self):
    #     '''Test walk_file level 2'''
    #     contents = path.walk_file("/analyses", level=2)
    #     expected = ['/analyses/study_design.txt',
    #                 '/analyses/deseq/results.csv',
    #                 '/analyses/fpkm/sample01/sample01.bam',
    #                 '/analyses/fpkm/sample01/sample01.count',
    #                 '/analyses/fpkm/sample02/sample02.bam',
    #                 '/analyses/fpkm/sample02/sample02.count',
    #                 '/analyses/reads/sample01.fq',
    #                 '/analyses/reads/sample02.fq',
    #                 ]
    #
    #     self.assertEqual(contents, expected)

    def test_unique_suffix_strict(self):
        file_list = ['~/data/sample/S1.txt',
                     '~/data/sample/S2.txt',
                     '~/data/sample/S3.txt',
                     '~/data/sample/S4.txt',
                     ]
        suffixes = path.unique_suffix(file_list)
        expected = list('1234')

        self.assertEqual(suffixes, expected)

    def test_unique_suffix_relaxed(self):
        file_list = ['~/data/sample/S1.txt',
                     '~/data/sample/S2.txt',
                     '~/data/sample/S3.txt',
                     '~/data/sample/S4.txt',
                     ]
        suffixes = path.unique_suffix(file_list,
                                      strict=False,
                                      suffix_terminator='/',
                                      )
        expected = ['S' + x for x in list('1234')]

        self.assertEqual(suffixes, expected)


class Test_common_directory(unittest.TestCase):

    def test_returns_expected_directory(self):
        files = ['foo/bar/name{}.txt'.format(i) for i in range(5)]
        common = path.common_directory(files)
        self.assertEqual(common, 'foo/bar')

    def test_returns_parent_directory_if_list_contains_a_dir(self):
        '''Test expected behaviour fails if given a directory'''
        files = ['foo/bar/name{}.txt'.format(i) for i in range(5)]
        common = path.common_directory(files + ['foo/bar'])
        self.assertEqual(common, 'foo')


class Test_makedirs(unittest.TestCase):

    def setUp(self):
        patcher = mock.patch('os.makedirs')
        self.mock_makedirs = patcher.start()
        self.addCleanup(patcher.stop)

        patcher = mock.patch('os.path.isdir', return_value=False)
        self.mock_isdir = patcher.start()
        self.addCleanup(patcher.stop)

        patcher = mock.patch.object(path, 'protect')
        self.mock_protect = patcher.start()
        self.addCleanup(patcher.stop)

    def test_raises_OSError_if_bad_path(self):
        self.mock_protect.side_effect = AssertionError()
        with self.assertRaises(OSError):
            path.makedirs('foo/bar')

    def test_raises_OSError_if_os_makedirs_does(self):
        self.mock_makedirs.side_effect = OSError()
        with self.assertRaises(OSError):
            path.makedirs('foo/bar')

    def test_ignores_os_makedirs_err_if_dir_exists(self):
        self.mock_makedirs.side_effect = OSError()
        self.mock_isdir.return_value = True
        path.makedirs('foo/bar')  # should not raise


class Test_avoid_overwrite(unittest.TestCase):

    def setUp(self):
        self.test_str = 'foo/bar.txt'

        def update_x(y):
            self.mock_ctime.x += 1
            return self.mock_ctime.x

        # patch sorting
        patcher = mock.patch('os.path.getctime')
        self.mock_ctime = patcher.start()
        self.mock_ctime.x = 0
        self.mock_ctime.side_effect = update_x
        self.addCleanup(patcher.stop)

        # make abspath return the given value
        patcher = mock.patch('os.path.abspath', side_effect=lambda x: x)
        self.mock_abspath = patcher.start()
        self.addCleanup(patcher.stop)

        # make all files look big
        patcher = mock.patch('os.path.getsize', return_value=10)
        self.mock_getsize = patcher.start()
        self.addCleanup(patcher.stop)

    def mock_path_exists(self, retval=True, side_effect=None):
        if side_effect:
            patcher = mock.patch('os.path.exists', side_effect=side_effect)
        else:
            patcher = mock.patch('os.path.exists', return_value=retval)
        self.addCleanup(patcher.stop)
        return patcher.start()

    # def mock_path_(self, retval=True):
    #     patcher = mock.patch('os.path.exists', retval=retval)
    #     self.addCleanup(patcher.stop)
    #     return patcher.start()

    def test_returns_protected_path_if_unique(self):
        self.mock_path_exists(retval=False)
        with mock.patch.object(path, 'protect', side_effect=lambda x: x) as m:
            unique = path.avoid_overwrite(self.test_str)
        m.assert_called_once_with(self.test_str)
        self.assertEqual(unique, self.test_str)

    def test_ignores_OSError_from_getsize(self):
        n_exist = 6
        return_list = [True] * (n_exist + 1) + [False]
        self.mock_getsize.side_effect = OSError()
        self.mock_path_exists(side_effect=return_list)
        with mock.patch('os.path.isfile', return_value=True) as m:
            unique_name = path.avoid_overwrite(self.test_str)
        self.assertEqual(unique_name, self.test_str + '.{}'.format(n_exist))

    def test_calls_exists_max_number_of_times_if_exists_True(self):
        mock_exists = self.mock_path_exists()
        n_max = 40
        path.avoid_overwrite(self.test_str, max_attempts=n_max)
        # +1 for protected call
        self.assertEqual(mock_exists.call_count, n_max + 1)

    def test_calls_isfile_x_times_if_allow_empty(self):
        self.mock_path_exists()
        n_max = 5
        with mock.patch('os.path.isfile', return_value=True) as m:
            path.avoid_overwrite(
                self.test_str, max_attempts=n_max, allow_empty=True)
        self.assertEqual(m.call_count, n_max)

    def test_does_not_call_isfile_if_not_allow_empty(self):
        self.mock_path_exists()
        n_max = 5
        with mock.patch('os.path.isfile', return_value=True) as m:
            path.avoid_overwrite(
                self.test_str, max_attempts=n_max, allow_empty=False)
        self.assertEqual(m.call_count, 0)

    def test_returns_nonexistant_file_with_number_added(self):
        n_exist = 6
        return_list = [True] * (n_exist + 1) + [False]
        m = self.mock_path_exists(side_effect=return_list)
        unique_name = path.avoid_overwrite(self.test_str)
        self.assertEqual(unique_name, self.test_str + '.{}'.format(n_exist))

    def test_returns_empty_file(self):
        n_exist = 6
        self.mock_getsize.side_effect = [10] * n_exist + [0]
        self.mock_path_exists()
        with mock.patch('os.path.isfile', return_value=True) as m:
            unique_name = path.avoid_overwrite(self.test_str)
        self.assertEqual(unique_name, self.test_str + '.{}'.format(n_exist))

    def test_returns_oldest_file(self):
        '''If max tries exists, return the oldest'''
        m = self.mock_path_exists()
        n_max = 5
        unique_name = path.avoid_overwrite(
            self.test_str, max_attempts=n_max, allow_empty=False)
        self.assertEqual(m.call_count, n_max + 1)
        self.assertEqual(unique_name, self.test_str + '.{}'.format(0))

if __name__ == "__main__":
    # import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
