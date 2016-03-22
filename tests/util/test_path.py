'''
Created on Apr 29, 2015

@author: biiremployee
'''

import unittest
from unittest import mock

from libpipe.util import path


import logging
log = logging.getLogger(__name__)


class TestPath(unittest.TestCase):

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

        pass

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

        pass

    def test__walk_file__Extension_suffix_string__2files(self):
        '''Test walk_file retrieve files given a suffix (> extension)'''
        contents = path.walk_file("/analyses",
                                  extension="_hits.bam")
        expected = [
            '/analyses/fpkm/sample01/tophat/accepted_hits.bam',
            '/analyses/fpkm/sample02/tophat/accepted_hits.bam',
        ]

        self.assertListEqual(contents, expected)

        pass

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

        pass

    def test__walk_file__Extension_list__NoFiles(self):
        '''Test walk_file retrieve all files'''
        contents = path.walk_file("/analyses",
                                  extension=['sam', 'xls'])
        expected = []

        self.assertListEqual(contents, expected)

        pass

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
        pass

    def test__walk_file__Pattern_lookbehind__2files(self):
        '''Test walk_file retrieve all files'''
        contents = path.walk_file("/analyses",
                                  pattern='(?<!tophat)\/\w*?\.bam')
        expected = [
            '/analyses/fpkm/sample01/sample01.bam',
            '/analyses/fpkm/sample02/sample02.bam',
        ]

        self.assertListEqual(contents, expected)

        pass

    def test__walk_file__Pattern_lookbehind2__2files(self):
        '''Test walk_file retrieve all files'''
        contents = path.walk_file("/analyses",
                                  pattern='(?<!accepted_hits)\.bam')
        expected = [
            '/analyses/fpkm/sample01/sample01.bam',
            '/analyses/fpkm/sample02/sample02.bam',
        ]

        self.assertListEqual(contents, expected)

        pass

    def test__walk_file__Pattern_lookahead__2files(self):
        '''Test walk_file retrieve all files'''
        contents = path.walk_file("/analyses",
                                  pattern='sample(?=\d+\.bam)')
        expected = [
            '/analyses/fpkm/sample01/sample01.bam',
            '/analyses/fpkm/sample02/sample02.bam',
        ]

        self.assertListEqual(contents, expected)

        pass

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

        pass

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

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
