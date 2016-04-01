
# import unittest
from unittest import mock


# from libpipe.util import path
from libpipe.type import index
from tests import base
from tests.type import test_base

import logging
log = logging.getLogger(__name__)


class IndexTypeTestCase(base.LibpipeTestCase):

    def setUp(self):
        # create our test index
        extns = {'.foo': 4, '.bar': 2}
        self.INDEX = index.IndexType.factory(extns)

        # setup example output from walk_file (or walk_safe)
        n_foo, n_bar = (4, 2)
        self.default_file_list = (
            ['{}.foo'.format(i) for i in range(n_foo)] +
            ['{}.bar'.format(i) for i in range(n_bar)]
        )

    def setup_mock_walk_safe(self, retval=None):
        patcher = mock.patch('libpipe.util.path.walk_safe')
        m = patcher.start()
        if not retval:
            retval = {'file': list(self.default_file_list)}
        m.return_value = retval
        self.addCleanup(patcher.stop)
        return m

    def setup_mock_walk_file(self, retval=None):
        patcher = mock.patch('libpipe.util.path.walk_file')
        m = patcher.start()
        if not retval:
            retval = list(self.default_file_list)
        m.return_value = retval
        self.addCleanup(patcher.stop)
        return m

    def setup_mock_check_prefix(self):
        patcher = mock.patch.object(index.IndexType, '_check_prefix')
        m = patcher.start()
        self.addCleanup(patcher.stop)
        return m

    def setup_mock_check_extns(self):
        patcher = mock.patch.object(index.IndexType, '_check_extns')
        m = patcher.start()
        self.addCleanup(patcher.stop)
        return m


#
# test base and factory
#

class TestIndexType__TypeBase(test_base.TestTypeBase):
    pass

    def setUp(self):
        super().setUp()
        self.TYPE = index.IndexType

        patcher = mock.patch.object(index.IndexType, '_check_prefix')
        patcher.start()
        self.addCleanup(patcher.stop)

        patcher = mock.patch.object(index.IndexType, '_check_extns')
        patcher.start()
        self.addCleanup(patcher.stop)


class TestIndexType(IndexTypeTestCase):

    '''Test that IndexType fulfills TypeBase tests'''

    def setUp(self):
        self.TYPE = index.IndexType
        extns = {'.foo': 4, '.bar': 2}
        self.INDEX = self.TYPE.factory(extns)

        self.setup_mock_check_extns()
        self.setup_mock_check_prefix()

    def test_str_index_returns_given_file(self):
        _index = self.INDEX('foo/bar')
        self.assertEqual(str(_index), 'foo/bar')


class TestIndexType_factory(base.LibpipeTestCase):

    def test_factory_sets_expected_extns_attr_to_empty_list(self):
        index_obj = index.IndexType.factory()
        self.assertEqual(index_obj.extns, [])

    def test_factory_sets_expected_counts_attr_to_empty_list(self):
        index_obj = index.IndexType.factory()
        self.assertEqual(index_obj.counts, [])

    def test_factory_list_sets_extns_in_alpha_order(self):
        expected = ['.rev.bt2', '.bt2', ]
        _index = index.IndexType.factory(expected)
        self.assertEqual(_index.extns, sorted(expected))

    def test_factory_multi_sets_extns_in_alpha_order(self):
        expected = ['.rev.bt2', '.bt2', ]
        _index = index.IndexType.factory(*expected)
        self.assertEqual(_index.extns, sorted(expected))

    def test_factory_sets_expected_attr_from_kwargs(self):
        extn_dict = {'.bt2': 4, '.rev.bt2': 2, '.foo': None, '.bar': None}

        # repeat check multiple times to ensure consistent
        for i in range(50):
            _index = index.IndexType.factory(**extn_dict)

            # extns
            expected = sorted(list(extn_dict.keys()))
            self.assertEqual(_index.extns, expected)

            # counts
            expected = [extn_dict[k] for k in sorted(list(extn_dict.keys()))]
            self.assertEqual(_index.counts, expected)

    def test_factory_dict_sets_extns_in_alpha_order(self):
        extn_dict = {'.bt2': 4, '.rev.bt2': 2, '.foo': None, '.bar': None}

        # repeat check multiple times to ensure consistent
        for i in range(50):
            _index = index.IndexType.factory(extn_dict)
            expected = sorted(list(extn_dict.keys()))
            self.assertEqual(_index.extns, expected)

    def test_factory_dict_sets_counts_in_extn_order(self):
        extn_dict = {'.bt2': 4, '.rev.bt2': 2, '.foo': None, '.bar': None}

        # repeat check multiple times to ensure consistent
        for i in range(50):
            _index = index.IndexType.factory(extn_dict)
            expected = [extn_dict[k] for k in sorted(list(extn_dict.keys()))]
            self.assertEqual(_index.counts, expected)

    def test_factory_sets_counts_to_1_if_not_given(self):
        extns = ['.rev.bt2', '.bt2', ]
        _index = index.IndexType.factory(*extns)
        expected = [1 for e in extns]
        self.assertEqual(_index.counts, expected)


#
#   Test index checking
#


class TestIndexType__check_extns(IndexTypeTestCase):

    def test_init_calls_walk_file_with_extensions(self):
        mock_walk = self.setup_mock_walk_file()
        _index = self.INDEX('foo/bar')
        expected = mock.call('foo/bar', extension=_index.extns, pattern=[])
        mock_walk.assert_has_calls([expected])

    def test_init_raises_ValueError_if_n_files_w_extn_ne_counts(self):
        mock_walk = self.setup_mock_walk_file()

        tests = [(4, 1), (4, 0), (0, 2), (6, 2)]
        for n_foo, n_bar in tests:
            mock_walk.return_value = (
                ['{}.foo'.format(i) for i in range(n_foo)] +
                ['{}.bar'.format(i) for i in range(n_bar)]
            )
            with self.assertRaises(ValueError):
                self.INDEX('test/dir')

    def test_name_set_as_pattern_for_walk_file(self):
        mock_walk = self.setup_mock_walk_file()
        with mock.patch('os.path.isdir', side_effect=[False, True]):
            _index = self.INDEX('path/to/index')
            mock_walk.assert_called_once_with(
                _index._path, extension=_index.extns, pattern=[_index._name])


class TestIndexType__check_prefix(IndexTypeTestCase):

    def test_prefix_and_name_are_None_if_dir_given(self):
        self.setup_mock_check_extns()
        with mock.patch('os.path.isdir', return_value=True):
            _index = self.INDEX('foo/bar')
        self.assertIsNone(_index._prefix)
        self.assertIsNone(_index._name)

    def test_init_sets_prefix_and_name_if_file_prefix_given(self):
        '''Test for detection of an incomplete file name'''

        self.setup_mock_check_extns()
        with mock.patch('os.path.isdir', side_effect=[False, True]):
            _index = self.INDEX('foo/bar')
        self.assertEqual(_index._path, 'foo')
        self.assertEqual(_index._prefix, 'foo/bar')
        self.assertEqual(_index._name, 'bar')

    def test_str_index_returns_given_prefix(self):
        self.setup_mock_check_extns()
        with mock.patch('os.path.isdir', side_effect=[False, True]):
            _index = self.INDEX('foo/bar')
        self.assertEqual(str(_index), 'foo/bar')

#     def test_sandbox(self):
#         _index = self.INDEX('foo/bar')
