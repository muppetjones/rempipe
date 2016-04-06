
import os.path
from unittest import mock


# from libpipe.util import path
from libpipe.type import base as _type
from libpipe.type import index
from tests import base
from tests.type import test_base

import logging
log = logging.getLogger(__name__)


class IndexTypeTestCase(base.LibpipeTestCase):

    def setUp(self):
        super().setUp()
        self.PARENT = index.IndexType
        self.FACTORY = index.factory
        self.CHILD = self.FACTORY('IndexSubType')

#
# test base and factory
#


class TestIndexMeta(base.LibpipeTestCase):

    def test_IndexMeta_inherits_from_TypeMeta(self):
        self.assertTrue(
            issubclass(index.IndexMeta, _type.TypeMeta),
            "IndexMeta is NOT subclassed from TypeMeta",
        )


class TestIndexType__TypeBase(test_base.TestTypeBase):
    pass

    def setUp(self):
        super().setUp()
        self.PARENT = index.IndexType
        self.FACTORY = index.factory
        self.CHILD = self.FACTORY('IndexSubType')

        patcher = mock.patch('libpipe.util.path.walk_file')
        patcher.start()
        self.addCleanup(patcher.stop)

        # patcher = mock.patch.object(index.IndexType, '_check_extns')
        # patcher.start()
        # self.addCleanup(patcher.stop)

    def test_factory_sets_child_name(self):
        child_class = self.FACTORY(name='Dolly')
        self.assertEqual(child_class.__name__, 'Dolly')

    def test_child_is_still_TypeMeta(self):
        self.assertIsInstance(self.CHILD, _type.TypeMeta)


class TestIndexFactory(IndexTypeTestCase):

    def test_factory_stores_extensions_on_child_class(self):
        extn = ['.foo', '.bar']
        child = self.FACTORY(extns=extn)
        self.assertEqual(child.extns, extn)

    def test_factory_stores_counts_on_child_class(self):
        extn = ['.foo', '.bar']
        counts = [3, 5]
        child = self.FACTORY(extns=extn, counts=counts)
        self.assertEqual(child.counts, counts)

    def test_factory_assumes_all_counts_eq_1_iff_not_given(self):
        extn = ['.foo', '.bar']
        counts = [1] * len(extn)
        child = self.FACTORY(extns=extn)
        self.assertEqual(child.counts, counts)

    def test_factory_raises_ValueError_if_n_counts_ne_n_extn(self):
        extn = ['.foo', '.bar']
        counts = [1, ]
        with self.assertRaises(ValueError):
            self.FACTORY(extns=extn, counts=counts)


class TestIndexSubType(IndexTypeTestCase):

    def setUp(self):
        super().setUp()
        self.child = self.FACTORY(
            name='TestIndex', extns=['.bar', '.foo'], counts=[2, 4])

        n_foo, n_bar = (4, 2)

        self.default_file_list = (
            ['{}.foo'.format(i) for i in range(n_foo)] +
            ['{}.bar'.format(i) for i in range(n_bar)]
        )

    #
    #   Mock setup
    #

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

    def setup_mock_check_extns(self, obj=None):
        obj = obj if obj else self.CHILD
        patcher = mock.patch.object(obj, '_check_extns')
        m = patcher.start()
        self.addCleanup(patcher.stop)
        return m

    #
    #   Tests
    #

    def test_is_subclassed_from_str(self):
        self.assertTrue(issubclass(self.child, str), 'Not subclassed from str')

    def test_init_with_no_args_returns_empty_string(self):
        self.assertEqual(self.child(), '')

    def test_eq_returns_True_for_given_str(self):
        self.setup_mock_check_extns(self.child)

        _index = self.child('foo/bar')
        self.assertEqual(_index, 'foo/bar')

    def test_init_sets_path_from_dirname_of_given(self):
        self.setup_mock_check_extns(self.child)
        path = 'foo/bar'
        index = self.child(path)
        self.assertEqual(index.path, os.path.dirname(path))

    def test_init_sets_name_from_basename_of_given(self):
        self.setup_mock_check_extns(self.child)
        path = 'foo/bar'
        index = self.child(path)
        self.assertEqual(index.name, os.path.basename(path))

    def test_init_calls_walk_file_with_extensions(self):
        mock_walk = self.setup_mock_walk_file()

        _index = self.child('foo/bar')

        expected = mock.call(
            'foo', extension=_index.extns, pattern=['bar'])
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
                self.child('test/dir')
