
import os.path

from functools import partial
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

        n_foo, n_bar = (4, 2)
        self.CHILD = self.FACTORY(
            name='TestIndex', extns=['.bar', '.foo'], counts=[n_bar, n_foo])

    #
    #   Mock setup
    #

    def create_file_list(self, extn='.foo', n=4):
        return ['i{}{}'.format(i, extn) for i in range(n)]

    def setup_mock_check_extns(self, obj=None):
        obj = obj if obj else self.CHILD
        patcher = mock.patch.object(obj, '_check_extns')
        m = patcher.start()
        self.addCleanup(patcher.stop)
        return m

    def setup_mock_walk_safe(self, files=[], dirs=[]):
        '''Mock util.path.walk_safe

        CRITICAL: Walk safe should a dict of ALL dir and files.
        '''
        patcher = mock.patch('libpipe.util.path.walk_safe')
        m = patcher.start()
        m.return_value = {'files': files, 'dirs': dirs}
        self.addCleanup(patcher.stop)
        return m

    def setup_mock_walk_file(self, retval=[]):
        '''Mock util.path.walk_file

        CRITICAL: Walk file should return [] if the files were NOT found.
            Otherwise, it returns files that MATCH the desired pattern.
        '''
        patcher = mock.patch('libpipe.util.path.walk_file')
        m = patcher.start()
        self.addCleanup(patcher.stop)
        m.return_value = retval
        return m

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

    def setUp(self):
        super().setUp()
        self.PARENT = index.IndexType
        self.FACTORY = partial(index.factory, extns=['.foo'])
        self.CHILD = self.FACTORY('IndexFooType')  # includes extns=['.foo']

        patcher = mock.patch.object(index.IndexType, '_check_extns')
        patcher.start()
        self.addCleanup(patcher.stop)

    def test_factory_sets_child_name(self):
        child_class = self.FACTORY(name='Dolly')
        self.assertEqual(child_class.__name__, 'Dolly')

    def test_child_is_still_TypeMeta(self):
        self.assertIsInstance(self.CHILD, _type.TypeMeta)


class TestIndexFactory(IndexTypeTestCase):

    def test_factory_return_is_subclassed_from_str_by_default(self):
        subtype = self.FACTORY(extns=list('ab'))
        self.assertTrue(issubclass(subtype, str), 'Not subclassed from str')

    def test_factory_object_is_str_instance_by_default(self):
        subtype = self.FACTORY(extns=list('ab'))
        self.assertIsInstance(subtype(), str)

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

    def test_factory_returns_class_that_is_subclass_of_parent(self):
        subtype = index.IndexMeta(
            'SubIndex', (index.IndexType, ), dict())
        subsub = self.FACTORY(extns=list('ab'), parent=subtype)
        self.assertTrue(
            issubclass(subsub, index.IndexType), 'Not an IndexType subclass')
        self.assertTrue(
            issubclass(subsub, subtype), 'Not a child subclass')

    def test_factory_raises_ValueError_if_n_counts_ne_n_extn(self):
        extn = ['.foo', '.bar']
        counts = [1, ]
        with self.assertRaises(ValueError):
            self.FACTORY(extns=extn, counts=counts)

    def test_factory_raises_ValueError_if_no_extns_given(self):
        with self.assertRaises(ValueError):
            self.FACTORY(extns=[])

    def test_facctory_raises_ValueError_if_parent_is_not_an_IndexType(self):
        with self.assertRaises(ValueError):
            self.FACTORY(extns=list('abc'), parent=int)


class TestIndexSubType(IndexTypeTestCase):

    #
    #   Tests
    #

    def test_init_with_no_args_returns_empty_string(self):
        self.assertEqual(self.CHILD(), '')

    def test_init_sets_path_from_dirname_of_given(self):
        self.setup_mock_check_extns(self.CHILD)
        path = 'foo/bar'
        index = self.CHILD(path)
        self.assertEqual(index.path, os.path.dirname(path))

    def test_init_sets_name_from_basename_of_given(self):
        self.setup_mock_check_extns(self.CHILD)
        path = 'foo/bar'
        index = self.CHILD(path)
        self.assertEqual(index.name, os.path.basename(path))

    def test_init_calls_walk_file_with_extensions(self):
        mock_walk = self.setup_mock_walk_file(
            self.create_file_list('.foo', 4) +
            self.create_file_list('.bar', 2)
        )

        _idx = self.FACTORY('BadIndex', extns=['.foo', '.bar'], counts=[4, 2])
        _idx = self.CHILD('hello/world')

        expected = mock.call(
            'hello', extension=_idx.extns, pattern=['world'])
        mock_walk.assert_has_calls([expected])

    def test_init_raises_ValueError_if_expected_extn_and_cnt_not_found(self):
        _idx = self.FACTORY('BadIndex', extns=['.foo', '.bar'], counts=[4, 2])
        mock_walk = self.setup_mock_walk_file()

        tests = [(4, 1), (4, 0), (0, 2), (6, 2)]
        for n_foo, n_bar in tests:
            mock_walk.return_value = (
                self.create_file_list('.foo', n_foo) +
                self.create_file_list('.bar', n_bar)
            )
            with self.assertRaises(ValueError):
                _idx('test/dir')


class TestIndexSubType__eq(IndexTypeTestCase):

    #
    #   __eq__
    #

    def test_eq_True_for_matching_str(self):
        self.setup_mock_check_extns(self.CHILD)

        _index = self.CHILD('foo/bar')
        self.assertEqual(_index, 'foo/bar')

    def test_eq_False_for_mismatched_str(self):
        self.setup_mock_check_extns(self.CHILD)

        _index = self.CHILD('foo/bar')
        self.assertNotEqual(_index, 'foo/bar2')

    def test_eq_calls_via_IndexType(self):
        self.setup_mock_check_extns(index.IndexType)

        # Does not pass if IndexType or IndexMeta patched...
        with mock.patch.object(self.CHILD, '__eq__') as m:
            _index = self.CHILD('foo/bar')
            _index == 'hi'
        m.assert_called_once_with('hi')

    def test_eq_False_for_diff_extn_same_path_any_type(self):
        self.setup_mock_check_extns(index.IndexType)

        _cls1 = self.FACTORY('IndexSubType1', extns=list('abc'))
        _cls2 = self.FACTORY('IndexSubType2', extns=list('abd'))
        _idx1 = _cls1('foo/bar')
        _idx2 = _cls2('foo/bar')
        self.assertNotEqual(_idx1, _idx2)

    def test_eq_False_for_same_extn_diff_path_any_type(self):
        self.setup_mock_check_extns(index.IndexType)

        _cls1 = self.FACTORY('IndexSubType1', extns=list('abc'))
        _cls2 = self.FACTORY('IndexSubType2', extns=list('abc'))
        _idx1 = _cls1('foo/bar')
        _idx2 = _cls2('foo/bar2')
        self.assertNotEqual(_idx1, _idx2)

    def test_eq_True_for_same_extn_same_path_any_type(self):
        self.setup_mock_check_extns(index.IndexType)

        _cls1 = self.FACTORY('IndexSubType1', extns=list('abc'))
        _cls2 = self.FACTORY('IndexSubType2', extns=list('abc'))
        _idx1 = _cls1('foo/bar')
        _idx2 = _cls2('foo/bar')
        self.assertEqual(_idx1, _idx2)


class TestIndexSubType__instancecheck(IndexTypeTestCase):

    def setUp(self):
        super().setUp()

        self.default_file_list = (
            self.create_file_list('.b', 4) + self.create_file_list('.f', 2)
        )
        self.CHILD = self.FACTORY(
            name='TypeCheckIndex', extns=['.b', '.f'], counts=[4, 2])

    #
    #   __instancecheck__
    #

    def test_returns_True_for_str_if_index_files_found(self):

        m = self.setup_mock_walk_file(self.default_file_list)
        self.assertTrue(
            isinstance('some/path/to/index', self.CHILD),
            'isinstance returns False instead of True'
        )

    def test_returns_False_for_str_if_index_files_not_found(self):
        self.setup_mock_walk_file([])
        self.assertFalse(
            isinstance('some/path/to/index', self.CHILD),
            'isinstance returns True instead of False'
        )

    def test_returns_True_for_factory_objects(self):
        _idx = self.FACTORY('SubIndex', extns=list('abc'))
        self.assertIsInstance(_idx(), index.IndexType)

    def test_returns_False_for_sister_objects(self):
        _idx1 = self.FACTORY('SubIndex', extns=list('abc'))
        _idx2 = self.FACTORY('SubIndex2', extns=list('abc'))
        self.assertNotIsInstance(_idx1(), _idx2)
        self.assertNotIsInstance(_idx2(), _idx1)

    def test_returns_True_for_sub_sub_objects(self):
        _idx1 = self.FACTORY('SubIndex', extns=list('abc'))
        _idx2 = self.FACTORY('SubIndex2', extns=list('abc'), parent=_idx1)
        self.assertIsInstance(_idx2(), index.IndexType)
        self.assertIsInstance(_idx2(), _idx2)

    def test_returns_false_for_child_types(self):
        _idx1 = self.FACTORY('SubIndex', extns=list('abc'))
        _idx2 = self.FACTORY('SubIndex2', extns=list('abc'), parent=_idx1)
        self.assertNotIsInstance(_idx1(), _idx2)
