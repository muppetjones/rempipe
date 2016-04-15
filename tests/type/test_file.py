
import copy

from functools import partial
from unittest import mock

# from libpipe.util import path
from libpipe.type import base as _type
from libpipe.type import file as _file
from tests import base
from tests.type import test_base

import logging
log = logging.getLogger(__name__)


class FileTypeTestCase(base.LibpipeTestCase):

    '''Define default setup and helper functions'''

    @classmethod
    def setUpClass(cls):
        cls.old_registry = copy.deepcopy(_file.FileType.registry)
        cls.old_children = copy.deepcopy(_file.FileType.children)

    @classmethod
    def tearDownClass(cls):
        # Reload saved registries after each test
        _file.FileType.registry = copy.deepcopy(cls.old_registry)
        _file.FileType.children = copy.deepcopy(cls.old_children)

    def setUp(self, child=True):
        super().setUp()

        self.META = _file.FileMeta
        self.PARENT = _file.FileType
        self.FACTORY = _file.factory
        if child:
            self.CHILD = self.FACTORY(
                name='TestFile', extns=['.bar', '.foo'])

    def tearDown(self):
        # Reload saved registries after each test
        _file.FileType.registry = copy.deepcopy(self.old_registry)
        _file.FileType.children = copy.deepcopy(self.old_children)


class TestFileMeta(base.LibpipeTestCase):

    '''Tests against the meta directly'''

    def test_FileMeta_inherits_from_TypeMeta(self):
        self.assertTrue(
            issubclass(_file.FileMeta, _type.TypeMeta),
            "FileMeta is NOT subclassed from TypeMeta",
        )


class TestFileType__TypeBase(base.LibpipeTestCase, test_base.TestTypeBase):

    '''Test against FileType directly (inherits from TypeBase tests)'''

    @classmethod
    def setUpClass(cls):
        cls.old_registry = copy.deepcopy(_file.FileType.registry)
        cls.old_children = copy.deepcopy(_file.FileType.children)

    @classmethod
    def tearDownClass(cls):
        # Reload saved registries after each test
        _file.FileType.registry = copy.deepcopy(cls.old_registry)
        _file.FileType.children = copy.deepcopy(cls.old_children)

    def setUp(self):

        super().setUp()

        # define variables used by TypeBase tests
        # -- note the partial factory declaration!
        self.META = _file.FileMeta
        self.PARENT = _file.FileType
        self.FACTORY = partial(_file.factory, extns=['.foo'])
        self.CHILD = self.FACTORY('FileFooType')  # includes extns=['.foo']

        # We don't need to check extensions here
        patcher = mock.patch('libpipe.type.file.FileType._check_extns')
        self.mock_check_extns = patcher.start()
        self.addCleanup(patcher.stop)

    def tearDown(self):
        # Reload saved registries after each test
        _file.FileType.registry = copy.deepcopy(self.old_registry)
        _file.FileType.children = copy.deepcopy(self.old_children)

    def test_child_is_still_TypeMeta(self):
        self.assertIsInstance(self.CHILD, _type.TypeMeta)


class TestFileType(FileTypeTestCase):

    def setUp(self):
        super().setUp(child=False)

    def test_class_wo_any_valid_extns_raises_TypeError(self):
        child = self.META('ChildWoExtn', (self.PARENT, ), dict())
        for cls in [self.PARENT, child]:
            with self.subTest(no_extn=cls):
                with self.assertRaises(TypeError):
                    cls('path/to/know.where')

    def test_class_wo_extns_accepts_any_immediate_child_extn(self):
        child = self.META('ChildWoExtn', (self.PARENT, ), dict())
        gchild1 = _file.factory(extns=['.foo'], parent=child)
        _file.factory(extns=['.bar'], parent=child)
        _file.factory(extns=['.baz'], parent=gchild1)

        valid_extns = ['.foo', '.bar']
        invalid_extns = ['.baz']
        for extn in valid_extns:
            with self.subTest(valid_extn=extn):
                # should not raise!
                child('path/to/file{}'.format(extn))
        for extn in invalid_extns:
            with self.subTest(invalid_extn=extn):
                with self.assertRaises(ValueError):
                    child('path/to/file{}'.format(extn))


class TestFileType__factory(FileTypeTestCase):

    '''Test against FileType factory'''

    #
    #   Factory--returned class
    #

    def test_factory_sets_child_name(self):
        child_class = self.FACTORY(name='Dolly', extns=list('ab'))
        self.assertEqual(child_class.__name__, 'Dolly')

    def test_factory_creates_custom_name_via_extns_if_not_given(self):
        # also tests removal of leading periods
        child_class = self.FACTORY(extns=['.a', 'b'])
        self.assertTrue(child_class.__name__.endswith('_a_b'))

    def test_factory_return_is_subclassed_from_str_by_default(self):
        subtype = self.FACTORY(extns=list('ab'))
        self.assertTrue(issubclass(subtype, str), 'Not subclassed from str')

    def test_factory_adds_extns_attribute_to_class(self):
        subtype = self.FACTORY(extns=list('ab'))
        self.assertEqual(subtype.extns, list('ab'))

    def test_factory_converts_extns_to_lowercase(self):
        subtype = self.FACTORY(extns=list('AB'))
        self.assertEqual(subtype.extns, list('ab'))

    def test_factory_class_is_subtype_of_given_parent(self):
        subtype = self.FACTORY(extns=list('ab'))
        subsub = self.FACTORY(extns=list('ac'), parent=subtype)
        self.assertTrue(issubclass(subsub, subtype), 'Parent not set')

    def test_factory_raises_TypeError_if_non_Meta_parent_given(self):

        types = [int, str, type('Faux', (object,), dict())]
        for test_type in types:
            with self.subTest(test_type=test_type):
                with self.assertRaises(TypeError):
                    self.FACTORY(extns=list('ab'), parent=test_type)

    def test_factory_raises_ValueError_if_no_extns_given(self):
        with self.assertRaises(ValueError):
            self.FACTORY()

    def test_factory_class_is_str_subclass_by_default(self):
        subtype = self.FACTORY(extns=list('ab'))
        self.assertTrue(issubclass(subtype, str), 'Factory class is not str')

    def test_raises_ValueError_if_name_already_exists_w_diff_extn(self):
        self.FACTORY(name='HelloWorld', extns=list('ab'))
        with self.assertRaises(ValueError):
            self.FACTORY(name='HelloWorld', extns=list('cd'))

    def test_returns_existing_class_if_same_name_and_extn_given(self):
        a = self.FACTORY(name='HelloWorld', extns=list('ab'))
        b = self.FACTORY(name='HelloWorld', extns=list('ab'))
        self.assertEqual(a, b)


class TestSubType(FileTypeTestCase):

    '''Test children types created by the factory'''

    def test_init_with_no_args_returns_empty_string(self):
        self.assertEqual(self.CHILD(), '')

    def test_object_is_set_to_given_value(self):
        obj = self.CHILD('foo/file.bar')
        self.assertEqual(obj, 'foo/file.bar')

    def test_init_calls_check_extn_with_given_value(self):

        with mock.patch.object(self.CHILD, '_check_extns') as m:
            self.CHILD('foo/file.bar')
        m.assert_called_once_with('foo/file.bar')

    def test_init_raises_ValueError_if_value_does_not_have_extn(self):
        with self.assertRaises(ValueError):
            self.CHILD('foo/file.42')

    def test_init_correctly_accepts_multi_period_extensions(self):
        self.CHILD.extns.append('.3.bt2')
        self.CHILD('foo/file.1.54.foo')  # does not raise -- matches .foo
        self.CHILD('foo/file.3.bt2')  # does not raise -- matches .3.bt2
        with self.assertRaises(ValueError):
            self.CHILD('foo/file.4.bt2')  # requires complete match to .3.bt2

    def test_check_extns_is_case_insensitive(self):
        self.CHILD('foo/file.FOO')  # should not raise

    def test_invalid_extn_err_is_informative(self):
        try:
            self.CHILD('foo/file.42')
        except ValueError as e:
            err_msg = str(e)

        test_msg = {
            '42': 'Msg does not contain bad extension',
            self.CHILD.__name__: 'Msg does not contain offended class name',
        }
        test_msg.update({
            v: 'Msg does not contain expected extension "{}"'.format(v)
            for v in self.CHILD.extns
        })
        for k, v in test_msg.items():
            with self.subTest(test=k):
                self.assertIn(k, err_msg, v)


class TestSubType__eq(FileTypeTestCase):

    '''Test children type equality

    Default __eq__ should have correct functionality.
    '''

    def test_obj_eq_given_str(self):
        expected = 'foo/file.foo'
        obj = self.CHILD(expected)
        self.assertTrue(obj == expected, 'Obj does not equal give str.')

    def test_obj_ne_to_same_prefix_str_w_other_valid_extn(self):
        expected = 'foo/file.bar'
        obj = self.CHILD(expected.replace('.bar', '.foo'))
        self.assertEqual(obj, expected.replace('.bar', '.foo'))
        self.assertTrue(obj != expected, 'Obj should not equal valid extn.')

    def test_sister_obj_ne_just_bc_same_extn(self):
        expected = 'foo/file.foo'
        obj1 = self.CHILD(expected)
        obj2 = self.CHILD(expected.replace('file', 'seq'))
        self.assertTrue(obj1 != obj2, 'Diff sister objects should be !=')

    def test_sister_obj_are_eq_if_same_str(self):
        expected = 'foo/file.foo'
        obj1 = self.CHILD(expected)
        obj2 = self.CHILD(expected)
        self.assertTrue(obj1 == obj2, 'Sister objects are not equal')

    def test_sister_obj_are_ne_if_same_prefix_w_valid_extn(self):
        expected = 'foo/file.foo'
        obj1 = self.CHILD(expected)
        obj2 = self.CHILD(expected.replace('.foo', '.bar'))
        self.assertTrue(obj1 != obj2, 'Sister objects are equal')

    def test_cousin_obj_are_eq_if_same_str(self):
        expected = 'foo/file.foo'
        obj1 = self.CHILD(expected)
        cls2 = self.FACTORY(name='Faux1', extns=['.foo', '.fq'])
        obj2 = cls2(expected)
        self.assertTrue(obj1 == obj2, 'Cousin objects are not equal')

    def test_cousin_obj_are_ne_if_same_prefix_w_valid_extn(self):
        expected = 'foo/file.foo'
        obj1 = self.CHILD(expected)
        cls2 = self.FACTORY(name='Faux1', extns=['.foo', '.fq'])
        obj2 = cls2(expected.replace('.foo', '.fq'))
        self.assertTrue(obj1 != obj2, 'Cousin objects are equal')


class TestSubType__instancecheck(FileTypeTestCase):

    '''Test children type isinstance'''

    def test_obj_is_self_type(self):
        obj = self.CHILD()
        self.assertIsInstance(obj, self.CHILD)

    def test_obj_is_base_parent_type(self):
        obj = self.CHILD()
        self.assertIsInstance(obj, self.PARENT)

    def test_obj_is_not_sister_type(self):
        obj = self.CHILD()
        cls2 = self.FACTORY(name='Faux1', extns=['.foo', '.fq'])
        self.assertNotIsInstance(obj, cls2)

    def test_grandchild_is_parent_types(self):
        # subtle check for inf recursion, too
        subsub = self.FACTORY(name='Faux', extns=['.42'], parent=self.CHILD)
        subobj = subsub()
        self.assertIsInstance(subobj, subsub)  # should not raise
        self.assertIsInstance(subobj, self.CHILD)  # should not raise
        self.assertIsInstance(subobj, self.PARENT)  # should not raise

    def test_str_w_valid_extn_returns_isinstance_true(self):
        self.assertIsInstance('hello_my_baby.foo', self.CHILD)

    def test_str_w_invalid_extn_returns_isinstance_false(self):
        self.assertNotIsInstance('hello_my_baby', self.CHILD)


class TestSubType__register(
        FileTypeTestCase, test_base.TestTypeMixin__register):

    def setUp(self):
        super().setUp()
        self.FACTORY = partial(self.FACTORY, extns=['.par'])


# __END__
