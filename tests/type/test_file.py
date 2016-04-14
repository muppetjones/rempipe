
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

    def setUp(self):
        super().setUp()

        self.PARENT = _file.FileType
        self.FACTORY = _file.factory

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


class TestFileType__TypeBase(test_base.TestTypeBase):

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

    def test_FileType_cannot_be_called_directly(self):
        with self.assertRaises(TypeError):
            _file.FileType('path/to/nowhere')


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


class TestSubType__register(FileTypeTestCase):

    def setUp(self):
        super().setUp()

        self.FACTORY = partial(self.FACTORY, extns=['.par'])

    def test_registry_shared_across_all_class(self):
        i1 = self.FACTORY(name='I_1')
        i2 = self.FACTORY(name='I_2')
        ii1 = self.FACTORY(name='II_1', parent=i1)
        ii2 = self.FACTORY(name='II_2', parent=i2)

        ii2.registry.update({'brave': 'world'})
        clses = [_file.FileType, i1, i2, ii1, ii2]
        for cls in clses:
            with self.subTest(class_test=cls):
                self.assertIn('brave', cls.registry)
                self.assertEqual(cls.registry['brave'], 'world')

    def test_child_class_added_to_registry(self):
        name = 'FauxPas'
        cls = self.FACTORY(name=name)
        name = name.lower()
        self.assertIn(name, cls.registry)
        self.assertEqual(cls.registry[name], cls)

    def test_grandchild_class_added_to_registry(self):
        name = 'FauxPas'
        cls = self.FACTORY(name='Meh')
        subcls = self.FACTORY(name=name, parent=cls)
        name = name.lower()
        self.assertIn(name, cls.registry)
        self.assertEqual(cls.registry[name], subcls)

    def test_existing_name_raises_ValueError(self):
        name = 'FauxPas'
        _file.FileMeta(name, (_file.FileType,), dict())
        with self.assertRaisesRegex(ValueError, '[Dd]uplicate'):
            _file.FileMeta(name, (_file.FileType,), dict())

    def test_only_FileMeta_types_added_to_registry(self):
        name = 'FauxPas'
        cls = self.FACTORY(name=name)
        self.assertNotIn('str', cls.registry)
        for clsobj in cls.registry.values():
            self.assertIsInstance(clsobj, _file.FileMeta)

    def test_child_registry_initialized_with_no_children(self):
        name = 'FauxPas'
        cls = self.FACTORY(name=name)
        name = name.lower()
        self.assertIn(name, cls.children)
        self.assertEqual(cls.children[name], [])

    def test_child_adds_self_name_to_immediate_parental_child_lists(self):
        child = self.FACTORY(name='ChildClass')
        gchild = self.FACTORY(name='GchildClass', parent=child)

        gchild_name = gchild.__name__.lower()
        child_name = child.__name__.lower()
        parent_name = _file.FileType.__name__.lower()

        self.assertIn(gchild_name, child.children[child_name])
        self.assertNotIn(gchild_name, child.children[parent_name])

    def test_non_file_type_not_stored_in_child_registry(self):
        name = 'FauxPas'
        cls = self.FACTORY(name=name)
        self.assertNotIn('str', cls.children)

        # test that all of the parents in the child registry
        # are also in the main registry ()
        children_keys = sorted(cls.children.keys())
        registry_keys = sorted(cls.registry.keys())
        self.assertEqual(registry_keys, children_keys)

    def test_get_children_returns_list_of_self_children_classes(self):
        subcls = self.FACTORY()
        subsub = self.FACTORY(name='letslipthedogsofwar', parent=subcls)

        self.assertEqual(self.PARENT.get_children(), [self.CHILD, subcls])
        self.assertEqual(subcls.get_children(), [subsub])
        self.assertEqual(subsub.get_children(), [])

# __END__
