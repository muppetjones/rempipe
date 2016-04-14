'''Test metaclass type

Note the following naming is used:
    meta > parent > child > obj
'''


import unittest

from tests import base

# from libpipe.util import path
from libpipe.type import base as _type

import logging
log = logging.getLogger(__name__)


class TestTypeBase():  # base.LibpipeTestCase):

    '''A set of restrictive tests for metaclass factories

    These tests are intentionally over-restrictive. Metaclasses
    are nasty, trick-sy little things. Ease up later.
    '''

    def print_summary(self):
        meta = self.META
        parent = self.PARENT
        factory = self.FACTORY
        child = self.CHILD
        child_obj = child('path/to/something')

        log.debug('meta:    {}'.format(type(meta)))
        log.debug('parent:  {}'.format(type(parent)))
        log.debug('factory: {}'.format(type(factory)))
        log.debug('child:   {}'.format(type(child)))
        log.debug('object:  {}'.format(type(child_obj)))

        log.debug('Parent is Meta?   {}'.format(isinstance(parent, meta)))
        log.debug('Child is Meta?    {}'.format(isinstance(child, meta)))
        log.debug('Object is Meta?   {}'.format(isinstance(child_obj, meta)))
        log.debug('Child is Parent?  {}'.format(isinstance(child, parent)))
        log.debug('Object is Parent? {}'.format(isinstance(child_obj, parent)))
        log.debug('Object is Child?  {}'.format(isinstance(child_obj, child)))

    def test_the_parent_is_instance_of_meta(self):
        self.assertIsInstance(self.PARENT, self.META)

    def test_factory_returns_instance_of_meta(self):
        self.assertIsInstance(self.FACTORY('ChildClass'), self.META)

    def test_factory_return_is_not_instance_of_parent(self):
        self.assertNotIsInstance(self.FACTORY('ChildClass'), self.PARENT)

    def test_child_obj_is_not_instance_of_meta(self):
        obj = self.CHILD('value')
        self.assertNotIsInstance(obj, self.META)

    def test_child_obj_is_instance_of_parent(self):
        obj = self.CHILD('value')
        self.assertIsInstance(obj, self.PARENT)

    def test_child_obj_is_instance_of_child(self):
        obj = self.CHILD('value')
        self.assertIsInstance(obj, self.CHILD)

    def test_factory_classes_are_unique(self):
        foo = self.FACTORY('Foo')
        bar = self.FACTORY('Bar')
        self.assertNotEqual(foo, bar)

    def test_subsequent_factory_calls_return_different_subclasses(self):
        foo = self.FACTORY('Foo')
        bar = self.FACTORY('Bar')
        self.assertEqual(type(foo), type(bar))  # Same type
        self.assertFalse(isinstance(bar, foo))  # Not inherited
        self.assertNotEqual(foo, bar)  # Not the same class

    def test_subseq_factory_classes_create_sister_objects(self):
        foo = self.FACTORY('Foo')
        bar = self.FACTORY('Bar')

        fobj = foo('foo')
        bobj = bar('bar')

        # NOTE: both obj are same super-type, i.e., IndexMeta
        self.assertTrue(isinstance(bobj, bar))
        self.assertTrue(isinstance(fobj, foo))
        self.assertFalse(isinstance(bobj, foo))
        self.assertFalse(isinstance(fobj, bar))

    def test_factory_sets_child_name(self):
        child_class = self.FACTORY(name='Dolly')
        self.assertEqual(child_class.__name__, 'Dolly')


class TestTypeMixin__register(object):

    def test_registry_shared_across_all_class(self):
        i1 = self.FACTORY(name='I_1')
        i2 = self.FACTORY(name='I_2')
        ii1 = self.FACTORY(name='II_1', parent=i1)
        ii2 = self.FACTORY(name='II_2', parent=i2)

        ii2.registry.update({'brave': 'world'})
        clses = [self.PARENT, i1, i2, ii1, ii2]
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
        self.META(name, (self.PARENT,), dict())
        with self.assertRaisesRegex(ValueError, '[Dd]uplicate'):
            self.META(name, (self.PARENT,), dict())

    def test_only_FileMeta_types_added_to_registry(self):
        name = 'FauxPas'
        cls = self.FACTORY(name=name)
        self.assertNotIn('str', cls.registry)
        for clsobj in cls.registry.values():
            self.assertIsInstance(clsobj, self.META)

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
        parent_name = self.PARENT.__name__.lower()

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

        parent_name = self.PARENT.__name__.lower()
        expected_parent = [
            self.old_registry[child]
            for child in self.old_children[parent_name]
        ] + [self.CHILD, subcls]
        self.assertEqual(self.PARENT.get_children(), expected_parent)
        self.assertEqual(subcls.get_children(), [subsub])
        self.assertEqual(subsub.get_children(), [])


class TestTypeBase_BaseOnly(base.LibpipeTestCase, TestTypeBase):

    def setUp(self):
        super().setUp()
        self.META = _type.TypeMeta
        self.PARENT = _type.TypeBase
        self.FACTORY = _type.factory
        self.CHILD = self.FACTORY('ChildClass')
