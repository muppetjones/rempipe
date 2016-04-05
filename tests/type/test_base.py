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


class TypeBaseTestCase(base.LibpipeTestCase):
    pass


class TestTypeBase(TypeBaseTestCase):

    @classmethod
    def setUpClass(cls):
        if cls is TestTypeBase:
            raise unittest.SkipTest("Skip BaseTest tests, it's a base class")
            super(TestTypeBase, cls).setUpClass()

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
        foo = self.FACTORY('First')
        bar = self.FACTORY('First')
        self.assertNotEqual(foo, bar)

    def test_factory_classes_create_unrelated_objects(self):
        foo = self.FACTORY('Foo')
        bar = self.FACTORY('Bar')

        fobj = foo('foo')
        bobj = bar('bar')

        self.assertIsInstance(bobj, bar)
        self.assertIsInstance(fobj, foo)
        self.assertNotIsInstance(bobj, foo)
        self.assertNotIsInstance(fobj, bar)

    def test_factory_sets_child_name(self):
        child_class = self.FACTORY(name='Dolly')
        self.assertEqual(child_class.__name__, 'Dolly')


class TestTypeBase_BaseOnly(TestTypeBase):

    def setUp(self):
        super().setUp()
        self.META = _type.TypeMeta
        self.PARENT = _type.TypeBase
        self.FACTORY = _type.factory
        self.CHILD = self.FACTORY('ChildClass')
