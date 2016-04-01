
import unittest

from tests import base

# from libpipe.util import path
from libpipe.type import base as _type

import logging
log = logging.getLogger(__name__)


class TypeBaseTestCase(base.LibpipeTestCase):

    @classmethod
    def setUpClass(cls):
        if cls is TypeBaseTestCase:
            raise unittest.SkipTest("Skip BaseTest tests, it's a base class")
        super(TypeBaseTestCase, cls).setUpClass()


class TestTypeBase(TypeBaseTestCase):

    def setUp(self):
        super().setUp()
        self.TYPE = _type.TypeBase

    def test_factory_returns_a_child_TypeBase_class(self):

        parent = self.TYPE()
        child_class = self.TYPE.factory()
        child = child_class()

        self.assertNotIsInstance(parent, child_class)
        self.assertIsInstance(child, self.TYPE)

    def test_factory_child_classes_are_unique(self):
        first_child = self.TYPE.factory()
        first_child_obj = first_child()

        second_child = self.TYPE.factory()
        self.assertIsInstance(first_child_obj, first_child)
        self.assertNotIsInstance(first_child_obj, second_child)

    def test_factory_sets_child_name(self):
        child_class = self.TYPE.factory(name='Dolly')
        self.assertEqual(child_class.__name__, 'Dolly')
