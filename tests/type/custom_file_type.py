
from tests import base


class TestSuperFileTypeTestCase(base.LibpipeTestCase):

    @classmethod
    def setUpClass(cls):

        cls.VALID_EXTNS = list(set(
            extn.lower()
            for child in cls.CHILDREN
            for extn in child.VALID_EXTNS
        ))
        cls.INVALID_EXTNS = list(set(
            extn.lower()
            for child in cls.CHILDREN
            for extn in child.INVALID_EXTNS
            if extn.lower() not in cls.VALID_EXTNS
        ))


class TestFileTypeMixin(object):

    def test_type_is_meta(self):
        self.assertEqual(type(self.TYPE), self.META)

    def test_subclassed_to_parent(self):
        self.assertTrue(
            issubclass(self.TYPE, self.PARENT),
            'Not subclassed to parent',
        )

    def test_accepts_valid_extensions(self):
        for extn in self.VALID_EXTNS:
            with self.subTest(extn=extn):
                # should not raise!
                self.TYPE('path/to/file{}'.format(extn))

    def test_raises_ValueError_for_invalid_extensions(self):
        for extn in self.INVALID_EXTNS:
            with self.subTest(extn=extn):
                with self.assertRaises(ValueError):
                    self.TYPE('path/to/file{}'.format(extn))
