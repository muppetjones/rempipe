
import unittest


class TestSandbox(unittest.TestCase):

    def test_dynamic_class_attribute_recall(self):

        class DynamicClass(object):

            def __init__(self):
                self.val = 'a'

            def dyn_meth1(self):
                return self.val

            def dyn_meth2(self):
                return lambda: self.val

        dc = DynamicClass()
        dm1 = dc.dyn_meth1()
        dm2 = dc.dyn_meth2()
        dm3 = dc.dyn_meth1

        dc.val = 'b'

        self.assertEqual(dc.val, 'b')  # positive control (self)
        self.assertNotEqual(dm1, 'b')  # negative control (immutable string)
        self.assertEqual(dm2(), 'b')  # lambda test
        self.assertEqual(dm3(), 'b')  # dynamic function test
