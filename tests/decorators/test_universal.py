
import inspect


from textwrap import dedent
from unittest import TestCase
from unittest.mock import MagicMock

import libpipe.decorators.universal as unidec
# import remsci.lib.decorators.base as decor
# from remsci.lib.decorators import universal_function_decorator

import logging
log = logging.getLogger(__name__)


call_tester = MagicMock()


class TestDecoratorFactoryOnFunctions(TestCase):

    # Setup
    # -------------------------------------------------------------------------

    def setUp(self):
        call_tester.reset_mock()

    # Helper methods
    # -------------------------------------------------------------------------

    def get_wrapped_function(self):
        @unidec.universal_function_decorator
        def test_descriptor(f, args, kwargs):
            '''this is the decorated decorator'''
            call_tester.start()
            rv = f(*args, **kwargs)
            call_tester.end()
            return rv

        @test_descriptor
        def wrapped_function(*args, **kwargs):
            '''this is the wrapped function'''
            call_tester.wrapped_call(*args, **kwargs)
            return "that's the way I need it"

        return wrapped_function

    # Positive Tests
    # -------------------------------------------------------------------------

    def test_factory_called_correctly(self):
        wf = self.get_wrapped_function()
        wf()

        self.assertEqual(1, call_tester.start.call_count)
        self.assertEqual(1, call_tester.end.call_count)

    def test_function_called_correctly(self):
        wf = self.get_wrapped_function()
        wf("hello", "world", answer=42)
        call_tester.wrapped_call.assert_called_once_with(
            "hello", "world", answer=42)

    def test_basic_decorator_preservation_on_function(self):
        wf = self.get_wrapped_function()

        self.assertEqual('wrapped_function', wf.__name__)
        self.assertEqual('this is the wrapped function', wf.__doc__)

    def test_basic_decorator_function_return_value(self):
        wf = self.get_wrapped_function()
        rv = wf('any', 'way', 'you', 'want', 'it')

        self.assertEqual(rv, "that's the way I need it")

    def test_basic_decorator_preserves_inspect_values(self):
        wf = self.get_wrapped_function()

        inspect_src = dedent(inspect.getsource(wf)).rstrip().lstrip()
        expected_src = dedent("""
            @test_descriptor
            def wrapped_function(*args, **kwargs):
                '''this is the wrapped function'''
                call_tester.wrapped_call(*args, **kwargs)
                return "that's the way I need it"
        """).rstrip().lstrip()

        self.assertEqual(inspect_src, expected_src)


class TestDecoratorFactoryOnMethods(TestCase):

    # Setup
    # -------------------------------------------------------------------------

    def setUp(self):
        call_tester.reset_mock()

    # Helper methods
    # -------------------------------------------------------------------------

    def get_wrapped_method(self):
        @unidec.universal_function_decorator
        def test_descriptor(f, args, kwargs):
            '''this is the decorated decorator'''
            call_tester.start()
            rv = f(*args, **kwargs)
            call_tester.end()
            return rv

        class Classy(object):

            @test_descriptor
            def wrapped_method(self, *args, **kwargs):
                '''this is the wrapped method'''
                call_tester.wrapped_call(*args, **kwargs)
                return "that's the way I need it"

        return Classy().wrapped_method

    # Positive Tests
    # -------------------------------------------------------------------------

    def test_factory_called_correctly(self):
        wm = self.get_wrapped_method()
        wm()

        self.assertEqual(1, call_tester.start.call_count)
        self.assertEqual(1, call_tester.end.call_count)

    def test_function_called_correctly(self):
        wm = self.get_wrapped_method()
        wm("hello", "world", answer=42)
        call_tester.wrapped_call.assert_called_once_with(
            "hello", "world", answer=42)

    def test_basic_decorator_preservation_on_method(self):
        wm = self.get_wrapped_method()

        self.assertEqual('wrapped_method', wm.__name__)
        self.assertEqual('this is the wrapped method', wm.__doc__)

    def test_basic_decorator_function_return_value(self):
        wm = self.get_wrapped_method()
        rv = wm('any', 'way', 'you', 'want', 'it')

        self.assertEqual(rv, "that's the way I need it")

    def test_basic_decorator_preserves_inspect_values(self):
        wm = self.get_wrapped_method()

        inspect_src = dedent(inspect.getsource(wm)).rstrip().lstrip()
        expected_src = dedent("""
            @test_descriptor
            def wrapped_method(self, *args, **kwargs):
                '''this is the wrapped method'''
                call_tester.wrapped_call(*args, **kwargs)
                return "that's the way I need it"
        """).rstrip().lstrip()

        self.assertEqual(inspect_src, expected_src)
