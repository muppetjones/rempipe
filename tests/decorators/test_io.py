import builtins
import io

import unittest
from unittest.mock import MagicMock, patch, mock_open

from libpipe.decorators.io import file_or_handle
from libpipe.util import path

import logging
log = logging.getLogger(__name__)


call_tester = MagicMock()


class TestSetup(unittest.TestCase):

    '''Setup mock_open for tests'''

    def setUp(self):
        self.mopen = self.setup_mock_write()

    def setup_mock_write(self):
        patcher = patch.object(builtins, 'open',
                               mock_open(),
                               create=True)
        m = patcher.start()
        self.addCleanup(patcher.stop)
        return m


class CommonTests(object):

    '''Tests common to 'file_or_handle' and 'file_or_handle_write' decorators

    NOTE: file_or_handle_write is not included in this distribution

    Children MUST implement:
        * get_decorated_function(self, *args, **kwargs)

    Included methods:
        * validate_open_called_with

    Includes tests:
        * [+] test_open_is_not_called_if_filehandle_given
        * [+] test_string_handle_successfully_written_to
        * [-] test_invalid_handle_object_raises_TypeError
        * [-] test_write_on_read_rasies_UnsupportedOperation
    '''

    # Helper methods
    # -------------------------------------------------------------------------

    def get_decorated_function(self, *args, **kwargs):
        raise NotImplementedError

    def validate_open_called_with(self, mode, *args, **kwargs):
        df = self.get_decorated_function(mode)
        filename = path.protect('~/test.txt')
        df(filename, *args, **kwargs)
        self.mopen.assert_called_once_with(filename, mode)

    # Positive
    # -------------------------------------------------------------------------

    def test_open_is_not_called_if_filehandle_given(self):
        df = self.get_decorated_function('w')

        filename = path.protect('~/test.txt')
        with open(filename, 'hahaha') as fh:
            df(fh)

        # Checking that a function is NOT called is very difficult.
        # So, we instead pass nonsensical data ('hahaha' in this case),
        # into the function we want to check and assert that it's only
        # called the ONE time the test calls it AND with the nonsensical
        # arguments.
        # NOTE: the filename is passed to the function, so we need to
        #       check the args NOT passed ('hahaha')
        self.mopen.assert_called_once_with(filename, 'hahaha')

    def test_string_handle_successfully_written_to(self):
        df = self.get_decorated_function('w')

        with io.StringIO() as sh:
            df(sh, 'writing to string!')
            sh.seek(0)
            result = sh.read()

        self.assertEqual(result, 'writing to string!')

    # Negative Tests
    # -------------------------------------------------------------------------

    def test_invalid_handle_object_raises_TypeError(self):
        class BadHandle(object):
            pass

        bh = BadHandle()
        df = self.get_decorated_function('w')

        # Raises here due to utility.path expecting a string.
        # Ideally, would fail due to an inability to open the object;
        # however, open would also raise a TypeError.
        with self.assertRaises(TypeError):
            df(bh, 'writing to string!')

    def test_write_on_read_raises_UnsupportedOperation(self):
        df = self.get_decorated_function('w')

        filename = path.protect('~/test.txt')
        with open(filename, 'r') as fh:

            # Sort-of-a-cheat. 'mock_open' will allow 'write' to be called
            #   on a file handle opened as read, so we have to force
            #   the error (which is what would be called anyway).
            # Interestingly, the unmocked handle is correctly interpreted
            #   as a handle (meaning it -does- have 'write').
            # NOTE: Could add 'writeable' check to decorator or check for
            #   a string type instead of 'write', but the current
            #   implementation works fine.
            fh.write.side_effect = io.UnsupportedOperation('not writable')
            with self.assertRaisesRegex(
                    io.UnsupportedOperation, 'not writable'):
                df(fh, 'writing to string!')


class TestFileOrHandle(TestSetup, CommonTests):

    '''Test 'file_or_handle' decorator

    Implements:
        * get_decorated_function(mode)

    Tests:
        * test_open_called_with_correct_mode
        * test_open_defaults_to_read_mode
    '''

    # Helper methods
    # -------------------------------------------------------------------------

    def get_decorated_function(self, mode):
        @file_or_handle(mode=mode)
        def decorated_function(fh, line=''):
            if line:
                fh.write(line)
            else:
                fh.read()
        return decorated_function

    # Positive Tests
    # -------------------------------------------------------------------------

    @unittest.skip('Patch python first--mock_issue18622_2')
    def test_open_called_with_correct_mode(self):
        for mode in list('wr+'):
            # self.mopen = self.setup_mock_write()
            self.mopen.reset_mock()
            with self.subTest(mode=mode):
                self.validate_open_called_with(mode)

    def test_open_defaults_to_read_mode(self):
        @file_or_handle
        def decorated_function(fh):
            fh.read()

        filename = path.protect('~/test.txt')
        decorated_function(filename)

        self.mopen.assert_called_once_with(filename, 'r')


class TestFileOrHandle_METHOD(TestFileOrHandle):

    '''Test 'file_or_handle_write' decorator on class methods'''

    def get_decorated_function(self, mode):
        class TempClass(object):

            @file_or_handle(mode=mode)
            def decorated_function(self, fh, line=''):
                if line:
                    fh.write(line)
                else:
                    fh.read()

        return TempClass().decorated_function

    def get_instance_with_decorated_function(self, *args, **kwargs):
        class ClassWithDecoratedMethod(object):

            def __init__(self):
                self.filename = None

            @file_or_handle
            def decorated_function(self, fh, line=''):
                fh.write(line)
        return ClassWithDecoratedMethod()

    def test_decorator_saves_filename_to_instance_if_hasattr(self):
        obj = self.get_instance_with_decorated_function()
        filename = path.protect('~/test.txt')
        obj.decorated_function(filename)

        self.assertEqual(obj.filename, filename)
