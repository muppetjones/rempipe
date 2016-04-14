'''Test the ability to cause cmd input to "fall_through" to output

NOTE: The input attribute is set on a Cmd instance during linking ONLY.
    Thus, the following tests attempt to simulate this by setting input
    on each instance using lambda.
'''


import unittest
from unittest import mock

from libpipe.cmd.base import CmdBase
from libpipe.decorators.cmd import fall_through
from libpipe.decorators.cmd import universal_fall_through

import logging
log = logging.getLogger(__name__)


class TestUniversalFallThroughDecorator(unittest.TestCase):

    '''Test the universal_fall_through decorator

    Due to the implementation, it really doesn't need to be universal.
    It should only ever be used on bound class methods.
    '''

    def setUp(self):
        patcher = mock.patch(
            'libpipe.cmd.attr.CmdAttributes',
        )
        self.mock_attr = patcher.start()
        self.mock_attr().defaults = {}
        self.addCleanup(patcher.stop)

        class FakeCmd(CmdBase):

            attr = self.mock_attr()

            @universal_fall_through
            def output(self):
                return ['hi', ]

        self.CMD = FakeCmd

    def test_universal_fallthrough_decorator(self):
        cmd = self.CMD()
        cmd.input = lambda: list('abc')
        expected = cmd.input() + ['hi', ]
        self.assertEqual(cmd.output(), expected)

    def test_raises_AttributeError_if_no_input_attribute(self):
        cmd = self.CMD()
        with self.assertRaises(AttributeError):
            cmd.output()


class TestBasicFallThroughDecorator(unittest.TestCase):

    def setUp(self):
        patcher = mock.patch(
            'libpipe.cmd.attr.CmdAttributes',
        )
        self.mock_attr = patcher.start()
        self.mock_attr().defaults = {}
        self.addCleanup(patcher.stop)

        class FakeCmd(CmdBase):

            attr = self.mock_attr()

            def output(self):
                return ['hi', ]

        self.CMD = FakeCmd

    def test_basic_fall_through_decorator(self):
        cmd = self.CMD()
        cmd.input = lambda: list('abc')
        cmd.output = fall_through(cmd.output)
        expected = cmd.input() + ['hi', ]
        self.assertEqual(cmd.output(), expected)

    def test_raises_AttributeError_if_no_input_attribute(self):
        cmd = self.CMD()
        cmd.output = fall_through(cmd.output)
        with self.assertRaises(AttributeError):
            cmd.output()
