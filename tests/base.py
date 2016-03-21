import builtins
import unittest

from unittest import mock


class PipeBaseTest(unittest.TestCase):

    def setup_mock_write(self):
        patcher = mock.patch.object(
            builtins, 'open',
            mock.mock_open(),
            create=True)
        m = patcher.start()
        self.addCleanup(patcher.stop)
        return m

    def setup_mock_read(self, read_data):
        patcher = mock.patch.object(
            builtins, 'open',
            mock.mock_open(read_data=read_data),
            create=True)
        m = patcher.start()
        self.addCleanup(patcher.stop)
        m.return_value.__iter__ = lambda s: s
        m.return_value.__next__ = lambda s: s.readline()
        return m
