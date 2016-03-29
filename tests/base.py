import builtins
import unittest

from unittest import mock


class LibpipeTestCase(unittest.TestCase):

    def setup_mock_write(self, create=False, side_effect=None):
        if create:
            m = mock.mock_open()
        else:
            try:
                m = self._mock_open
            except AttributeError:
                m = mock.mock_open()
                self._mock_open = m

        patcher = mock.patch.object(
            builtins, 'open', m, create=True)
        m = patcher.start()
        self.addCleanup(patcher.stop)
        if side_effect:
            m.side_effect = side_effect
        return m

    def setup_mock_read(self, read_data, side_effect=None):
        patcher = mock.patch.object(
            builtins, 'open',
            mock.mock_open(read_data=read_data),
            create=True)
        m = patcher.start()
        self.addCleanup(patcher.stop)
        if side_effect:
            m.side_effect = side_effect
        m.return_value.__iter__ = lambda s: s
        m.return_value.__next__ = lambda s: s.readline()
        return m
