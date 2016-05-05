
from unittest import mock

from libpipe.cmd import interface
from libpipe.pipe import base
from libpipe.pipe import batch
from tests.base import LibpipeTestCase  # includes read and write mock


class BatchPipeTestCase(LibpipeTestCase):

    pass

#
#   Pipe
#


class TestBatchPipe(BatchPipeTestCase):

    '''Test that BatchPipe still acts like a pipe'''

    def test_inherits_from_Pipe(self):

        pipe = batch.BatchPipe()
        self.assertIsInstance(pipe, base.Pipe)

    def test_init_calls_super_init(self):
        self.fail()

    def test_init_raises_ValueError_if_not_given_input_dict(self):
        '''Test that init complains if not given kwarg input as dict'''

        tests = ['filename.txt', ['a.txt', 'b.txt'], 0]
        for test in tests:
            with self.subTest(input=test):
                with self.assertRaises(ValueError):
                    batch.BatchPipe(input=test)

#
#   CmdInterface
#


class TestBatchPipe_CmdInterface(BatchPipeTestCase):

    '''Test that the BatchPipe implements the CmdInterface as expected'''

    def test_inherits_from_cmd_interface(self):
        pipe = batch.BatchPipe()
        self.assertIsInstance(pipe, interface.CmdInterface)

    def test_link_raises_TypeError_when_called(self):
        pipe = batch.BatchPipe()
        m = mock.MagicMock('libpipe.cmd.base.CmdBase')
        with self.assertRaises(TypeError):
            pipe.link(m)
