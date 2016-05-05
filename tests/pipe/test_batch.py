
from libpipe.cmd import interface
from libpipe.pipe import base
from libpipe.pipe import batch
from tests.base import LibpipeTestCase  # includes read and write mock


class BatchPipeTestCase(LibpipeTestCase):

    pass

#
#   Pipe
#


class TestBatchPipe_Pipe(BatchPipeTestCase):

    '''Test that BatchPipe still acts like a pipe'''

    def test_inherits_from_Pipe(self):

        pipe = batch.BatchPipe()
        self.assertIsInstance(pipe, base.Pipe)

#
#   CmdInterface
#


class TestBatchPipe_CmdInterface(BatchPipeTestCase):

    '''Test that the BatchPipe implements the CmdInterface as expected'''

    def test_inherits_from_cmd_interface(self):
        pipe = batch.BatchPipe()
        self.assertIsInstance(pipe, interface.CmdInterface)
