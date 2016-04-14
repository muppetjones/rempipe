

from libpipe.type import file as _file
from libpipe.type import seq

from tests import base


class TestSeqType(base.LibpipeTestCase):

    def test_SeqType_type_is_FileMeta(self):
        self.assertEqual(type(seq.SeqType), _file.FileMeta)


class TestFastqType(base.LibpipeTestCase):

    def test_FastqType_is_SeqType(self):
        self.assertTrue(issubclass(seq.FastqType, seq.SeqType))
