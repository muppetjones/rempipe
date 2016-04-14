

from libpipe.type import file as _file
from libpipe.type import seq

from tests import base


class TestSeqType(base.LibpipeTestCase):

    def test_SeqType_type_is_FileMeta(self):
        self.assertEqual(type(seq.SeqType), _file.FileMeta)


class TestFastqType(base.LibpipeTestCase):

    def test_FastqType_is_SeqType(self):
        self.assertTrue(issubclass(seq.FastqType, seq.SeqType))

    def test_FastqType_type_is_FileMeta(self):
        self.assertEqual(type(seq.FastqType), _file.FileMeta)

    def test_FastqType_accepts_valid_extensions(self):
        valid_extns = ['.fastq', '.fq', '.FASTQ', '.FQ', '.FastQ']
        for extn in valid_extns:
            with self.subTest(extn=extn):
                # should not raise!
                seq.FastqType('path/to/file{}'.format(extn))

    def test_FastqType_raises_ValueError_for_invalid_extensions(self):
        valid_extns = ['.fasq', '.fa', '.FASTA', '.bam']
        for extn in valid_extns:
            with self.subTest(extn=extn):
                with self.assertRaises(ValueError):
                    seq.FastqType('path/to/file{}'.format(extn))
