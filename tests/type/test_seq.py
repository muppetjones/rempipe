

from libpipe.type import file as _file
from libpipe.type import seq

from tests import base
from tests.type import custom_file_type as cft

import logging
log = logging.getLogger(__name__)


#
#   Base types
#


class TestFastqType(base.LibpipeTestCase, cft.TestFileTypeMixin):

    TYPE = seq.FastqType
    META = _file.FileMeta
    PARENT = seq.SeqType

    VALID_EXTNS = ['.fastq', '.fq', '.FASTQ', '.FQ', '.FastQ']
    INVALID_EXTNS = ['.fasq', '.fa', '.FASTA', '.bam']


class TestFastaType(base.LibpipeTestCase, cft.TestFileTypeMixin):

    TYPE = seq.FastaType
    META = _file.FileMeta
    PARENT = seq.SeqType

    VALID_EXTNS = ['.fasta', '.fa', '.fna', '.seq', '.faa', '.ffn']
    INVALID_EXTNS = ['.fastq', '.fq', '.gff', '.bam']


class TestBamType(base.LibpipeTestCase, cft.TestFileTypeMixin):

    TYPE = seq.BamType
    META = _file.FileMeta
    PARENT = seq.SeqMapType

    VALID_EXTNS = ['.bam', '.BAM']
    INVALID_EXTNS = ['.fastq', '.fq', '.gff', '.sam']


class TestSamType(base.LibpipeTestCase, cft.TestFileTypeMixin):

    TYPE = seq.SamType
    META = _file.FileMeta
    PARENT = seq.SeqMapType

    VALID_EXTNS = ['.sam', '.SAM']
    INVALID_EXTNS = ['.fastq', '.fq', '.gff', '.bam']

#
#   Super types
#


class TestSeqType(cft.TestSuperFileTypeTestCase, cft.TestFileTypeMixin):

    TYPE = seq.SeqType
    META = _file.FileMeta
    PARENT = _file.FileType
    CHILDREN = [TestFastqType, TestFastaType]

    @classmethod
    def setUpClass(cls):
        super().setUpClass()


class TestSeqMapType(cft.TestSuperFileTypeTestCase, cft.TestFileTypeMixin):

    TYPE = seq.SeqMapType
    META = _file.FileMeta
    PARENT = _file.FileType
    CHILDREN = [TestSamType, TestBamType]

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
