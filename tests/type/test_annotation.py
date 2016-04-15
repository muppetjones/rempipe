
from libpipe.type import annotation as anno
from libpipe.type import file as _file

from tests import base
from tests.type import custom_file_type as cft

import logging
log = logging.getLogger(__name__)


class TestGtfType(base.LibpipeTestCase, cft.TestFileTypeMixin):

    TYPE = anno.GtfType
    META = _file.FileMeta
    PARENT = anno.AnnotationType

    VALID_EXTNS = ['.gtf', '.GTF']
    INVALID_EXTNS = ['.gff', '.gff3', '.FASTA', '.bam']


class TestGffType(base.LibpipeTestCase, cft.TestFileTypeMixin):

    TYPE = anno.GffType
    META = _file.FileMeta
    PARENT = anno.AnnotationType

    VALID_EXTNS = ['.gff', '.GFF', '.gff3', '.gff2']
    INVALID_EXTNS = ['.gtf', '.fa']


class TestAnnoType(cft.TestSuperFileTypeTestCase, cft.TestFileTypeMixin):

    TYPE = anno.AnnotationType
    META = _file.FileMeta
    PARENT = _file.FileType
    CHILDREN = [TestGffType, TestGtfType]

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
