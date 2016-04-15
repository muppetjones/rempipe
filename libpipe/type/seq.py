
from libpipe.type import file as _file


class SeqType(_file.FileType):
    pass

FastqType = _file.factory(
    name='FastqType',
    parent=SeqType,
    extns=['.fq', '.fastq'],
)

FastaType = _file.factory(
    name='FastaType',
    parent=SeqType,
    extns=['.fa', '.fasta', '.seq', '.fna', '.faa', '.ffn'],
)


class SeqMapType(_file.FileType):
    pass

BamType = _file.factory(
    name='BamType',
    parent=SeqMapType,
    extns=['.bam', ],
)

SamType = _file.factory(
    name='SamType',
    parent=SeqMapType,
    extns=['.sam', ],
)
