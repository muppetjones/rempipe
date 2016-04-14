
from libpipe.type import file as _file


class SeqType(_file.FileType):
    pass

FastqType = _file.factory(
    name='FastqType',
    parent=SeqType,
    extns=['.fq', '.fastq'],
)
