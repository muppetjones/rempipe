
from libpipe.type import file as _file


class AnnotationType(_file.FileType):
    pass

'''Generic Transfer Format (GTF)'''
GtfType = _file.factory(
    name='GtfType',
    parent=AnnotationType,
    extns=['.gtf', ],
)

'''Generic Feature Format (GFF)'''
GffType = _file.factory(
    name='GffType',
    parent=AnnotationType,
    extns=['.gff', '.gff3', '.gff2'],
)
