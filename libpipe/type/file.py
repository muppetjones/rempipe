
from libpipe.type import base


def factory(name='SubFileType', extns=[]):
    subcls = FileMeta(name, (FileType, ), dict())
    return subcls


class FileMeta(base.TypeMeta):
    pass


class FileType(metaclass=FileMeta):

    def __init__(self, value):
        pass
