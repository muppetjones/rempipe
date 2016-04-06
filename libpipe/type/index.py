
import os.path

from libpipe.type import base
from libpipe.util import path


import logging
log = logging.getLogger(__name__)


def factory(name='IndexSubType', extns=[], counts=[]):
    if not counts:
        counts = [1] * len(extns)  # assume only one of each extn required!
    elif counts and len(extns) != len(counts):
        msg = 'Unequal num of extns ({}) vs counts ({})'.format(
            len(extns), len(counts))
        raise ValueError(msg)
    _subcls = type(name, (str, IndexType, ), dict(extns=extns, counts=counts))
    return _subcls


class IndexMeta(base.TypeMeta):
    pass


class IndexType(metaclass=IndexMeta):

    '''Alignment index type, e.g., hisat2 index

    Stores and validates a given path as an alignment index.

    Attributes:
        extns: A list of expected extensions for the index.
        counts: A list of expected file counts for each extension.
        path: The directory where the index is found.
        name: The basename of the index.

    Examples:
        # For a genome with a hisat2 index named GRCh38p5, i.e., GRCh38p5.*.ht2
        index_class = IndexType.factory({'.ht2': 8}, name="Hisat2Index")
        index_obj = index_class("path/to/genome/GRCh38p5")
        print(index_obj)
        print(index_obj.path)
        print(index_obj.name)

        > path/to/genome/GRCh38p5
        > path/to/genome
        > GRCh38p5
    '''

    def __init__(self, value=None):

        if value:
            self.path, self.name = os.path.split(value)
            self._check_extns()

    def _check_extns(self):
        files = path.walk_file(
            self.path, extension=self.extns, pattern=[self.name])

        errs = []
        for extn, cnt in zip(self.extns, self.counts):
            matched = [f for f in files if f.endswith(extn)]
            if len(matched) != cnt:
                errs.append('{}: {} v {}'.format(extn, len(matched), cnt))

        if errs:
            msg = 'Bad index: {}'.format('; '.join(errs))
            raise ValueError(msg)
