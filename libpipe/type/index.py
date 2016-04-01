
import os.path

from libpipe.type import base
from libpipe.util import path


import logging
log = logging.getLogger(__name__)


class IndexType(str, base.TypeBase):

    '''Alignment index type, e.g., hisat2 index

    Stores and validates a given path as an alignment index.

    Attributes:
        extns: A list of expected extensions for the index.
        counts: A list of expected file counts for each extension.
        _path: The directory where the index is found.
        _prefix: The "full" path to the index (path + name)
        _name: The basename of the index.

    Examples:
        # For a genome with a hisat2 index named GRCh38p5, i.e., GRCh38p5.*.ht2
        index_class = IndexType.factory({'.ht2': 8}, name="Hisat2Index")
        index_obj = index_class("path/to/genome/GRCh38p5")
        print(index_obj._path)
        print(index_obj._prefix)
        print(index_obj._name)

            > path/to/genome
            > path/to/genome/GRCh38p5
            > GRCh38p5
    '''

    counts = []
    extns = []

    def __init__(self, _path=None):
        self._path = _path
        self._prefix = None
        self._name = None

        self._check_prefix()
        self._check_extns()  # raise ValueError if unexpected file count!

    #
    #   magic methods
    #

    def __str__(self):
        return self._prefix or self._path

    def __eq__(self, other):
        return str(self) == other

    #
    #   Check methods
    #

    def _check_prefix(self):
        try:
            _is_full = os.path.isdir(self._path)
            _is_partial = os.path.isdir(os.path.dirname(self._path))
            if not _is_full and _is_partial:
                self._prefix = self._path
                self._path, self._name = os.path.split(self._prefix)
        except TypeError:
            pass  # probably given None

    def _check_extns(self):
        pattern = [self._name] if self._name is not None else []
        _files = path.walk_file(
            self._path, extension=self.extns, pattern=pattern)

        errs = []
        for extn, cnt in zip(self.extns, self.counts):
            _found = [f for f in _files if f.endswith(extn)]
            if len(_found) != cnt:
                errs.append('{}: {} v {}'.format(extn, len(_found), cnt))
        if errs:
            msg = 'Bad index: {}'.format('; '.join(errs))
            raise ValueError(msg)

    #
    #   fatory
    #

    def factory(*args, name=None, **kwargs):
        '''Create an IndexType with a defined set of extensions.

        Define and IndexType that expects a given number of the expected
        extensions.

        Example:
            index_class = IndexType.factory('.fq', '.fastq')
        Example:
            index_class = IndexType.factory(['.fq', '.fastq'])
        Example:
            index_class = IndexType.factory({'.bt2': 4, '.rev.bt2': 2})
        Example:
            index_class = IndexType.factory(**{'.bt2': 4, '.rev.bt2': 2})
        '''

        # TODO(sjbush): consider making a full metaclass

        _counts = []
        _extns = []

        if len(args) > 1:
            _extns = sorted(list(args))
        elif len(args) == 1:
            _extns = args[0]

        if kwargs:
            _extns = kwargs

        if isinstance(_extns, dict):
            _extns_sorted = sorted(list(_extns.keys()))
            _counts = [_extns[v] for v in _extns_sorted]
            _extns = _extns_sorted
        else:
            try:
                _extns = sorted(_extns)
            except TypeError:
                pass  # we don't care if it's not iterable

        if not _counts:
            _counts = [1] * len(_extns)

        name = name or 'IndexSubType'
        _subcls = type(name, (IndexType, ), dict(
            extns=_extns,
            counts=_counts
        ))

        return _subcls
