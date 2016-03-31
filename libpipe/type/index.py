

from libpipe.type import base


import logging
log = logging.getLogger(__name__)


class IndexType(base.TypeBase):

    def factory(*args, **kwargs):
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

        counts = None
        extns = None

        if len(args) > 1:
            extns = sorted(list(args))
        elif len(args) == 1:
            extns = args[0]

        if kwargs:
            extns = kwargs

        if isinstance(extns, dict):
            _extns = sorted(list(extns.keys()))
            counts = [extns[v] for v in _extns]
            extns = _extns
        else:
            try:
                extns = sorted(extns)
            except TypeError:
                pass  # we don't care if it's not iterable

        class IndexSubType(IndexType):
            expected_extns = extns
            expected_counts = counts

        return IndexSubType
