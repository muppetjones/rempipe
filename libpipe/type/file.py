
import os.path

from libpipe.type import base

import logging
log = logging.getLogger(__name__)


def factory(name='SubFileType', extns=[], parent=None):

    if not extns:
        msg = 'No extensions given. Try factory(extns=[".fa", ...])'
        raise ValueError(msg)

    if not parent:
        parent = FileType

    if not type(parent) == FileMeta:
        msg = 'Parent type is not FileMeta'
        raise TypeError(msg)

    subcls = FileMeta(name, (parent, str), dict(
        extns=extns,
    ))
    return subcls


class FileMeta(base.TypeMeta):

    def __instancecheck__(self, instance):
        _check = super(self.__class__, self).__instancecheck__(instance)

        if not _check and type(instance) == str:
            try:
                self._check_extns(instance)
                _check = True
            except ValueError:
                pass  # already false...don't bother
        return _check


class FileType(metaclass=FileMeta):

    def __init__(self, value=''):

        if self.__class__ == FileType:
            msg = 'Cannot call base FileType directly.'
            raise TypeError(msg)

        if value:
            # Allowing an empty value matches other types, e.g.,
            # -- int() --> 0
            # -- str() --> ''
            self._check_extns(value)

    @classmethod
    def _check_extns(cls, value):
        '''Check given value for the presence of an expected extension.

        Instead of using os.path.splitext, use str.endswith. This allows
        for simpler checking given more complex extension types,
        e.g., *.rev.bt2. The downside is looping through the extensions
        instead of checking if the extension is in the set.

        Alternate method: Split by periods. More fragile.

        Arguments:
            value: The value to check.
        Returns:
            None.
        Raises:
            ValueError if a valid extension not found.
        '''

        for extn in cls.extns:
            if value.endswith(extn):
                return  # found a valid extension
        msg = 'Invalid extension (file: {}) for type {} [{}]'.format(
            os.path.basename(value), cls.__name__, ', '.join(cls.extns))
        raise ValueError(msg)
