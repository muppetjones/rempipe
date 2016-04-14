
import os.path

from libpipe.type import base

import logging
log = logging.getLogger(__name__)


def factory(name=None, extns=[], parent=None):

    if not extns:
        msg = 'No extensions given. Try factory(extns=[".fa", ...])'
        raise ValueError(msg)

    if not name:
        name = 'FileType_{}'.format('_'.join(
            extn if not extn.startswith('.') else extn[1:] for extn in extns
        ))

    if not parent:
        parent = FileType

    if not type(parent) == FileMeta:
        msg = 'Parent type is not FileMeta'
        raise TypeError(msg)

    try:
        subcls = FileMeta(name, (parent, str), dict(
            extns=extns,
        ))
    except ValueError:
        # Could be a class with the same name exists
        # -- if so, and the extensions are the same, just return the original
        possible_dup = FileType.registry[name.lower()]
        if extns == possible_dup.extns:
            subcls = possible_dup
        else:
            raise
    return subcls


class FileMeta(base.TypeMeta):

    def __init__(mcl, name, parents, dct):

        # add to registry -- check for existing first
        class_id = name.lower()

        try:
            if class_id in mcl.registry:
                msg = 'Duplicate FileType found: {}'.format(class_id)
                raise ValueError(msg)
        except AttributeError:
            mcl.registry = {}  # init main registry
        except ValueError:
            raise

        mcl.registry[class_id] = mcl

        # add self (or init) children registry
        try:
            mcl.children[class_id] = []
        except AttributeError:
            mcl.children = {class_id: []}  # init child registry

        # add self to parent registry
        for parent in parents:
            if type(parent) != FileMeta:
                continue
            parent_id = parent.__name__.lower()
            mcl.children[parent_id].append(class_id)

        super(FileMeta, mcl).__init__(name, parents, dct)

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

    @classmethod
    def get_children(cls):
        names = cls.children[cls.__name__.lower()]
        return [cls.registry[name] for name in names]
