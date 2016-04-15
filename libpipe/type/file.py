
import os.path

from libpipe.type import base

import logging
log = logging.getLogger(__name__)


def factory(name=None, extns=[], parent=None):

    if not extns:
        msg = 'No extensions given. Try factory(extns=[".fa", ...])'
        raise ValueError(msg)
    else:
        # store as lowercase!
        extns = [extn.lower() for extn in extns]

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
            mcl.registry[class_id]
        except AttributeError:
            mcl.registry = {}  # init main registry
        except KeyError:
            pass  # class is not registered. Good!
        else:
            msg = 'Duplicate FileType found: {}'.format(class_id)
            raise ValueError(msg)

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

        if value:
            # Allowing an empty value matches other types, e.g.,
            # -- int() --> 0
            # -- str() --> ''
            try:
                self._check_extns(value)
            except AttributeError:
                # no extensions set--try children
                self._check_child_extns(value)
                # msg = 'Cannot call {} directly--no extension set'.format(
                #     self.__class__.__name__)
                # raise TypeError(msg)

    @classmethod
    def _check_extns(cls, value):
        '''Check given value for the presence of an expected extension.

        Performs a case insensitive check for the presence of an
        expected extension at the end of the given value.

        Instead of using os.path.splitext, use str.endswith. This allows
        for simpler checking given more complex extension types,
        e.g., *.rev.bt2. The downside is looping through the extensions
        instead of checking if the extension is in the set.

        Arguments:
            value: The value to check.
        Returns:
            None.
        Raises:
            ValueError if a valid extension not found.
        '''

        value = value.lower()  # case insensitive!
        for extn in cls.extns:
            if value.endswith(extn):
                return  # found a valid extension
        msg = 'Invalid extension (file: {}) for type {} [{}]'.format(
            os.path.basename(value), cls.__name__, ', '.join(cls.extns))
        raise ValueError(msg)

    @classmethod
    def _check_child_extns(cls, value):
        checked_extns = []
        for child in cls.get_children():
            try:
                child._check_extns(value)
            except ValueError:
                checked_extns.extend(child.extns)
            except AttributeError:
                pass  # child doesn't have extns either
            else:
                return  # we have a match!

        if not checked_extns:
            # No valid extensions at all!
            msg = 'Cannot call {} directly--no extensions set'.format(
                cls.__name__)
            raise TypeError(msg)
        else:
            msg = 'Invalid extension (file: {}) for type {} [{}]'.format(
                os.path.basename(value),
                cls.__name__, ', '.join(checked_extns)
            )
            raise ValueError(msg)

    @classmethod
    def get_children(cls):
        names = cls.children[cls.__name__.lower()]
        return [cls.registry[name] for name in names]

# __END__
