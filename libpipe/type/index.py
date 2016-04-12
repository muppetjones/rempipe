
import os.path

from libpipe.type import base
from libpipe.util import path


import logging
log = logging.getLogger(__name__)


def factory(name='IndexSubType', extns=[], counts=[], parent=None):
    '''Create an IndexType class

    Create and return an IndexMeta type, subclassed to str and parent.

    Arguments:
        name: A string indicating the name that the new class should have.
        extns: A list of string extensions expected, e.g., ['.txt'].
        counts: The number of each extension expected. Defaults to 1 each.
    Returns:
        The new subclass.
    Raises:
        ValueError if no extensions given
        or if the number of counts does not match the number of given extns.
    '''

    # Verify lineage
    if not parent:
        parent = IndexType
    if not issubclass(parent, IndexType):
        msg = 'Parent must be an IndexType'
        raise ValueError(msg)

    if not extns:
        msg = 'No extensions given. Try factory(extns=[".fa", ...])'
        raise ValueError(msg)

    if not counts:
        counts = [1] * len(extns)  # assume only one of each extn required!

    # Check for existing type
    # -- return the type if found
    # -- raise ValueError if same name, different extensions
    if name.lower() in IndexType.registry:
        existing = IndexType.registry[name.lower()]
        if existing.extns == extns and existing.counts == counts:
            return existing

        msg = 'IndexType with name {} already exists'.format(name)
        raise ValueError(msg)

    if counts and len(extns) != len(counts):
        msg = 'Unequal num of extns ({}) vs counts ({})'.format(
            len(extns), len(counts))
        raise ValueError(msg)

    _subcls = type(name, (parent, str), dict(extns=extns, counts=counts))
    return _subcls


class IndexMeta(base.TypeMeta):

    def __init__(mcl, name, parents, dct):
        '''Add Registry and Children to new classes

        Adds a simple registry {name: cls} for all sub classes.
        Also creates a child lookup {name: [child name ...]}.
        '''

        try:
            class_id = name.lower()
            mcl.registry[class_id] = mcl
        except AttributeError:
            mcl.registry = {class_id: mcl}

        # New type might not be registered yet
        try:
            mcl.children[class_id] = []
        except AttributeError:
            mcl.children = {class_id: []}

        for parent in parents:
            if not type(parent) == IndexMeta:
                continue  # only add children lookup to IndexMeta

            # add new class to parent registry
            try:
                parent_id = parent.__name__.lower()
                mcl.children[parent_id].append(class_id)
            except KeyError:
                mcl.children[parent_id] = [class_id]

        super(IndexMeta, mcl).__init__(name, parents, dct)

    def __instancecheck__(self, inst):
        '''Customize return value for `isinstance`

        If the given instance can be initialized to the type, we
        return True. IndexType restrictively checks for presence
        of expected number of specific files with a given extension.
        '''

        _check = any(
            self.__subclasscheck__(c)
            for c in {type(inst), inst.__class__}
            if c != str
        )

        # type(inst) should return inst.__class__,
        # so not str unless actually a str
        if not _check and type(inst) == str:
            try:
                self(inst)
                _check = True  # if we can instantiate, then True!
            except ValueError:
                pass  # already False--just leave it

        return _check


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

    def __eq__(self, other):
        '''Test for equality with an IndexType object

        Arguments:
        other: The object to compare against.
        Returns:
        True if str values equal, False otherwise.
        For IndexTypes, compares extns and counts, too.
        '''

        try:
            return sorted(self.extns) == sorted(other.extns) \
                and super().__eq__(other)
        except AttributeError:
            return super().__eq__(other)

    def __ne__(self, other):
        return not self.__eq__(other)

    @classmethod
    def get_children(cls):
        return cls.children[cls.__name__.lower()]

    def _check_extns(self):
        files = path.walk_file(
            self.path, extension=self.extns, pattern=[self.name])

        if not self.extns:
            msg = 'No extensions set'
            raise ValueError(msg)

        errs = []
        for extn, cnt in zip(self.extns, self.counts):
            matched = [f for f in files if f.endswith(extn)]
            if len(matched) != cnt:
                errs.append('{}: {} v {}'.format(extn, len(matched), cnt))

        if errs:
            msg = 'Bad index: {}'.format('; '.join(errs))
            raise ValueError(msg)
