
import os.path

from libpipe.type import base
from libpipe.util import path


import logging
log = logging.getLogger(__name__)


def factory(name='IndexSubType', extns=[], counts=[], parent=None):
    '''Create an IndexType class

    Create and return an IndexMeta type, subclassed to str and parent.
    By subclassing to string, instances of the new class may be treated
    as path strings--adding new functionality (check for index files)
    without require a complete rewrite or specialization of how path
    strings are handled.

    Arguments:
        name: A string indicating the name that the new class should have.
        extns: A list of string extensions expected, e.g., ['.txt'].
        counts: The number of each extension expected. Defaults to 1 each.
    Returns:
        The new subclass.
    Raises:
        ValueError if no extensions given or if the number of counts
        does not match the number of given extns.
        TypeError if the given parent is not of type IndexMeta.

    Examples:
        # For a genome with a hisat2 index named GRCh38p5, i.e., GRCh38p5.*.ht2
        from libpipe.type import index

        # create the new index class
        name = 'Hisat2Index'
        extensions = ['.ht2']
        number_of_files_per_extension = [8]
        index_class = index.factory(
            name=name, extns=extensions, counts=number_of_files_per_extension)

        # create an instance of the class
        index_obj = index_class("~/genomes/human/GRCh38p5")
        print(index_obj)
        print(index_obj.path)
        print(index_obj.name)

        > ~/genomes/human/GRCh38p5
        > ~/genomes/human
        > GRCh38p5
    '''

    # Verify lineage
    if not parent:
        parent = IndexType
    if not type(parent) == IndexMeta:
        msg = 'Parent must be an IndexType'
        raise TypeError(msg)

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

    '''Define a metaclass for handling file indexes

    (Most?) sequence aligners require users to create a multi-file index
    of the reference genome. For the TopHat suite, these index files
    generally have the name <prefix>.N.<extn> and are identified (named)
    by their prefix, making simple identification via extension difficult.

    IndexMeta defines a new type that checks for the presence of the
    index files when given a prefix. IndexType is the base type. Each
    new IndexType is stored in a registry dict (by name), along with
    a dict lookup by lineage (__init__). This allows the base IndexType
    to check if a given prefix matches with a known index and, eventually,
    will allow for automatic detection of possible index types.

    IndexMeta classes will check for the presence of these index files
    under two conditions:
        * isinstance--iff a string is given.
        * object instantiation

    Example:
        # You have hisat2 index for the human genome named GRCh38p5 and
        # stored at the directory path ~/genomes/human.
        # The Hisat2Index is defined in libpipe.cmd.align with Hisat2Cmd,
        # and expects 8 files with the '.ht2' extension

        from libpipe.cmd import align
        index_prefix = '~/genomes/human/GRCh38p5'

        # check a string path prefix for index files
        isinstance(index_prefix, align.Hisat2Index)  # True!

        # create a index object
        hisat2index_obj = align.Hisat2Index(index_prefix)
    '''

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

        First, naive subclass check.
        IFF naive check failed AND the given instance is a string,
        then check for presence of expected index files (self or child).

        That last bit is crucial--IndexType does NOT define extensions,
        which means you can declare other types that don't set extensions
        and rely on the children checking.

        Arguments:
            inst: The object to check.
        Returns:
            True if inst is an instantiated IndexMeta type object or
            if a valid index was found that matches self OR any of
            self's children.
        '''

        # naive type check
        _check = any(
            self.__subclasscheck__(c)
            for c in {type(inst), inst.__class__}
            if c != str
        )

        # type(inst) should return inst.__class__,
        # so not str unless actually a str
        if not _check and type(inst) == str:
            try:
                self._check_extns(inst)
                _check = True  # if we can instantiate, then True!
            except ValueError:
                pass  # already False--just leave it

        return _check


class IndexType(metaclass=IndexMeta):

    '''A new type for recognizing indices stored across several files.

    Stores and validates a given path as an alignment index.
    New index types should be created using the factory--not
    by subclassing.

    Attributes:
        extns*: A list of expected extensions for the index.
        counts*: A list of expected file counts for each extension.
        path: The directory where the index is found.
        name: The basename of the index.

        * extn and counts are set via the factory.
    '''

    def __init__(self, value=None):

        if self.__class__ == IndexType:
            # No direct calls to base class!
            # -- this is clunky (code smell)
            # -- TODO(sjbush): Function that decides WHICH index is needed.
            msg = 'Cannot create IndexType object directly'
            raise TypeError(msg)

        if value:
            self.path, self.name = os.path.split(value)
            self._check_extns(value)

    def __eq__(self, other):
        '''Test for equality with an IndexType object

        Arguments:
            other: The object to compare against.
        Returns:
            True if str values equal, False otherwise.
            For IndexTypes, compares extns and counts, too.
        '''

        try:
            return self.extns \
                and sorted(self.extns) == sorted(other.extns) \
                and super().__eq__(other)
        except AttributeError:
            return super().__eq__(other)

    def __ne__(self, other):
        return not self.__eq__(other)

    @classmethod
    def get_children(cls):
        return [
            cls.registry[name.lower()]
            for name in cls.children[cls.__name__.lower()]
        ]

    @classmethod
    def _check_extns(cls, prefix):
        '''Checks the given path for existing files with expected extensions

        Check that the expected index files actually exist. If they don't,
        we complain. This is more than just a check that the path has
        the right extension (it shouldn't have one)--we check for files
        the file system.

        Arguments:
            prefix: The path prefix. E.g., 'dir/gen' for 'dir/gen.1.bt2'.
        Returns:
            None.
        Raises:
            ValueError if files with expected extension not found or
            if extensions not set.
        '''

        # NOTE: Perhaps this method should raise by default
        errs = []
        _path, name = os.path.split(prefix)

        # Find files with given extension and pattern
        # -- but, extns may not be set.
        try:
            files = path.walk_file(
                _path, extension=cls.extns, pattern=[name])

        except AttributeError:
            # No extensions, but maybe a child can match!
            msg = 'No extns defined for Index {}'.format(
                cls.__name__)
            errs = [msg]
            # raise ValueError(msg)  # Not yet--check the children first!

        else:
            # We have files with expected extensions
            # -- ensure we have as many as expected
            errs = []
            for extn, cnt in zip(cls.extns, cls.counts):
                matched = [f for f in files if f.endswith(extn)]
                if len(matched) != cnt:
                    errs.append('{}: {} v {}'.format(extn, len(matched), cnt))

        finally:
            if errs:
                # Something went wrong--no match for current object
                # -- try our children (return at first success)
                for child in cls.get_children():
                    try:
                        child._check_extns(prefix)  # recurse for each child
                    except ValueError:
                        pass  # one child failed...no big deal...keep checking
                    else:
                        return  # Found a match! (No ValueError raised)

                msg = 'Bad index: {}'.format('; '.join(errs))
                raise ValueError(msg)

        return
