'''Define a series of file types

Hierarchy:
    TypeMeta > TypeBase > [subtype] > [subtype object]
    `factory` returns the subtype class, which is used
    to create subtype objects.

Examples:
    import libpipe.type.base as _type
    SeqType = _type.factory('SeqType', extn=['.fq', '.fa'])

    try:
        seqfile = SeqType('path/to/seq.txt')
    except ValueError:
        msg = 'Bad type'
        raise ValueError(msg)
'''

# NOTE: Each subtype is NOT a TypeBase, but a TypeMeta!

import logging
log = logging.getLogger(__name__)


def factory(name=None):
    name = name if name else 'SubType'
    _subcls = type(name, (TypeBase, ), dict())
    return _subcls


class TypeMeta(type):
    pass

    # def __new__(cls, name, parents, dct):
    #     # if 'class_id' not in dct:
    #     #     dct['class_id'] = name.lower()
    #
    #     return super(TypeMeta, cls).__new__(cls, name, parents, dct)
    #
    # def __init__(cls, name, bases, dct):
    #     if not hasattr(cls, 'registry'):
    #         cls.registry = {}  # this is base class--create empty registry
    #     else:
    #         # derived class -- add to registry
    #         type_id = name.lower()
    #         cls.registry[type_id] = cls
    #
    #     super(TypeMeta, cls).__init__(name, bases, dct)


class TypeBase(metaclass=TypeMeta):

    def __init__(self, _value):
        pass
