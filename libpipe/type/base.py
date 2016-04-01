

class TypeBase(object):

    @classmethod
    def factory(cls, name=None):
        name = name or 'SubType'
        _subcls = type(name, (TypeBase, ), dict())
        return _subcls
