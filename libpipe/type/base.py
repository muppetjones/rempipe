

class TypeBase(object):

    @classmethod
    def factory(cls):
        class SubType(TypeBase):
            pass
        return SubType
