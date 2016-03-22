
from functools import wraps
from inspect import signature


class _ObjectProxy(object):

    '''Basic (very) object proxy

    This class allows the universal decorators to behave like the
    object they are decorating
    '''

    def __init__(self, f):
        self.f = f
        try:
            self.__name__ = f.__name__
            self.__doc__ = f.__doc__
        except AttributeError:
            pass

    @property
    def __class__(self):
        return self.f.__class__

    def __getattr__(self, name):
        return getattr(self.f, name)


class _UniversalDecoratorWrapper_Methods(_ObjectProxy):

    '''Universal decorator boiler plate (methods only)

    This decorator class sets up a decorator for a given method.
    '''

    def __init__(self, f, instance, decorator):
        super(_UniversalDecoratorWrapper_Methods, self).__init__(f)
        self.decorator = decorator
        self.instance = instance

    def __call__(self, *args, **kwargs):
        # allows for decorators that don't care about instance
        sig = signature(self.decorator)
        if len(sig.parameters) == 3:
            return self._call_no_instance(*args, **kwargs)

        if self.instance is None:
            instance, args = args[0], args[1:]
            wrapped = functools.partial(self.f, instance)
            return self.decorator(wrapped, instance, args, kwargs)
        return self.decorator(self.f, self.instance, args, kwargs)

    def _call_no_instance(self, *args, **kwargs):
        return self.decorator(self.f, args, kwargs)


class _UniversalDecoratorWrapper_Functions(_ObjectProxy):

    '''Universal decorator boiler plate (functions, pass-thru for methods)

    This decorator class sets up a decorator for a given function.
    If the function is bound (i.e., a method), the descriptor is handled
    via __get__, so we pass the decorator through.
    '''

    def __init__(self, f, decorator):
        super(_UniversalDecoratorWrapper_Functions, self).__init__(f)
        self.decorator = decorator

    def __get__(self, instance, owner):
        '''Decorated bound method! Pass through!'''
        f = self.f.__get__(instance, owner)
        return _UniversalDecoratorWrapper_Methods(f, instance, self.decorator)

    def __call__(self, *args, **kwargs):
        try:  # handle decorators that don't care about instance
            return self.decorator(self.f, None, args, kwargs)
        except TypeError:
            return self.decorator(self.f, args, kwargs)


def universal_function_decorator(decorator):
    '''Delcare a function or method to be a decorator

    This decorator is intended to simplify creation of new decorators.
    Use for:
        * any function
        * most methods (not @classmethod or @staticmethod)
    Features:
        * Keeps information about original function
            * __name__ and __doc__ are updated
            * `inspect` values are updated
            * tracks instance (if given)
        * Allows for keyword arguments in
    Known problems:
        * not usable with @classmethod or @staticmethod

    Example:
        @universal_function_decorator
        def my_decorator(f, args, kwargs):
            # NOTICE the missing asterisks in the decorator declaration!!
            return f(*args, **kwargs)

        @my_decorator
        def my_function(*args, **kwargs):
            pass

    Example:
        @universal_function_decorator
        def my_decorator(f, instance, args, kwargs):
            # NOTICE the missing asterisks in the decorator declaration!!
            return f(*args, **kwargs)

        class MyClass(object):
            @my_decorator(mode='w')
            def my_function(*args, **kwargs):
                pass
    '''

    @wraps(decorator)
    def _ufd(f):
        return _UniversalDecoratorWrapper_Functions(f, decorator)
    return _ufd
