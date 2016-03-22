
from functools import partial

from libpipe.decorators.universal import universal_function_decorator
from libpipe.util import path

import logging
log = logging.getLogger(__name__)


def file_or_handle(func=None, *, mode='r'):
    '''Decorator for accepting a file or handle

    Ensures a file-like object is passed to the wrapped function.

    NOTE: Sets `filename` attribute on instance, if exists.

    Arguments:
        mode (keyword, default='r'): See FileIO for options.

    Examples:
        @file_or_handle
        def read_some_data(fh):
            return fh.read()

        # 'mode' must be given as a keyword
        @file_or_handle(mode='w')
        def write_some_data(fh, data):
            fh.write(data)
    '''

    # allow decorator to be called without arguments
    if func is None:
        return partial(file_or_handle, mode=mode)

    @universal_function_decorator
    def _foh(func, instance, args, kwargs):
        fh, args = args[0], args[1:]
        if hasattr(fh, 'write') or fh is None:
            return func(fh, *args, **kwargs)

        else:
            protected_path = path.protect(fh)
            with open(protected_path, mode) as handle:
                try:
                    instance.filename = protected_path
                except AttributeError:
                    pass
                return func(handle, *args, **kwargs)

    return _foh(func)
