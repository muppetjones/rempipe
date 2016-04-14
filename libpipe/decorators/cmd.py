
from libpipe.decorators.universal import universal_function_decorator


def universal_fall_through(func):
    '''Causes cmd.output to include cmd.input

    NOTE: This method only seems to work using the @decorator syntax.

    Decorate the output function of a CmdBase object to cause
    any input to be included in the output.

    Example:
        class SomeCmd(CmdBase):

            @universal_fall_through
            def output(self):
                return ['list', ]
    '''

    @universal_function_decorator
    def _fall_through(func, instance, args, kwargs):
        return instance.input() + func()
    return _fall_through(func)


def fall_through(func):
    '''Causes cmd.output to include cmd.input

    NOTE: This simpler version works (only) without the @decorator syntax.

    Wrap the output function of a CmdBase object to cause
    any input to be included in the output.

    Example:
        self.output = fall_through(self.output)
    '''
    def _fall_through():
        try:
            return func.__self__.input() + func()
        except AttributeError as e:
            msg = 'Fall through decorator: Try linking the object'
            raise AttributeError(msg) from e
    return _fall_through
