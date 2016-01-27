
class RempipeError(Exception):

    '''Exception for throwing dynamic, pre-defined error messages

    NOTE: These exceptions are only meant for internal use. Children should
            also inherit from a specific exception type.
    NOTE: Children MUST inherit from this class BEFORE the exception type
            or the custom message will NOT work.

    Attributes:
        code    (str)   The string identifier for the message.
        insert  (list)  A list of values to insert into the error message.

    '''
    ERRMSG = {}

    def __init__(self, code, *args, details=[], ** kwargs):
        try:
            msg = self.ERRMSG[code].format(*details)
        except IndexError:
            msg = self.ERRMSG[code]  # missing expected values
        except KeyError:
            msg = code  # unknown error code
        finally:
            super().__init__(msg, *args, **kwargs)
