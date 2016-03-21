
import abc
import time


import logging
log = logging.getLogger(__name__)


class CmdInterface(metaclass=abc.ABCMeta):
    '''Define an interface for command-like objects

    Abstract methods:
        cmd: Return the bash-executable string of the command.
        link: Pair current 'output' with given command's 'input'.
            Should return the given command to allow chaining.
        output: Return a list of command output. Should be string paths
            in *most* cases.
    '''
    @abc.abstractmethod
    def cmd(self):
        pass

    @abc.abstractmethod
    def output(self):
        pass

    @abc.abstractmethod
    def link(self):
        pass


class CmdBase(CmdInterface):

    '''An abstract class for creating command-line objects

    Create a command line string for executing a given program. Children
    must define certain class attributes that characterize the program
    requirements. Certain methods may be overridden to allow specific
    handling and preparation of given arguments.

    Attributes:
        args: A list of positional arguments.
        flags: A list of flag-only arguments.
        kwargs: A dict of flag-value pairs.

        timestamp: A timestamp string. Default: init time. Useful for
            creating a unique time id that is syncronized within a pipe.

        strict: A bool indicating that only known arguments should be used.
        complain: A bool indicating whether or not strict violations should
            raise an exception (usually debugging only).

    Class methods for override (in call order)
            NOTE: The 'cmd' method calls several functions that may be
                overriden to customize how requirements are checked.
                None of these methods are implemented in the base class
        _pre_req(): Use to setup any last minute details before
            requirements are checked, such as handling the
            results of linking.
        _post_req(): Use to define additional requirements indirectly
            related to successful command execution. This should raise
            an exception if conditions are not met. Use sparingly.
        _pre_cmd(): Final modification before creating the executable string,
            e.g., such as setting up redirect.



    Example:
        # Flags may be given as positional arguments
        # -- i.e., we assume that ANY args starting with '-' are flags
        args = ['file.txt', '--verbose']
        kwargs = {'-o': 'output.txt'}

        # or as a list to kwargs
        # -- YAGNI?
        args = ['file.txt']
        kwargs = {
            'flags': ['--verbose', ],
            '-o': 'output.txt',
        }

        cmd = Cmd(*args, **kwargs)
    '''

    def __init__(
            self, *args,
            complain=False, flags=[], strict=True, timestamp=None,
            **kwargs
    ):

        # rename so we can modify with impunity
        # (and avoid unnecessary if statements)
        _args = list(args)
        _flags = None  # this will get changed later
        _kwargs = dict.copy(kwargs)

        # save a timestamp, given or generated
        self.timestamp = (timestamp if timestamp
                          else time.strftime("%y%m%d-%H%M%S"))

        # expected flags and kwargs
        # -- should we save this to avoid computing everytime
        #    or keep to allow it to change? (not likely, probably shouldn't)
        expected_flags = [
            flg[0] for flg in self.attr.args
            if flg[0] and not isinstance(flg[0], int)
        ]

        # separate flags from args
        _flags = [arg for arg in _args if arg.startswith('-')]
        _args = [arg for arg in _args if arg not in _flags]
        _flags.extend(flags)

        # remove unknown flags and kwargs if strict (default)
        if strict:
            if complain:
                _unk_flags = [v for v in _flags if v not in expected_flags]
                _unk_kw = [k for k in kwargs.keys() if k not in expected_flags]
                if _unk_flags or _unk_kw:
                    msg = 'Unknown arguments: {}'.format(
                        ', '.join(_unk_kw + _unk_flags))
                    raise ValueError(msg)
            _flags = [v for v in _flags if v in expected_flags]
            _kwargs = {k: v for k, v in kwargs.items() if k in expected_flags}

        self.args = _args  # deep copy list--but not elements
        self.flags = _flags
        self.kwargs = dict.copy(self.attr.defaults)
        self.kwargs.update(_kwargs)

        self.strict = strict
        self.complain = complain

    #
    #   Public access
    #

    def cmd(self):
        '''Generate the command string

        Build the command string using the given arguments/inputs.
        Performs a number of processing steps before actually creating the
        command string:
            1. Match up the given input with arguments (uses attr.req_type).
            2. Process arguments as needed BEFORE checking requirements.
                > via '_pre_req', defined per child class.
                > e.g., ensure correct argument order.
            3. Check that requirements are met (via attr.req_[kw]args).
            4. Check for any additional requirements.
                > via '_post_req', defined by child.
                > e.g., ensure a complete Bowtie index exists.
            5. Process arguments as needed AFTER checkeing requirements.
                > via '_pre_cmd', as defined per child class.
                > e.g., defining a redirect based on args.
            6. Generate the command string.
            7. Wrap the command string (if given).

        Arguments:
            *args, **kwargs: Catch (and ignore) anything passed in.
            wrap: See main docstring. Will wrap command in a BASH-style
                assignment to a variable with the given name.
        Return:
            A BASH executable string.
        Raises:
            TypeError if link input is bad/unusable or if wrong file type.
            AttributeError if bad kwargs.
            IndexError if bad args.
        '''

        self._pre_req()
        self._post_req()
        self._check_requirements()
        self._pre_cmd()

    def help(self):
        '''Return usage text'''
        raise NotImplementedError()

    def link(self):
        pass

    #
    #   To override
    #   NOTE: Don't forget about "output"
    #

    def _pre_req(self):
        pass  # must be defined by child (but don't care if not set)

    def _post_req(self):
        pass  # must be defined by child (but don't care if not set)

    def _pre_cmd(self):
        pass  # must be defined by child (but don't care if not set)

    #
    #   "Protected" methods
    #   NOTE: Don't override unless ABSOLUTELY NECESSARY!
    #

    def _check_requirements(self):
        '''Ensure all defined requirements are met.

        For each given argument, ensure that all requirements defined
        within our CmdAttributes are met.
        '''

        try:
            self.__check_args()
        except IndexError as e:
            log.error(str(e))
            raise

        try:
            self.__check_kwargs()
        except KeyError as e:
            log.error(str(e))
            raise

        try:
            self.__check_type()
        except ValueError as e:
            log.error(str(e))
            raise

    #
    #   Private methods
    #   NOTE: These are unavailable in the children.
    #

    def __check_args(self):
        '''Ensure non-type, positional argument requirements are met.

        Only one requirement for now: Do we have the expected number of
        positional arguments?

        NOTE: Does not check for a maximum number of positional args.

        Raises:
            IndexError on requirement failure.
        '''

        n_missing = self.attr.req_args - len(self.args)
        if n_missing > 0:
            msg = 'Missing positional arguments: {} of {}'.format(
                n_missing, self.attr.req_args)
            raise IndexError(msg)

    def __check_kwargs(self):
        '''Ensure non-type, keyword arg  requirements are met.

        Check that the given keyword args are present, if required, and
        meet the expected relational requirements, i.e., AND, XOR.

        Raises:
            KeyError on failure
        '''

        simple = [kw for kw in self.attr.req_kwargs if isinstance(kw, str)]
        missing = [kw for kw in simple if kw not in self.kwargs]

        # compound requirement -- require all or none of args in the tuple
        compound = [kw for kw in self.attr.req_kwargs if kw not in simple]
        log.debug(compound)
        for cmpd in compound:

            # tuple indicates all or nothing
            if isinstance(cmpd, tuple):
                missed = [kw for kw in cmpd if kw not in self.kwargs]
                if len(missed) != len(cmpd):  # remember, nothing is fine!
                    missing.append('AND({})'.format(', '.join(missed)))

            # a list indicates exactly one (must have at least one)
            elif isinstance(cmpd, list):
                found = [kw for kw in cmpd if kw in self.kwargs]
                if len(found) != 1:
                    missing.append('XOR({})'.format(', '.join(found)))

        if missing:
            msg = 'Missing required keyword arg: {}'.format(', '.join(missing))
            raise KeyError(msg)

    def __check_type(self):
        '''Ensure that pos and kw args match the expected type

        Raises:
            ValueError on failure
        '''

        pass
