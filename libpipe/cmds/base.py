import os.path
import time
from abc import ABCMeta

from libpipe.cmds.help import HelpCmd

import logging
log = logging.getLogger(__name__)


class BaseCmd(metaclass=ABCMeta):

    '''A base (abstract) class for handling command line calls.

    Defines a number of methods and attributes for creating command
    line, well, commands.

    NOTE: Argument requirements are checked in `cmd`.

    Class Attributes:
        NAME        A human readable string for identifying the command.
        INVOKE_STR  The string used to invoke the command.
        ARGUMENTS   A list of arguments defined by (<flag>, <type>, <descr>).
                    When <flag> is None, the argument is positional, and
                    when <type> is None, the argument is just a flag.
        DEFAULTS    Default values for keyword arguments. Make sure that
                    the flags have expected leading hyphens, e.g., '-f'.
        HELP_DICT   A dict containing elements for creating a HelpCmd object.
                    More of a reference for the class than anything.
                    Could be used to type check eventually.
        HELP        A HelpCmd object (initialized during help call)

        REQ_KWARGS  A list of the required keyword args, e.g., ['-f'].
                    A tuple indicates required together, e.g., [('-1', '-2')].
                    A list indicates XOR requirement, e.g., [['-U', '-1']].
        REQ_ARGS    Int denoting how many positional arguments are expected.
        REQ_TYPE    A list of kwarg flags and their expected extensions.
                    E.g., [ [('-f', ), ('.txt', )], ... ]

    Class Methods:
        __prepreq__()   Use to setup any last minute details before
                        requirements are checked, such as handling the
                        results of linking.
        __prepcmd__()   Use to setup any last minute details before cmd
                        is called, such as setting up redirect

    SEE ALSO:
        HelpCmd
    '''

    # command defining attributes
    NAME = None  # human readable name of command, e.g., 'samtools_sort'
    INVOKE_STR = None  # string to invoke command, e.g., 'samtools sort'
    ARGUMENTS = None  # see HELP_DICT['arguments'] for example
    DEFAULTS = {}  # starting values for kwargs

    # Help related attributes
    HELP_DICT = {  # help text
        'synopsis': 'cmd_invoke -f FILE [-n INT, -h] FILE ...',
        'description': 'A description of the command.',
        'arguments': [
            (None, 'FILE', 'Example input one'),
            (None, 'FILE', 'Example input two'),
            ('-f|--file', 'FILE', 'Example keyword argument'),
            ('-n', 'INT', 'Example keyword argument'),
            ('-v', None, 'Example flag argument'),
        ],
    }
    HELP = None  # instantiated at first call to 'help'

    # set argument requirements
    # NOTE: required kwargs is a list of the required argument flags
    # NOTE: required args is an int of the number of required args
    # NOTE: flags, by definition, are optional
    REQ_KWARGS = []
    REQ_ARGS = 0
    REQ_TYPE = []

    #
    #   Custom Exceptions
    #

    class CmdLinkError(ValueError):
        pass

    #
    #   Magic methods
    #

    def __init__(self, *args, **kwargs):
        '''Initialize command object

        Create a command object using given paramters.
        Flags, e.g., '-v', should be given as positional parameters.
        '''

        # save a timestamp, if passed
        try:
            self.timestamp = kwargs['timestamp']
            del kwargs['timestamp']
        except KeyError:
            self.timestamp = time.strftime("%y%m%d-%H%M%S")
            # self.timestamp = None

        # ensure the expected kwarg hyphens are in place
        # NOTE: Gives flags with len > 1 a double hyphen, e.g., '--',
        #       SO if it only takes a single hyphen MAKE SURE YOU PUT IT!
        # NOTE: we need to do this BEFORE checking requirements
        #       to ensure the given kwargs match the expected kwargs
        def ensure_hyphen(flag):
            if not flag.startswith('-'):
                fmt = '-{}' if len(flag) == 1 else '--{}'
                flag = fmt.format(flag)
            return flag
        kwargs = {ensure_hyphen(k): v for k, v in kwargs.items()}

        # separate flags from args
        # -- best option would be to set defaults in a custom init
        try:
            flags = [v for v in args if v.startswith('-')]
        except AttributeError:
            msg = 'Make sure to expand keyword args, e.g., Cmd(**kwargs)'
            raise AttributeError(msg)
        args = [v for v in args if v not in flags]

        # make sure we deep copy defaults and args
        self.redirect = None  # self.REDIRECT

        self.kwargs = {}
        self.kwargs.update(self.DEFAULTS)
        self.kwargs.update(kwargs)

        self.args = []
        self.args.extend(args)

        self.flags = []
        self.flags.extend(flags)

    def __str__(self):
        '''Returns the BASH executable string'''
        return self.cmd()

    #
    #   Property methods
    #

    @property
    def name(self):
        '''Name property. Do NOT override'''
        return self.NAME

    @property
    def invoke_str(self):
        '''Invocation string property. Do NOT override (unless necessary)'''
        return self.INVOKE_STR

    #
    #   Class methods
    #

    @classmethod
    def _filter_by_type(self, args, filt):
        filtered = [
            arg for arg in args
            if os.path.splitext(arg)[1] in filt
        ]
        return filtered

    @classmethod
    def _trubase(self, path_name):
        return os.path.splitext(os.path.basename(path_name))[0]

    #
    #   "Public" access
    #

    def output(self):
        '''Return list of created files.

        Must override. Should at the very least return the input to allow
        easier chaining.
        '''
        return None

    def link(self, cmd):
        '''Set dest input to self output

        Arguments:
        cmd     Another command object.
        Return:
        The given command object (for chaining)
        '''

        cmd.input = self.output
        return cmd

    def cmd(self, *args, **kwargs):
        '''Run command preprocessing and return command'''

        # run requirements prep, if provided by child
        # -- use for ensuring the requirements are met
        try:
            self.__prepreq__()
        except AttributeError:
            pass  # may not be set on child class

        # ensure we can create the command as expected
        try:
            self._check_requirements()
        except ValueError:
            log.error(self.help())
            raise

        # run command prep, if provided by child
        # -- use to ensure all data is ready for the command
        # -- NOTE: should require the necessary arguments
        try:
            self.__prepcmd__()
        except AttributeError:
            pass  # may not be set on child class

        return self.__cmd__(*args, **kwargs)

    def __prepreq__(self):
        try:
            self.__input__()  # should call self._get_input
        except AttributeError:
            return  # nothing to do

    def _get_input(self):
        '''Matches linked input to expected arguments. Must override.'''

        try:
            args = self.input()
        except AttributeError:
            return None
        except TypeError:
            raise TypeError('Bad link: "input" is not callable')

        return args

    def __cmd__(self, readable=True):
        '''Create BASH executable string.

        Arguments:
            readable (bool) Denotes whether command should be human readable
        Returns:
            A BASH executable string.
        Raises:
            ValueError if argument requirements are not met.
        '''

        # put the command on separate lines for human readable version
        sep = ' \\\n  ' if readable else ' '

        # make strings from given parameters
        flags = ' '.join(self.flags)  # flags don't need to be separated
        kwargs = sep.join(
            "{} {}".format(k, v)
            for k, v in sorted(self.kwargs.items())
        )
        args = sep.join(self.args)
        cmd_parts = filter(  # remove missing elements
            None, [self.invoke_str, flags, kwargs, args, self.redirect])
        return sep.join(cmd_parts)

    def help(self):
        '''Return usage text

        NOTE: Creates the HelpCmd object iff not found
        '''

        if self.HELP is None:
            try:
                self.HELP_DICT['arguments'] = self.ARGUMENTS
            except AttributeError:
                pass  # alias
            self.HELP_DICT['name'] = self.NAME
            self.HELP = HelpCmd(**self.HELP_DICT)
        return str(self.HELP)

    #
    #   "Private" access
    #

    def _check_requirements(self):
        '''Ensure all argument requirements are fulfilled

        Returns:
            True if all requirements fulfilled.
        Raises:
            ValueError with missing requirements if requirement not met
        '''

        # check for required kwargs
        try:
            missing = self._missing_kwargs()
        except ValueError:
            raise  # > 1 XOR option given
        else:
            if missing:
                raise ValueError('Missing arguments:\n\t{}\n'.format(
                    '\n\t'.join(str(m) for m in missing)
                ))

        # check for expected number of args
        if len(self.args) < self.REQ_ARGS:
            raise ValueError(
                'Missing {} of {} positional parameters; '.format(
                    self.REQ_ARGS - len(self.args), self.REQ_KWARGS
                ))

        # check for expected file types
        try:
            all_valid = self._check_type()
        except ValueError:
            raise  # Bad file type

        return

    def _check_type(self):
        '''Check specified kwargs to ensure given expected file type.

        NOTE: All extensions in REQ_TYPE should have leading periods.
        TODO: Update to use Django model to handle file types.
        '''

        for req in self.REQ_TYPE:
            flags, types = req
            # flags = [f for f in flags if f in self.kwargs]
            for f in flags:
                try:
                    extn = os.path.splitext(self.kwargs[f])[1]

                except KeyError:
                    try:
                        extn = os.path.splitext(self.args[f])[1]
                    except (TypeError, IndexError):
                        continue  # we don't have the flag, so skip check

                if extn not in types:
                    raise ValueError(
                        'Invalid extension "{}" for arg {}'.format(extn, f))
        return True

    def _missing_kwargs(self):
        '''Make sure that the given kwargs contain all required kwargs

        Checks all elements in REQ_KWARGS.
        > Strings are simple checks against the kwargs attribute.
        > Lists are XOR required keywords, e.g., one and only one req.
        > Tuples are AND required keywords, e.g., if one, then all.

        Example REQ_KWARGS:
            ['-f', ('-1', '-2'), ['-1', '-U']]

        Arguments:
            A kwargs dict.
        Returns:
            A list of missing kwargs.
        '''

        # simple requirement
        simple = [kw for kw in self.REQ_KWARGS if isinstance(kw, str)]
        missing = [kw for kw in simple if kw not in self.kwargs]

        # log.debug('in _missing_kwargs: {} vs {}'.format(
        #     self.REQ_KWARGS, self.kwargs))

        # compound requirement -- require all or none of args in the tuple
        compound = [kw for kw in self.REQ_KWARGS if kw not in simple]
        for cmpd in compound:

            # tuple indicates all or nothing
            if isinstance(cmpd, tuple):
                # log.debug('{}: AND'.format(cmpd))
                missed = [kw for kw in cmpd if kw not in self.kwargs]
                missing.extend(missed)

            # a list indicates one and only one
            elif isinstance(cmpd, list):
                # log.debug('{}: XOR'.format(cmpd))
                found = [kw for kw in cmpd if kw in self.kwargs]
                if not found:
                    missing.extend(cmpd)
                elif len(found) != 1:
                    raise ValueError(
                        'More than one arg given: {}'.format(cmpd))

        return missing
