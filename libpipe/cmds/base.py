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
                    A single requirement will have two to three elements:
                        <flags tuple>, <extension tuple>, [match any bool]
                    For example, the following type requirement would
                    naively define the input required by hisat:
                        REQ_TYPE = [
                            [('-1', '-2'), ('.fq', '.fa'), False],
                            [('-U', ),     ('.fq', '.fa'), False],
                            [('-S', ),     ('sam', )],
                        ]
                    The 'match_any' bool element is used only for matching
                    linked command input and is assumed True if ommitted.
                    If False, the option will not be matched unless the
                    given flags are matched exactly. In the above example,
                    a single linked fastq file would be set to '-U' and three
                    linked fastq files would not be used at all.
                    See '__match_input' and '_check_type' for more details.

    Public Methods:
        input()     Set after linking to another command. References
                    the linked command's output method.
        output()    Return a list of output files.
        link(CMD)   Link two commands: Set the input of CMD to this object's
                    output method.
                    Return the linked command to enable chaining.
        cmd()       Return the BASH-style command string.
                    All requirement checking is done here. If linked, the
                    input method is used to update defined arguments.
        help()      Returns a usage message. Uses class attributes NAME,
                    INVOKE_STR, ARGUMENTS, SYNOPSIS, and DESCRIPTION.


    Class Methods for override (in call order):
        _prepreq()  Use to setup any last minute details before
                    requirements are checked, such as handling the
                    results of linking.
                    NOTE: Call parent to ensure intended functionality.
        _additional_requirements():
                    Use to define additional requirements related
                    to successful command execution. This should raise
                    an exception if conditions are not met.
                    NOTE: This method is not implemented in BaseCmd.
        _prepcmd()  Use to setup any last minute details before cmd
                    is called, such as setting up redirect.
                    NOTE: This method is not implemented in BaseCmd.

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

    class CmdLinkError(TypeError):
        ERRMSG = {
            'input': 'Bad link: no input given or input not callable',
            'mismatch': 'Bad link: unexpected input type',
        }

        def __init__(self, key, *args, **kwargs):
            try:
                super().__init__(self.ERRMSG[key], *args, **kwargs)
            except KeyError:
                super().__init__(key, *args, **kwargs)

    class PositionalArgError(IndexError):
        pass

    POSITIONALARGERROR = {
        'missing': '{}: Missing {} of {} required positional parameters',
    }

    class KeywordArgError(AttributeError):
        ERRMSG = {
            'missing': '{}: Missing required keyword arguments: {}',
            'unknown': 'Unrecognized keyword argument given: "{}"',
        }

        def __init__(self, key, *args, val=[], ** kwargs):
            try:
                super().__init__(
                    self.ERRMSG[key].format(*val), *args, **kwargs)
            except KeyError:
                super().__init__(key, *args, **kwargs)

    class FileTypeError(TypeError):
        pass

    #
    #   Magic methods
    #

    def __init__(self, *args, strict=True, timestamp=None, ** kwargs):
        '''Initialize command object

        Create a command object using given paramters.
        Flags, e.g., '-v', should be given as positional parameters.
        '''

        # save a timestamp, given or generated
        self.timestamp = (timestamp if timestamp
                          else time.strftime("%y%m%d-%H%M%S"))

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

        if strict:
            # check for unknown kwargs
            # NOTE: We MUST have removed timestamp AND
            #       ensured hyphens BEFORE this
            try:
                known_flags = [v[0] for v in self.ARGUMENTS if v[0]]
            except TypeError:
                raise NotImplementedError(
                    'ARGUMENTS must be set on child class')
            unknown_flags = [
                k for k in flags + list(kwargs.keys())
                if k not in known_flags
            ]
            if unknown_flags:
                raise self.KeywordArgError('unknown', val=unknown_flags)

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

    def cmd(self, *args,  **kwargs):
        '''Run command preprocessing and return command'''

        # run requirements prep, if provided by child
        # -- use for ensuring the requirements are met
        self._prepreq()

        # ensure we can create the command as expected
        try:
            self.__check_requirements()
        except (
            self.KeywordArgError,
            self.PositionalArgError,
            self.FileTypeError
        ):
            log.error(self.help())
            raise

        # additional requirements defined by the child class
        # for the purpose of ensuring successful command execution
        try:
            self._additional_requirements()
        except AttributeError:
            pass  # ignore if it's not defined

        # run command prep, if provided by child
        # -- use to ensure all data is ready for the command
        # -- NOTE: should require the necessary arguments
        try:
            self._prepcmd()
        except AttributeError:
            pass  # may not be set on child class

        return self._cmd(*args, **kwargs)

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

    def _has_output(self):
        '''Checks for expected output from a command

        Arguments:
            cmd     Command object.
        Returns:
            True if ALL output is found.
            False otherwise
        '''

        for f in self.output():
            if not os.path.isfile(f):
                return False
        return True

    def _prepreq(self):
        try:
            self.__match_input()
        except AttributeError:
            return  # nothing to do

    def __match_input(self):
        '''Match linked inputs to keyword arguments according to REQ_TYPE

        Uses REQ_TYPE to match linked input to a keyword argument or to
        positional arguments IFF REQ_ARGS > 1, e.g.,
            input: ['a.txt', 'a.fq', 'b.fq']
            REQ_ARGS = 1
            REQ_TYPE = [[('-b', '-a'), ('.fq',)]]  <-- notice the order!
        will result in the following argument conditions:
            kwargs: {'-b': 'a.fq', '-a': 'b.fq'}
            args: ['a.txt']

        Return:
            None
        Raises:
            Attribute error if no link or input is found.
            TypeError if linked input is not callable.
        '''

        # Ensure we have valid, linked input
        try:
            args = self.input()
            n_args = len(args)
        except AttributeError:
            return None
            # raise AttributeError('Attribute "input" is not set (no link)')
        except TypeError:
            raise self.CmdLinkError('input')
        # else:
            # if n_args == 0:
            #     raise self.CmdLinkError('no_input')

            # match with required type
        for req_type in self.REQ_TYPE:
            try:
                flag_list, type_list, match_any = req_type
            except ValueError:
                flag_list, type_list = req_type
                match_any = True

            # find args with matching file type and update kwargs
            filtered = self._filter_by_type(args, type_list)
            if match_any or len(filtered) == len(flag_list):
                self.kwargs.update(dict(zip(flag_list, filtered)))
            else:
                continue  # we don't want to remove if we didn't use!!

            # remove used
            args = [a for a in args if a not in filtered]

        # if a positional argument was defined in REQ_TYPES,
        # it is now stored in kwargs -- move it to args and clean kwargs
        pos_kw = [k for k in self.kwargs.keys() if isinstance(k, int)]
        self.args.extend([self.kwargs[pos] for pos in sorted(pos_kw)])
        for pos in pos_kw:
            del self.kwargs[pos]  # delete instead of dict comp

        # otherwise, use as required args IFF expected
        if len(self.args) < self.REQ_ARGS:
            n_still_required = self.REQ_ARGS - len(self.args)
            self.args.extend(args[:n_still_required])

        # finally, if we haven't used ANY of the input, raise CmdLinkError
        else:
            if len(args) == n_args:
                raise self.CmdLinkError('mismatch')

        return

    def _cmd(self, verbose=True):
        '''Create BASH executable string.

        Arguments:
            readable (bool) Denotes whether command should be human readable
        Returns:
            A BASH executable string.
        Raises:
            ValueError if argument requirements are not met.
        '''

        # put the command on separate lines for human readable version
        sep = ' \\\n  ' if verbose else ' '

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

    def __check_requirements(self):
        '''Ensure all argument requirements are fulfilled

        Returns:
            True if all requirements fulfilled.
        Raises:
            KeywordArgError if invalid keyword(s) given.
            PositionalArgError if missing positional arguments.
            FileTypeError if given files do not match expected type.
        '''

        # check for required kwargs
        try:
            missing = self.__check_kwargs()
        except self.KeywordArgError:
            raise  # > 1 XOR option given
        else:
            if missing:
                raise self.KeywordArgError('missing', val=[
                    self.name, ', '.join(str(m) for m in missing)
                ])

        # check for expected number of args
        if len(self.args) < self.REQ_ARGS:
            msg = self.POSITIONALARGERROR['missing'].format(
                self.__class__.__name__,
                self.REQ_ARGS - len(self.args), self.REQ_ARGS,
            )
            raise self.PositionalArgError(msg)

        # check for expected file types
        try:
            self.__check_type()
        except self.FileTypeError:
            raise  # Bad file type

        return

    def __check_type(self):
        '''Check specified kwargs to ensure given expected file type.

        NOTE: All extensions in REQ_TYPE should have leading periods.
        TODO: Update to use Django model to handle file types.

        Returns:
            True if expected file types match given files.
        Raises:
            TypeError if a given argument has wrong file type.
        '''

        for req in self.REQ_TYPE:
            try:
                # match_any used for matching only, not checking
                flags, types, match_any = req
            except ValueError:
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
                    msg = 'Invalid file type "{}" given for argument {}'
                    raise self.FileTypeError(msg.format(extn, f))
        return True

    def __check_kwargs(self):
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

        # compound requirement -- require all or none of args in the tuple
        compound = [kw for kw in self.REQ_KWARGS if kw not in simple]
        for cmpd in compound:

            # tuple indicates all or nothing
            if isinstance(cmpd, tuple):
                missed = [kw for kw in cmpd if kw not in self.kwargs]
                # log.debug('{}: AND {}'.format(cmpd, missed))
                if len(missed) != len(cmpd):  # remember, nothing is fine!
                    missing.extend(missed)

            # a list indicates one and only one
            elif isinstance(cmpd, list):
                # log.debug('{}: XOR'.format(cmpd))
                found = [kw for kw in cmpd if kw in self.kwargs]
                if not found:
                    missing.extend(cmpd)
                elif len(found) != 1:
                    raise self.KeywordArgError(
                        'More than one arg given: {}'.format(cmpd))

        return missing
