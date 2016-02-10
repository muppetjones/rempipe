import os.path
import time
import abc
import re

from remsci.lib.decorators import universal_function_decorator

from libpipe.cmds.help import HelpCmd
from libpipe.utility.exceptions import RempipeError

import logging
log = logging.getLogger(__name__)


class CmdInterface(metaclass=abc.ABCMeta):

    @abc.abstractmethod
    def output(self):
        pass

    @abc.abstractmethod
    def link(self, cmd):
        pass

    @abc.abstractmethod
    def cmd(self):
        pass

    # @abc.abstractmethod
    # def help(self):
    #     pass


def universal_fall_through_decorator(func):
    '''A fall through decorator using the UFD

    NOTE: This only seems to work using the @decorator syntax.
    '''

    @universal_function_decorator
    def _fall_through_decorator(func, instance, args, kwargs):
        if instance is None:
            instance = func.__self__
        return instance.input() + func()
    return _fall_through_decorator(func)


def fall_through_decorator(func):
    '''A fall through decorator

    NOTE: This simpler version works without the @decorator syntax.
    Example:
        self.func = fall_through_decorator(self.func)
    '''
    def _fall_through_decorator():
        return func.__self__.input() + func()
    return _fall_through_decorator


class CmdAttributes(object):

    '''Store attributes related to handling command line calls.

    Attributes
        name        A human readable string for identifying the command.
        synopsis    [optional] A brief synopsis of the command.
        description [optional] A description of the command.
        invoke_str  The string used to invoke the command.
        arguments   A list of arguments defined by (<flag>, <type>, <descr>).
                    When <flag> is None, the argument is positional, and
                    when <type> is None, the argument is just a flag.
                    Examples:
                        (None, 'FILE',  'The input file'),
                        ('-h', None,    'Display the help text'),
                        ('-o', 'FILE',  'The output file')
        defaults    Dict of default values for keyword arguments.
                    Flags have expected leading hyphens, e.g., '-f'.

        req_args    Int denoting how many positional arguments are expected.
        req_kwargs  A list of the required keyword args, e.g., ['-f'].
                    A tuple indicates required together, e.g., [('-1', '-2')].
                    A list indicates XOR requirement, e.g., [['-U', '-1']].
        req_type    A list of kwarg flags and their expected extensions.
                    A single requirement will have two to three elements:
                        <flags tuple>, <extension tuple>, [bool match only]
                    For example, the following type requirement would
                    naively define the input required by hisat:
                        req_type = [
                            [('-1', '-2'), ('.fq', '.fa'), True],
                            [('-U', ),     ('.fq', '.fa'), True],
                            [('-S', ),     ('sam', )],
                        ]
                    The 'match_only' bool element is used only for matching
                    linked command input and is assumed False if ommitted.
                    If True, the option will not be matched unless the
                    given flags are matched exactly. In the above example,
                    a single linked fastq file would be set to '-U' and three
                    linked fastq files would not be used at all.
                    See '_match_input_with_args' & '_check_type' for details.
        n_priority_args
                    The number of arguments to put before keyword arguments
                    and flags. Set to -1 to put ALL args first. Default: 0.
        allow_bash_var
                    A bool indicating whether the command should allow
                    BASH-stype variables. If False, any variables will be
                    stripped: ${var_name} => var_name. Default: False.
                    It is not recommended to allow BASH-style variables
                    unless it is known that they will be used.

    Example:
        hisat_attr = CmdAttributes(
            name='hisat2',
            invoke_str='hisat2'
        )

    '''

    def __init__(
        self,
        name='', synopsis='', description='',
        invoke_str='',
        arguments=[], defaults={},
        req_args=0,
        req_kwargs=[],
        req_type=[],
        n_priority_args=0,
        allow_bash_var=False,
    ):
        if not invoke_str:
            invoke_str = name

        # non-mutable
        self.name = name
        self.synopsis = synopsis
        self.description = description
        self.invoke_str = invoke_str

        # ensure deep copies
        try:
            self.arguments = arguments[:]
        except TypeError:
            log.debug(arguments)
            raise
        self.defaults = dict.copy(defaults)
        self.req_args = req_args
        self.req_kwargs = req_kwargs[:]
        self.req_type = req_type[:]
        self.n_priority_args = n_priority_args
        self.allow_bash_var = allow_bash_var

        if not arguments:
            raise ValueError('No arguments given')

    def duplicate(self, **kwargs):
        '''Deep copy the current object and return the copy.

        Arguments:
            **kwargs: Values used to update the duplicate.
        Return:
            A new CmdAttributes object
        '''

        new_attr = CmdAttributes(**self.__dict__)
        filt_dict = {k: v for k, v in kwargs.items() if k in self.__dict__}
        new_attr.__dict__.update(filt_dict)

        return new_attr

    def get_types(self, flag):
        for rt in self.req_type:
            if flag in rt[0]:
                return rt[1]
        return None


class BaseCmd(CmdInterface):

    '''An abstract class for handling command line calls.

    Create a command line string for executing a given program. Children
    must define certain class attributes that characterize the program
    requirements. Certain methods may be overridden to allow specific
    handling and preparation of given arguments.

    TODO:           Implement command attributes using the Strategy pattern.
    TODO:           Allow ARGUMENTS to accept actual types to allow for
                    better type checking (instead of just strings)

    Attributes:
        strict: A boolean indicating whether to use strict arg handling.
            If true, cmd and init will raise for unknown kwargs and
            unused link input. Default: True.
        timestamp: A string indicating the time at initialization.
        wrap: An optional string variable name. If set, the command will be
            wrapped in a BASH-style assignment to a var with the given name.
            Default: None.
        fall_through: A boolean indicating whether or not the given input
            should be passed as ("fall through") output. Default: False.

    Interface:
        output(): Return a list of output files
        link(CMD): Link the given CMD to use the current CMD's output.
            Returns the given CMD to allow chaining.
        cmd(): Return the BASH-style command string. The majority
            of requirement checking is done here to allow for
            more dynamic creation of the command string.

    Class Attributes:
        See CmdAttributes class.

    Class Methods for override (in call order):
        _prepreq(): Use to setup any last minute details before
            requirements are checked, such as handling the
            results of linking.
            NOTE: Call parent to ensure intended functionality.
        _additional_requirements(): Use to define additional requirements
            related to successful command execution. This should raise
            an exception if conditions are not met.
            NOTE: This method is not implemented in BaseCmd.
        _prepcmd(): Use to setup any last minute details before cmd
            is called, such as setting up redirect.
            NOTE: This method is not implemented in BaseCmd.

    Custom Exceptions:
        BaseCmd defines a number of custom exceptions with pre-defined
        error messages. These messages are meant for internal use ONLY
        and can be caught w/ specific Python exceptions, e.g., ValueError.
        Using custom exceptions over custom messages helps to ensure that
        the error type remains constant for a given type of error:
            CmdLinkError -> TypeError
            PositionalArgError -> IndexError
            KeywordArgError -> AttributeError
            FileTypeError -> TypeError
            ConfigError -> ValueError

    Custom Configuration:
        Some commands may allow for the use of files that contain
        configuration data for running the command. See the specific
        child for use case.

    SEE ALSO:
        HelpCmd, CmdAttributes
    '''

    # NOTE on prefix commands:
    #     [SJB] I had once considered allowing for a sort-of sub command
    #     whose input was used directly by the following command, e.g.,
    #         k=$(velvetk <args>)
    #         velveth out_dir $k <args>
    #     However, the obvious issue is linking the prefix command to the
    #     current command. The linking proves fairly complicated with
    #     several use cases, particularly linking variable names, position
    #     index, or keyword.
    #
    #     Instead, I think there are two simpler options:
    #         1. Wrapping a command. Pass in a variable name, and wrap the
    #             command in a bash variable assignment: X="$(<cmd>)"
    #         2. Fall through. A command with fall through will pass ALL
    #             input into output (while optionally adding its own output).
    #
    #     Pros:
    #     > Less specialization needed. Wrapping is simple and globally
    #         configurable. Fall through is obvious.
    #     > Presumably, the argument that would use the prefix may still
    #         be used normally (given the value directly rather than relying
    #         upon the prefix). This simplifies that.
    #     > Presumably, commands that would use a prefix are very closely
    #         tied to the prefix--these are not global or common situations.
    #         As such, handling via a specific pipe would be better.
    #         (Pass responsibility up the chain)
    #
    #     Cons:
    #     > Requires more specific attention for each particular command.
    #         Loss of modularity.

    #
    #   Class attributes
    #

    # Children must set 'attr' to CmdAttributes instance
    attr = None

    #
    #   Custom Exceptions
    #

    class CmdLinkError(RempipeError, TypeError):
        ERRMSG = {
            'input': 'Bad link: no input given or input not callable',
            'mismatch': 'Bad link: output not matched to input',
        }

    class ArgError(RempipeError, ValueError):
        ERRMSG = {
            'type': 'Bad argument type. "{} != {}"',
        }

    class PositionalArgError(RempipeError, IndexError):
        ERRMSG = {
            'missing': 'Missing {} of {} required positional parameters',
            'str_only': ('Positional args should be strings. ' +
                         '(Did you expand **kwargs?)'),
            'unknown': 'Unknown positional parameters.',
        }

    class KeywordArgError(RempipeError, AttributeError):
        ERRMSG = {
            'missing': '{}: Missing required keyword arguments: {}',
            'unknown': 'Unrecognized keyword argument given: "{}"',
            'xor': 'More than 1 XOR kwarg given: {}',
        }

    class FileTypeError(RempipeError, TypeError):
        ERRMSG = {
            'missing': 'Expected file or file type not found',
        }

    class ConfigError(RempipeError, ValueError):
        ERRMSG = {
            'illegal': 'Illegal characters in config file: {}',
        }

    #
    #   Magic methods
    #

    def __init__(
            self, *args,
            strict=True, timestamp=None, wrap=None, fall_through=False,
            ** kwargs):
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
            raise self.PositionalArgError('str_only')
        args = [v for v in args if v not in flags]

        if strict:
            # check for unknown kwargs
            # NOTE: We MUST have removed timestamp AND
            #       ensured hyphens BEFORE this
            to_check = flags + list(kwargs.keys())
            self._check_for_unknown_flags(to_check)
        self.strict = strict

        # make sure we deep copy defaults and args
        self.redirect = None
        self.kwargs = dict.copy(self.attr.defaults)
        self.kwargs.update(kwargs)
        self.args = args[:]
        self.flags = flags[:]

        # setup wrap
        self.wrap = wrap  # wrap command in bash assignment statement

        # setup fall through
        if fall_through:
            self.output = fall_through_decorator(self.output)

    def __str__(self):
        '''Returns the BASH executable string'''
        return self.cmd()

    #
    #   Property methods
    #

    @property
    def name(self):
        '''Name property. Do NOT override'''
        return self.attr.name

    @property
    def invoke_str(self):
        '''Invocation string property. Do NOT override (unless necessary)'''
        return self.attr.invoke_str

    #
    #   "Public" access
    #

    @classmethod
    def help(cls):
        '''Return usage text

        NOTE: Creates the HelpCmd object iff not found
        '''
        try:
            return str(cls._help)
        except AttributeError:
            cls._help = HelpCmd(**cls.attr.__dict__)
            return str(cls._help)

    @abc.abstractmethod
    def output(self):
        '''Return list of created files.

        Children must override.
        Should at the very least return the input to allow easier chaining.
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

    def cmd(self, *args, wrap=None, **kwargs):
        '''Generate the command string

        Build the command string using the given arguments/inputs.
        Performs a number of processing steps before actually creating the
        command string:
            1. Match up the given input with arguments (uses attr.req_type).
            2. Process arguments as needed BEFORE checking requirements.
                > via '_prepreq', defined per child class.
                > e.g., ensure correct argument order.
            3. Check that requirements are met (via attr.req_[kw]args).
            4. Check for any additional requirements.
                > via '_additional_requirements', defined by child.
                > e.g., ensure a complete Bowtie index exists.
            5. Process arguments as needed AFTER checkeing requirements.
                > via '_prepcmd', as defined per child class.
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

        try:
            self._match_input_with_args()
        except self.CmdLinkError:
            raise  # TypeError!

        # run requirements prep, if provided by child
        # -- use for ensuring the requirements are met
        try:
            self._prepreq()
        except NotImplementedError:
            pass  # may not be set on child class

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
        except NotImplementedError:
            pass  # ignore if it's not defined

        # run command prep, if provided by child
        # -- use to ensure all data is ready for the command
        # -- NOTE: should require the necessary arguments
        try:
            self._prepcmd()
        except NotImplementedError:
            pass  # may not be set on child class

        cmd_str = self._cmd(*args, **kwargs)

        # wrap command
        # -- from command takes precedence from init
        if not wrap:
            wrap = self.wrap
        if wrap:
            cmd_str = '{}="$({})"'.format(wrap, cmd_str)

        # join all current parts with a new line character
        return cmd_str

    #
    #   To override
    #

    def _prepreq(self):
        raise NotImplementedError('Implement on child')

    def _prepcmd(self):
        raise NotImplementedError('Implement on child')

    def _additional_requirements(self):
        raise NotImplementedError('Implement on child')

    #
    #   "Private" access
    #

    def _has_output(self):
        '''Checks for existence of specified output files.

        Arguments:
            None.
        Returns:
            True if ALL output files exist.
            False otherwise.
        '''

        for f in self.output():
            if not os.path.isfile(f):
                return False
        return True

    def _match_input_with_args(self):
        '''Match linked inputs to arguments according to attr.req_type

        NOTE: If NOT strict AND the required number of positional args has
            not been met, linked input will be assigned to args. Any
            positional arguments set in attr.req_type will be set
            regardless of attr.req_args if the expected type is found.

        Uses attr.req_type to match linked input to a keyword or
        positional argument, e.g.,
            input: ['a.txt', 'a.fq', 'b.fq']
            attr.req_args = 1
            attr.req_type = [[('-b', '-a'), ('.fq',)]]  <-- notice the order!
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
        except TypeError:
            raise self.CmdLinkError('input')

        # match with required type
        for req_type in self.attr.req_type:
            try:
                flag_list, type_list, match_only = req_type
            except ValueError:
                flag_list, type_list = req_type
                match_only = False

            # find args with matching file type and update kwargs
            filtered = self._filter_by_type(args, type_list)
            if not match_only or len(filtered) == len(flag_list):
                self.kwargs.update(dict(zip(flag_list, filtered)))
            else:
                continue  # we don't want to remove if we didn't use!!

            # remove used
            args = [a for a in args if a not in filtered]

        # if a positional argument was defined in attr.req_types,
        # it is now stored in kwargs -- move it to args and clean kwargs
        pos_kw = [k for k in self.kwargs.keys() if isinstance(k, int)]
        self.args.extend([self.kwargs[pos] for pos in sorted(pos_kw)])
        for pos in pos_kw:
            del self.kwargs[pos]  # delete instead of dict comp

        # otherwise, use as required args IFF expected
        # i.e., attr.req_args > 1
        # BUT only if NOT strict
        if not self.strict and len(self.args) < self.attr.req_args:
            n_still_required = self.attr.req_args - len(self.args)
            self.args.extend(args[:n_still_required])

        # finally, if we haven't used ANY of the input, raise CmdLinkError
        # BUT, only if we're in strict
        elif len(args) == n_args:
            if self.strict:
                raise self.CmdLinkError('mismatch')
            else:
                log.warning(self.CmdLinkError.ERRMSG['mismatch'])

        return

    def _cmd(self, *args, verbose=True, strict=True, **kwargs):
        '''Create BASH executable string.

        Arguments:
            verbose: A bool indicating whether the command should be easier
                to read. Currently, that means splitting the command with
                new lines.
            strict: CRITICAL: Double check this. Should use from init?
                A bool indicating whether or not to use strict cleaning
                of the command string before returning.
        Returns:
            A BASH executable string.
        '''

        # put the command on separate lines for human readable version
        sep = ' \\\n  ' if verbose else ' '

        # make strings from given parameters
        flags = ' '.join(self.flags)  # flags don't need to be separated
        kwargs = sep.join(
            "{} {}".format(k, v)
            for k, v in sorted(self.kwargs.items())
        )

        # setup priority arguments (args that should come BEFORE kwargs)
        priority = (len(self.args) if self.attr.n_priority_args == -1
                    else self.attr.n_priority_args)
        priority_args = sep.join(self.args[:priority])
        args = sep.join(self.args[priority:])

        # redirect may be a tuple ('>', 'logfile.log') or string
        if not isinstance(self.redirect, str) and self.redirect is not None:
            redirect = ' '.join(self.redirect)
        else:
            redirect = self.redirect

        # NOTE: redirect WILL contain unsafe chars ('>' or '|', likely)
        #       We should account for this later...or prevent redirect
        #       from modification
        cmd_parts = filter(  # remove missing elements
            None, [self.invoke_str, priority_args, flags, kwargs, args])

        cmd = sep.join(cmd_parts)
        cmd_safe = self._unsafe_char_protect(cmd, strict=strict)

        if redirect:
            cmd_safe = '{}{}{}'.format(cmd_safe, sep, redirect)

        return cmd_safe

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
                raise self.KeywordArgError('missing', details=[
                    self.name, ', '.join(str(m) for m in missing)
                ])

        # check for expected number of args
        if len(self.args) < self.attr.req_args:
            raise self.PositionalArgError('missing', details=[
                self.attr.req_args - len(self.args), self.attr.req_args,
            ])

        # check for expected file types
        # NOTE: Do this LAST to avoid unnecessary errors
        try:
            self.__check_type()
        except self.FileTypeError:
            raise  # Bad file type

        return

    def __check_type(self):
        '''Check specified kwargs and args to ensure given expected data type.

        NOTE: All extensions in attr.req_type should have leading periods.
        TODO: Update to use Django model to handle file types.

        Returns:
            True if expected file types match given files.
        Raises:
            ValueError if not expected data type
            FileTypeError if a given argument has wrong file type.
        '''

        for req in self.attr.req_type:
            try:
                # match_any used for matching only, not checking
                flags, types, match_any = req
            except ValueError:
                flags, types = req

            # convert "flags" to values
            # NOTE: the "flags" may be arg indices, too
            vals = self.__get_values_from_flags(flags)
            if not vals:
                continue  # no relevant values

            # check for explicit types, e.g., int
            try:
                self.__check_basic_type(vals, types)
            except self.ArgError:
                raise  # invalid type!

            # Expected error: 'str' object not callable
            # check vals for file type
            except TypeError:
                try:
                    self.__check_file_type(vals, types)
                except AttributeError:
                    # the given object is probably not a string
                    # -- mising 'rfind'
                    raise

        return True

    def __get_values_from_flags(self, flags):
        '''Return a list of values for the given flags

        Arguments:
            flags: A list of kwarg flags or arg positions
        Return:
            A list of values for the given flags
        '''

        indices = [f for f in flags if isinstance(f, int)]
        flags = [f for f in flags if f not in indices]
        vals = [
            self.kwargs[flg]
            for flg in flags
            if flg in self.kwargs
        ]
        vals.extend(
            self.args[index] for index in indices
            if index < len(self.args)
        )

        return vals

    def __check_basic_type(self, vals, types):
        '''Check given values have expected data type

        Check that given values can be coerced into specific types.
        Compares string value after conversion to ensure no information loss.
        For example, int(0.5) is a valid statement that returns 0, but
        if the expected type is int, 0.5 is NOT a valid value.
        Unfortunately, the value may not be given as an int, so we cannot
        rely on isinstance. Therefore, we convert the final result to a
        string and compare the string of the original value:
            str(int(0.5)) == str(val)
        And raise a ValueError if they're not equal.

        NOTE: Profile efficiency of this function.
        NOTE: Will not handle multiple types well. For example,
            [int, float] will raise for [1.1].

        Arguments:
            vals: A list of values to check
            types: A list of type callables, eg., int
        Return:
            None
        Raises:
            ArgError [ValueError] if the value is not of the expected type.
            TypeError if the type is not callable.
        '''

        for val in vals:
            for typ in types:
                try:
                    if str(typ(val)) != str(val):
                        raise ValueError
                except ValueError:
                    raise self.ArgError('type', details=[val, typ])
                except TypeError:
                    raise

    @classmethod
    def __check_file_type(cls, vals, types):
        '''Check given values to ensure given expected file type.

        NOTE: All extensions in attr.req_type should have leading periods.
        TODO: Update to use Django model to handle file types.

        Arguments:
            vals: A list of argument values
            types: A list of file extensions
        Returns:
            None
        Raises:
            FileTypeError [TypeError] if a given val has wrong file type.
        '''

        for v in vals:
            try:
                extn = os.path.splitext(v)[1]
            except:
                raise

            if extn not in types:
                msg = 'Invalid file type "{}" given for argument {}'
                raise cls.FileTypeError(msg.format(extn, v))

    def __check_kwargs(self):
        '''Make sure that the given kwargs contain all required kwargs

        Checks all elements in attr.req_kwargs.
        > Strings are simple checks against the kwargs attribute.
        > Lists are XOR required keywords, e.g., one and only one req.
        > Tuples are AND required keywords, e.g., if one, then all.

        Example attr.req_kwargs:
            ['-f', ('-1', '-2'), ['-1', '-U']]

        Arguments:
            A kwargs dict.
        Returns:
            A list of missing kwargs.
        '''

        # simple requirement
        simple = [kw for kw in self.attr.req_kwargs if isinstance(kw, str)]
        missing = [kw for kw in simple if kw not in self.kwargs]

        # compound requirement -- require all or none of args in the tuple
        compound = [kw for kw in self.attr.req_kwargs if kw not in simple]
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
                    raise self.KeywordArgError('xor', details=cmpd)

        return missing

    #
    #   Private access class methods
    #

    @classmethod
    def _filter_by_type(self, args, filt):
        '''Filter a list of strings based on a list of extension strings

        Arguments:
            args: A list of file path strings
            filt: A list of file extension strings
        Return:
            A list of file path strings that match the given extensions.
        Example:
            self._filter_by_type(['a.txt', 'b.fq'], ['.fq'])
        '''
        filtered = [
            arg for arg in args
            if os.path.splitext(arg)[1] in filt
        ]
        return filtered

    @classmethod
    def _trubase(self, path_name):
        '''A shortcut for splitting the path and the extension'''

        return os.path.splitext(os.path.basename(path_name))[0]

    @classmethod
    def _unsafe_char_protect(cls, input_str, strict=True):
        '''Removes unsafe characters from a string

        Checks for unicode control characters and unsafe BASH chars and
        removes them from the given string. If not strict, unsafe BASH
        characters are escaped instead. Uses pre-compiled regex for speed.

        Arguments:
            input_str: The string to clean.
            strict: A bool indicating whether to remove or escape unsafe
                BASH characters.
        Return:
            The clean string.
        '''

        try:
            # protect BASH variable declarations
            if cls.attr.allow_bash_var:
                input_str = cls.rx_bash_var.sub(
                    r'BASHOPEN\g<var>BASHCLOSE', input_str)

            # remove unicode control characters
            input_str = cls.rx_illegal_char.sub('', input_str)

            # remove BASH unsafe characters
            replacement = 'ESCCHAR\g<unsafe>' if not strict else ''
            input_str = cls.rx_unsafe_char.sub(replacement, input_str)
            if replacement:
                input_str = input_str.replace(replacement[:7], '\\')

            # reinstate bash variable declarations
            if cls.attr.allow_bash_var:
                input_str = input_str.replace('BASHOPEN', '${')
                input_str = input_str.replace('BASHCLOSE', '}')
            return input_str.rstrip()

        except AttributeError:
            # compiile and create these regex only once
            cls.__init_bash_var_regex()
            cls.__init_unsafe_regex()
            cls.__init_illegal_char_regex()
            return cls._unsafe_char_protect(input_str, strict=strict)

    @classmethod
    def __init_bash_var_regex(cls):
        '''Compile the regex for bash variable interpolation
        NOTE: The braces are intentionally required.
        '''
        cls.rx_bash_var = re.compile('\$\{(?P<var>\w+)\}')

    @classmethod
    def __init_unsafe_regex(cls):
        '''Compile the unicode control char regex'''
        unsafe = re.escape(';&|><*?`$(){}[]!#')
        unsafe = r'(?P<unsafe>[{}]\s?)\s*'.format(unsafe)
        cls.rx_unsafe_char = re.compile(unsafe)

    @classmethod
    def __init_illegal_char_regex(cls):
        '''Compile the unsafe BASH char regex'''
        # From a blog post for removing illegal ASCII control characters
        # http://chase-seibert.github.io/blog/2011/05/20/stripping-control-characters-in-python.html
        RE_XML_ILLEGAL = (
            u'([\u0000-\u0008\u000b-\u000c\u000e-\u001f\ufffe-\uffff])' +
            u'|' +
            u'([%s-%s][^%s-%s])|([^%s-%s][%s-%s])|([%s-%s]$)|(^[%s-%s])' %
            (chr(0xd800), chr(0xdbff), chr(0xdc00), chr(0xdfff),
             chr(0xd800), chr(0xdbff), chr(0xdc00), chr(0xdfff),
             chr(0xd800), chr(0xdbff), chr(0xdc00), chr(0xdfff),
             ))
        cls.rx_illegal_char = re.compile(RE_XML_ILLEGAL)

    @classmethod
    def _check_for_unknown_flags(cls, flags):
        '''Check for unknown flags and kwargs

        Arguments:
            flags: A list of flags, e.g., ['-f', '--verbose'].
        Return:
            None.
        Raises:
            KeywordArgError if an unknown flag is found.
        '''
        try:
            known_flags = [v[0] for v in cls.attr.arguments if v[0]]
        except AttributeError:
            raise NotImplementedError(
                'attr.arguments must be set on child class')
        unknown_flags = [
            k for k in flags
            if k not in known_flags
        ]
        if unknown_flags:
            raise cls.KeywordArgError('unknown', details=unknown_flags)

    @classmethod
    def _ensure_file_and_extension(cls, file_name, extn_list):
        '''Check that a given file or file prefix has an extension

        Check that the filename ends with one of the expected extensions.
        If not, assume the file_name is actually a file prefix and look
        for existing files matching the prefix + extension.

        Arguments:
            file_name   The file name or file prefix.
            extn_list   A list of expected file endings.
        Return:
            The filename with the valid extension (unchanged if already
            had an expected extension).
        Raises:
            FileTypeError if no extension given and no existing file found.
        '''

        # check that filename ends with given list
        # INSTEAD of checking filename extension against list
        # -- provides some flexiblity, e.g., checking for '_1.fasta'
        file_extn = [
            extn for extn in extn_list
            if file_name.endswith(extn)
        ]

        # the file has one of the expected extensions
        if file_extn:
            return file_name

        # check for presence of file with expected extension
        # -- assume the given file_name is a prefix!
        # -- update with first found match
        # NOTE: LOGIC ERROR! If multiple exist, e.g., both .bed and .gff,
        #       only the first found will be used.
        for extn in extn_list:
            if os.path.isfile(file_name + extn):
                file_extn = extn
                # self.kwargs['-bed'] = bed_file + extn
                break

        # update the file name
        try:
            file_name = file_name + file_extn

        # file_extn not updated in loop (still a list from before)
        except TypeError:
            raise cls.FileTypeError('missing')

        # otherwise, return the new filename
        else:
            return file_name
