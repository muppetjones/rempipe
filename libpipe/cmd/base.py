'''Define the base class for command-line objects'''

import os.path
import time


from libpipe.cmd.interface import CmdInterface

import logging
log = logging.getLogger(__name__)


class CmdBase(CmdInterface):

    '''An abstract class for creating command-line objects

    Create a command line string for executing a given program. Children
    must define certain class attributes that characterize the program
    requirements. Certain methods may be overridden to allow specific
    handling and preparation of given arguments.

    TODO(sjbush): Define file classes (via Django) to handle file type reqs.
    TODO(sjbush): Consolidate and standardize Exceptions and their messages.
    TODO(sjbush): Define 'run' method to execute cmds directly (versus
        relying on a pipe to write to a script).

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

        self.redirect = None
        self.args = _args  # deep copy list--but not elements
        self.flags = _flags
        self.kwargs = dict.copy(self.attr.defaults)
        self.kwargs.update(_kwargs)

        self.strict = strict
        self.complain = complain

    #
    #   Property & magic methods
    #

    def __str__(self):
        return self.cmd(verbose=False)

    @property
    def name(self):
        '''Name property. Do NOT override'''
        return self.attr.name

    @property
    def invoke(self):
        '''Invocation string property. Do NOT override (unless necessary)'''
        return self.attr.invoke

    #
    #   Public access
    #

    def cmd(self, **kwargs):
        '''Generate the command string

        Build the command string using the given arguments/inputs.
        Performs a number of processing steps before actually creating the
        command string:
            1. Match up the given input with arguments (uses attr.req_types).
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

        TODO(sjbush): Implement 'wrap' argument to encase the command as
            a BASH-style variable, e.g., wrap='k' would yield 'k=$(<cmd>)'

        Return:
            A BASH executable string.
        Raises:
            TypeError if link input is bad/unusable or if wrong file type.
            AttributeError if bad kwargs.
            IndexError if bad args.
        '''

        try:
            self._match_input_with_args()
        except AttributeError:
            pass  # not linked...Strange, but ignore
        except TypeError:
            raise  # bad link...let someone know

        self._pre_req()
        self._post_req()
        self._check_requirements()
        self._pre_cmd()
        cmd_str = self._cmd(**kwargs)
        return cmd_str

    def help(self):
        '''Return usage text'''
        raise NotImplementedError()

    def link(self, next_cmd):
        '''Set next_cmd input to self output

        Arguments:
            next_cmd: Another command object.
        Return:
            The given command object
        '''

        next_cmd.input = self.output
        return next_cmd

    def run(self):
        raise NotImplementedError('Cmd cannot be run directly (yet)')

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
            if self.complain:
                log.error(str(e))
            raise

        try:
            self.__check_kwargs()
        except KeyError as e:
            if self.complain:
                log.error(str(e))
            raise

        try:
            self.__check_types()
        except ValueError as e:
            if self.complain:
                log.error(str(e))
            raise

    def _cmd(self, *args, verbose=True, **kwargs):
        '''Create BASH executable string.

        Arguments:
            verbose: A bool indicating whether the command should be easier
                to read. Currently, that means splitting the command with
                new lines.
        Returns:
            A BASH executable string.
        '''

        arg_sep = ' \\\n  ' if verbose else ' '

        # make strings from given parameters
        _args = arg_sep.join(self.args)
        flags = ' '.join(self.flags)  # flags don't need to be separated
        _kwargs = arg_sep.join(
            "{}{}{}".format(k, self.attr.flag_sep, v)
            for k, v in sorted(self.kwargs.items())
        )

        # handle priority args, e.g., velvet
        # -- NOT implemented
        priority_args = None

        # Filter main command elements
        # NOTE: Redirect is a special case--add later
        cmd_parts = filter(  # remove missing elements
            None, [self.invoke, priority_args, flags, _kwargs, _args])

        cmd_str = arg_sep.join(cmd_parts)

        # redirect may be a tuple ('>', 'logfile.log') or string
        # TODO(sjbush): rework 'redirect'
        if not isinstance(self.redirect, str) and self.redirect is not None:
            redirect = ' '.join(self.redirect)
        else:
            redirect = self.redirect
        if redirect:
            cmd_str = '{}{}{}'.format(cmd_str, arg_sep, redirect)

        return cmd_str

    def _match_input_with_args(self):
        '''Match linked inputs to args based on type requirement

        NOTE: If NOT strict AND the required number of positional args has
            not been met, linked input will be assigned to args. Any
            positional arguments set in attr.req_types will be set
            regardless of attr.req_args if the expected type is found.

        Uses attr.req_types to match linked input to a keyword or
        positional argument, e.g.,
            input: ['a.txt', 'a.fq', 'b.fq']
            attr.req_args = 1
            attr.req_types = [[('-b', '-a'), ('.fq',)]]  <-- notice the order!
        will result in the following argument conditions:
            kwargs: {'-b': 'a.fq', '-a': 'b.fq'}
            args: ['a.txt']

        Return:
            None
        Raises:
            ValueError if no link or input is found.
            TypeError if linked input is not callable.
        '''

        try:
            linked_input = self.input()
            n_input = len(linked_input)
        except AttributeError as e:
            if self.complain:
                log.error(str(e))
            raise  # not linked
        except TypeError as e:
            if self.complain:
                log.error(str(e))
            raise TypeError('Bad link')  # not callable
        except Exception as e:
            log.debug(str(e))
            raise

        # We have input...match it up
        # BUT FIRST -- restore from what we "had"
        try:
            self.args = list(self._original_args)
            self.kwargs = dict.copy(self._original_kwargs)
        except AttributeError:
            self._original_args = list(self.args)
            self._original_kwargs = dict.copy(self.kwargs)

        # match with required type
        for req_type in self.attr.req_types:
            log.debug(req_type)
            try:
                flag_list, type_list, exact_match = req_type
            except ValueError:
                flag_list, type_list = req_type
                exact_match = False

            # find args with matching file type and update kwargs
            # -- but only update iff:
            #    a) we're not looking for an exact match, OR
            #    b) we do want an exact match, AND we got it!
            # -- also remove input already matched
            filtered = self._filter_by_type(linked_input, type_list)
            if not exact_match or len(filtered) == len(flag_list):
                self.kwargs.update(dict(zip(flag_list, filtered)))
                linked_input = [a for a in linked_input if a not in filtered]

        # if a positional argument was defined in attr.req_types,
        # it is now stored in kwargs -- move it to args and clean kwargs
        pos_kw = [k for k in self.kwargs.keys() if isinstance(k, int)]
        self.args.extend([self.kwargs[pos] for pos in sorted(pos_kw)])
        for pos in pos_kw:
            del self.kwargs[pos]  # delete instead of dict comp

        # check that some of the input was used
        # -- but we only really care if we're strict or complaining
        if len(linked_input) == n_input:
            msg = 'Unused linked input: {}'.format(', '.join(linked_input))
            if self.strict:
                log.warning(msg)
            if self.complain:
                raise ValueError(msg)
        return

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
        return

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
        for cmpd in compound:

            # tuple indicates all or nothing
            if isinstance(cmpd, tuple):
                missed = [kw for kw in cmpd if kw not in self.kwargs]
                if missed and len(missed) != len(cmpd):  # all missing is fine!
                    missing.append('AND({})'.format(', '.join(missed)))

            # a list indicates exactly one (must have at least one)
            elif isinstance(cmpd, list):
                found = [kw for kw in cmpd if kw in self.kwargs]
                if len(found) != 1:
                    missing.append('XOR({})'.format(', '.join(found)))

        if missing:
            msg = 'Missing required keyword arg: {}'.format(', '.join(missing))
            raise KeyError(msg)
        return

    def __check_types(self):
        '''Check specified kwargs and args to ensure given expected data type.

        NOTE: All extensions in attr.req_types should have leading periods.

        Check arguments for types as defined in CmdAttributes:
            (arg_id, ...), (<type>, ...), [match exact]
        Where 'arg_id' is a kw flag or position index, and type is a
        string file extension or a callable type, e.g., int.

        TODO(sjbush): Enable better parg matching BEFORE checking. The current
            two arg swap is inefficient; plus, would add robustness.

        Raises:
            ValueError if not expected data type
        '''

        for req in self.attr.req_types:
            try:
                arg_ids, types, match_any = req
            except ValueError:
                # match_any used for matching only, not checking...ignore it
                arg_ids, types = req

            # get argument values (pos or kw)
            vals = self.__get_arg_values(arg_ids)
            if not vals:
                continue  # no relevant values

            # check for argument type requirements
            bad_vals = [
                str(val) for val in vals if not self.__check_type(val, types)]

            # basic mercy--try swapping args (if relevant)
            # -- first check if we need to and we have exactly to args
            # -- then make sure that at least one of the args is positional
            # (don't bother swapping back if didn't work!)
            if bad_vals and len(self.args) == 2:
                if any(isinstance(x, int) for x in arg_ids):
                    self.args.reverse()
                    vals = self.__get_arg_values(arg_ids)
                    bad_vals = [
                        str(val) for val in vals
                        if not self.__check_type(val, types)
                    ]

            if bad_vals:
                msg = 'Bad values for [{}]: {}'.format(
                    ', '.join(str(i) for i in arg_ids),
                    ', '.join(bad_vals),
                )
                raise ValueError(msg)

        return

    def __get_arg_values(self, arg_ids):
        '''Return a list of values for the given args

        Arguments:
            arg_ids: A list of kwarg flags or arg positions
        Return:
            A list of values for the given flags
        '''

        indices = [v for v in arg_ids if isinstance(v, int)]
        flags = [v for v in arg_ids if v not in indices]
        vals = [self.kwargs[f] for f in flags if f in self.kwargs]
        vals.extend(self.args[i] for i in indices if i < len(self.args))

        return vals

    #
    #   Class methods
    #   NOTE: Some of these are unavailable to children
    #

    @classmethod
    def __check_type(cls, val, types):
        '''Check that the given value matches the type requirement

        Checks for basic type first, then a file check.

        Args:
            val: The value to test
            types: The types to try
        Returns:
            True if type matched; false otherwise.
        '''

        try:
            try:
                cls.__check_basic_type(val, types)
            except TypeError:
                cls.__check_file_type(val, types)
        except ValueError:
            return False
        return True

    @classmethod
    def __check_basic_type(cls, val, types):
        '''Check given value has expected data type

        Two checks:
            1. value can be coerced into expected type.
            2. coercion does not result in data loss.
        For example, the pair (0.5, int) would pass (1), but not (2).

        TODO(sjbush): Run profiling.

        Arguments:
            val: The value to check
            types: A list of type callables, eg., int
        Raises:
            ValueError if the value is not of the expected type.
            TypeError if the type is not callable.
        '''
        possible_type_error = None
        valid = False
        for _type in types:
            try:
                if str(_type(val)) != str(val):
                    raise ValueError('second check')
            except ValueError as e:
                err = e
            except TypeError as e:
                possible_type_error = e
            else:
                valid = True
                break
        if not valid:
            if possible_type_error:
                raise possible_type_error
            else:
                raise err
        return

    @classmethod
    def __check_file_type(cls, val, types):
        '''Check given value to ensure given expected file type.

        Arguments:
            val: The value to check
            types: A list of file extensions
        Raises:
            ValueError on failure
        '''
        extn = os.path.splitext(val)[1]
        if extn not in types:
            raise ValueError()
        return

    @classmethod
    def _filter_by_type(self, vals, types):
        '''Filter a list of strings based on a list of extension strings

        Arguments:
            vals: A list of file path strings
            types: A list of file extension strings
        Return:
            A list of file path strings that match the given extensions.
        Example:
            self._filter_by_type(['a.txt', 'b.fq'], ['.fq'])
        '''

        file_type = [_type for _type in types if isinstance(_type, str)]
        basic_type = [_type for _type in types if _type not in file_type]

        filtered = [
            val for val in vals
            if type(val) in basic_type
            or (isinstance(val, str) and os.path.splitext(val)[1] in file_type)
        ]
        return filtered
