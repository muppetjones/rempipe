import os.path
from abc import ABCMeta

import logging
log = logging.getLogger(__name__)


class BaseCmd(metaclass=ABCMeta):

    # command defining attributes
    NAME = None  # human readable name of command, e.g., 'samtools_sort'
    INVOKE_STR = None  # string to invoke command, e.g., 'samtools sort'
    DEFAULTS = {}  # starting values for kwargs
    HELP_DICT = {  # help text
        'synopsis': 'cmd_invoke -f FILE [-n INT, -h] FILE ...',
        'description': 'A description of the command.',
        'positional arguments': [
            ('FILE', 'Example input one'),
            ('FILE', 'Example input two'),
        ],
        'keyword arguments': [
            ('-f|--file', 'FILE', 'Example keyword argument'),
            ('-n', 'INT', 'Example keyword argument'),
        ],
        'flags': [
            ('-v', 'Example flag argument'),
        ],
    }

    # these parameters should be set during init
    # NOTE: redirect (& possibly flags) should be set at a class level
    FLAGS = []
    KWARGS = {}
    ARGS = []
    REDIRECT = ''

    # set argument requirements
    # NOTE: required kwargs is a list of the required argument flags
    # NOTE: required args is an int of the number of required args
    # NOTE: flags, by definition, are optional
    REQ_KWARGS = []
    REQ_ARGS = 0

    def __init__(self, *args, **kwargs):
        '''Initialize command object

        Create a command object using given paramters.
        Flags, e.g., '-v', should be given as positional parameters.
        '''

        # ensure the expected hyphens are in place
        # NOTE: we need to do this BEFORE checking requirements
        #       to ensure the given kwargs match the expected kwargs
        def ensure_hyphen(flag):
            if not flag.startswith('-'):
                fmt = '-{}' if len(flag) == 1 else '--{}'
                flag = fmt.format(flag)
            return flag
        kwargs = {ensure_hyphen(k): v for k, v in kwargs.items()}

        # separate flags from args
        flags = [v for v in args if v.startswith('-')]
        args = [v for v in args if v not in flags]

        # check for required kwargs
        try:
            missing = self._check_kwargs(kwargs)
            if missing:
                raise ValueError('Missing arguments:\n\t{}\n'.format(
                    '\n\t'.join(missing)
                ))

            # check for expected number of args
            if len(args) < self.REQ_ARGS:
                raise ValueError(
                    'Missing {} of {} positional parameters; '.format(
                        self.REQ_ARGS - len(args), self.REQ_KWARGS
                    ))
        except AttributeError:
            pass  # required_kwargs or required_args not set

        # make sure we deep copy defaults and args
        self.redirect = self.REDIRECT

        self.kwargs = {}
        self.kwargs.update(self.DEFAULTS)
        self.kwargs.update(kwargs)

        self.args = []
        self.args.extend(args)

        self.flags = []
        self.flags.extend(flags)

    @property
    def name(self):
        '''Name property. Do NOT override'''
        return self.NAME

    @property
    def invoke_str(self):
        '''Invocation string property. Do NOT override (unless necessary)'''
        return self.INVOKE_STR

    def __str__(self):
        return self.cmd()

    @property
    def output(self):
        '''Return list of created files.

        Must override. Should at the very least return the input to allow
        easier chaining.
        '''
        return None

    @classmethod
    def _trubase(self, path_name):
        return os.path.splitext(os.path.basename(path_name))[0]

    def _check_kwargs(self, kwargs):
        '''Make sure that the given kwargs contain all required kwargs

        Arguments:
            A kwargs dict.
        Returns:
            A list of missing kwargs.
        '''
        simple = [kw for kw in self.REQ_KWARGS if isinstance(kw, str)]
        compound = [kw for kw in self.REQ_KWARGS if kw not in simple]

        missing = [kw for kw in simple if kw not in kwargs]

        # compound requirement -- require all or none of args in the tuple
        for cmpd in compound:
            missed = [kw for kw in cmpd if kw not in kwargs]
            if missed and len(missed) != len(cmpd):
                missing.extend(missed)

        return missing

    def cmd(self, readable=True):
        sep = ' \\\n  ' if readable else ' '

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

        # 'example', 'args', 'kwargs', 'flags']
        help_order = [
            'synopsis', 'description',
            'positional arguments', 'keyword arguments', 'flags']

        sec_txt = ["NAME\n\t{}\n".format(self.NAME), ]

        for sec in help_order:
            if isinstance(self.HELP_DICT[sec], list):
                txt = '\n\t'.join([
                    '\t'.join(t) for t in self.HELP_DICT[sec]])
            else:
                txt = self.HELP_DICT[sec]
            sec_txt.append('{}\n\t{}\n'.format(
                sec.upper(), txt
            ))

        sec_txt = "\n".join(sec_txt)
        return sec_txt
