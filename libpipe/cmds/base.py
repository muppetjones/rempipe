from abc import ABCMeta

import logging
log = logging.getLogger(__name__)


class BaseCmd(metaclass=ABCMeta):

    # flags and redirect options must be set on a per software basis
    flags = []
    redirect = ''

    # required options must also be set per sofware
    required_kwargs = {}
    required_args = 0

    def __init__(self, *args, **kwargs):

        # ensure the expected hyphens are in place
        # NOTE: we need to do this BEFORE checking requirements
        #       to ensure the given kwargs match the expected kwargs
        def ensure_hyphen(flag):
            if not flag.startswith('-'):
                fmt = '-{}' if len(flag) == 1 else '--{}'
                flag = fmt.format(flag)
            return flag
        kwargs = {ensure_hyphen(k): v for k, v in kwargs.items()}

        # check for required kwargs
        try:
            missing = self._check_kwargs(kwargs)
            if missing:
                raise ValueError('Missing arguments:\n\t{}\n'.format(
                    '\n\t'.join(missing)
                ))

            # check for expected number of args
            if len(args) < self.required_args:
                raise ValueError(
                    'Missing {} of {} positional parameters; '.format(
                        self.required_args - len(args), self.required_args
                    ))
        except AttributeError:
            pass  # required_kwargs or required_args not set

        self.bin = self.bin_name
        try:
            self.sub = self.sub_cmd
        except AttributeError:
            pass  # not all commands have a sub command
        self.kwargs = {}
        self.kwargs.update(self.defaults)
        self.kwargs.update(kwargs)
        self.args = args

    def __str__(self):
        return self.cmd()

    @property
    def name(self):
        return self.bin_name

    @property
    def output(self):
        '''Return list of created files. Define on a software level'''
        return None

    def _check_kwargs(self, kwargs):
        '''Make sure that the given kwargs contain all required kwargs

        Arguments:
            A kwargs dict.
        Returns:
            A list of missing kwargs.
        '''
        simple = [kw for kw in self.required_kwargs if isinstance(kw, str)]
        compound = [kw for kw in self.required_kwargs if kw not in simple]

        missing = [kw for kw in simple if kw not in kwargs]

        # compound requirement -- require all or none of args in the tuple
        for cmpd in compound:
            missed = [kw for kw in cmpd if kw not in kwargs]
            if missed and len(missed) != len(cmpd):
                missing.extend(missed)

        return missing

    def cmd(self):
        try:
            command = "{} {}".format(self.bin, self.sub)
        except AttributeError:
            command = self.bin
        flags = ' '.join(self.flags)
        kwargs = ' '.join(
            "{} {}".format(k, v)
            for k, v in self.kwargs.items()
        )
        args = ' '.join(self.args)
        return '{} {} {} {} {}'.format(
            command,
            flags,
            kwargs,
            args,
            self.redirect,
        )

    def help(self):

        try:
            print(self.description)
        except AttributeError:
            pass

        try:
            print("Attributes")
            print("\n".join(
                "\t{}\t{}".format(k, v)
                for k, v in self.attributes.items()
            ))
        except AttributeError:
            pass
