import logging
log = logging.getLogger(__name__)


class BaseCmd(object):

    def __init__(self, *args, **kwargs):

        # check for required kwargs
        try:
            missing_req = []
            for req in self.required_kwargs:
                if isinstance(req, tuple):
                    missing_dual_req = [r for r in req if r not in kwargs]
                    if missing_dual_req:
                        missing_req.append(req)
                elif req not in kwargs:
                    missing_req.append(req)
            if missing_req:
                raise ValueError(
                    'Missing required parameters: {}'.format(missing_req))

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
        log.debug(self.args)

        # flags and redirect options must be set on a per software basis
        self.flags = []
        self.redirect = ''

        # required options must also be set per sofware
        self.required_kwargs = {}
        self.required_args = 0

    def __str__(self):
        return self.cmd()

    @property
    def name(self):
        return self.bin_name

    @property
    def output(self):
        '''Return list of created files. Define on a software level'''
        return None

    def cmd(self):
        try:
            command = "{} {}".format(self.bin, self.sub)
        except AttributeError:
            command = self.bin
        flags = ' '.join(self.flags)
        kwargs = ' '.join(
            "-{} {}".format(k, v)
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
