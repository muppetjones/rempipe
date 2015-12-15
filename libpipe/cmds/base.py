import logging
log = logging.getLogger(__name__)


class BaseCmd(object):

    def __init__(self, *args, **kwargs):

        # check for required kwargs
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
                'Missing positional parameters; expected {}, given {}'.format(
                    self.required_args, len(args)
                ))

        self.bin = self.bin_name
        self.kwargs = self.defaults
        self.kwargs.update(kwargs)
        self.args = args

        # flags and redirect options must be set on a per software basis
        self.flags = []
        self.redirect = ''

    def __str__(self):
        return self.cmd()

    def cmd(self):
        flags = ' '.join(self.flags)
        kwargs = ' '.join(
            "-{} {}".format(k, v)
            for k, v in self.defaults.items()
        )
        args = ' '.join(self.args)
        return '{} {} {} {} {}'.format(
            self.bin,
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

    def output(self):
        '''Return list of created files. Define on a software level'''
        return None
