

import logging
log = logging.getLogger(__name__)


class CmdAttributes(object):
    '''Store attributes for executing a program on the command line

    NOTE: This class should NOT be called directly. Use libpipe.cmd.Cmd.
    NOTE: Be sure to read the expected formats!

    Store various attributes directly related to calling a particular
    program from the command line. Given values are deep copied to prevent
    unexpected changes.

    The attributes have a specific format to enable JSON import and export.
    Although detailed, it *should* be fairly intuitive. At a later point,
    it may make sense to implement an 'argparse' style interface.

    TODO: Enable JSON import and export.
    TODO: Enable usage printing on error.
    TODO: Enable natural language processing for 'req_kwargs'
    TODO: Enable custom cmd attribute population based on help string.

    Attributes:
        invoke: [REQUIRED] A string indicating how to call the program.
            Ex: 'tophat2'.
            NOTE: Be sure to include the path if necessary.
        args: [REQUIRED] A list of tuples indicating the expected parameters,
            where each element should have the following format:
                (<flag>, <type>, <description>)
            A <flag> value of None or an integer indicates a positional
            argument; an integer should be used if the argument value needs
            to be checked (see 'req_types').
            A <type> value of None indicates a flag option. The value
            will only be used for documentation purposes.
            Likewise, the <description> is only used for documentation.
            Examples:
                (0, 'FILE',  'The input file'),
                ('-h', None,    'Display the help text'),
                ('-o', 'FILE',  'The output file')

        name: A human readable string of the program name. Ex: 'Tophat 2'.
            Set to the invoke string if not given.
        synopsis: A string detailing usage.
        description: A string describing the program.

        req_pargs: An integer indicating no. of required positional args.
        req_kwargs: A list of flags for required keyword args. List elements
            should be strings (individual flags), tuples (mutually required),
            or lists (exclusively required).
            Example:
                # define basic bowtie / hisat requirements
                # -x is required
                # if -1, then -2 required
                # if -1 or -U, then not the other
                [ '-x', ('-1', '-2'), ['-1', '-U'], ]
        req_types: A list of args and their expected extensions.
            A single requirement will have 2-3 elements:
                <flags>, <types>, [exact_match]
            where <flags> is a tuple listing the string flags or integer
            positions, <types> is a tuple of extension strings OR callable
            types. The <exact_match> element is optional (default: False)
            and indicates whether the each of the given flags must meet
            the requirements.

            NOTE: File extensions should include the period.

            For example, requirements defined as follows:
                [('-1', '-2', '-U'), ('.fq', '.fastq'), ]
            will require that any of the following have a FASTQ extension.
                [('-1', '-2'),       ('.fq', '.fastq'), True],
            will require that BOTH '-1' and '-2' have a FASTQ extension.
            In addition, this restricts input to exactly 2 FASTQ files.
            Example (hisat/bowtie):
                # Accept either one (-U) or two (-1, -2) FASTQ files,
                # but no more.
                req_type = [
                    [('-U', ),     ('.fq', '.fastq'), True],
                    [('-1', '-2'), ('.fq', '.fastq'), True],
                    [('-S', ),     ('sam', )],
                ]
        defaults: A dict containing the default values to use for a given
            keyword argument. Should NOT accept positional or flag-only args,
            but not checked.
    '''

    def __init__(self, *args, strict=True, **kwargs):

        required = ['invoke', 'args']
        expected = [
            'name', 'synopsis', 'description',  # help related
            'req_args', 'req_kwargs', 'req_types',  # argument checking
            'defaults', 'flag_sep',
        ]

        # Ensure at least the basic information was given
        missing = [arg for arg in required if arg not in kwargs]
        if missing:
            msg = 'Missing required args: {}'.format(', '.join(missing))
            raise ValueError(msg)

        # Ensure the defaults only provide values for known arguments
        # -- More of a ID10T/typo check.
        try:
            known = [arg[0] for arg in kwargs['args'] if arg is not None]
            unknown = [
                k for k in kwargs['defaults']
                if k not in known
            ]
            if unknown:
                msg = 'Defaults given for unknown args: {}'.format(
                    ', '.join(unknown))
                raise ValueError(msg)
        except KeyError:
            pass  # don't require 'defaults' as input

        # Define allowable attributes
        # -- pseudo deep copy dicts and lists (not their elements)
        # -- Code smell? It isn't pretty.
        allowed_attr = {
            k: (dict.copy(v) if isinstance(v, dict) else
                v[:] if isinstance(v, list) else
                v)
            for k, v in kwargs.items()
            if not strict or k in required or k in expected
        }

        # set as attributes
        self.__dict__.update(allowed_attr)

        if 'name' not in kwargs:
            self.name = kwargs['invoke']

        if 'flag_sep' not in kwargs:
            self.flag_sep = ' '
