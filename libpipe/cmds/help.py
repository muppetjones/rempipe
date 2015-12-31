
from textwrap import dedent

import logging
log = logging.getLogger(__name__)


class HelpCmd(object):

    '''Store help information for command line software.

    TODO: Pull from a parser

    Arguments expected in (flag, input, description) format, where a `None`
    in index 0 or 1 denotes a positional or flag only argument, respectively.
    '''

    TEMPLATE = dedent('''
        NAME{sep}{name}

        SYNOPSIS{sep}{synopsis}

        DESCRIPTION{sep}{description}

        POSITIONAL ARGUMENTS{sep}{pos}

        KEYWORD ARGUMENTS{sep}{kw}

        FLAGS{sep}{flag}

        EXAMPLES{sep}{examples}
    ''')

    def __init__(self,
                 name=None, synopsis=None, description=None,
                 arguments=None, examples=None,
                 ):
        self.dict = {
            'name':  name,
            'synopsis':  synopsis,
            'description':  description,
            'examples':  examples,
            'pos': None,
            'flag': None,
            'kw': None,
        }
        try:
            self.dict.update({
                'pos':  [a for a in arguments if a[0] is None],
                'flag':  [a for a in arguments if a[1] is None],
                'kw':  [
                    a for a in arguments
                    if a[0] is not None and a[1] is not None
                ],
            })
        except TypeError:
            pass  # most likely "'NoneType' object is not iterable"

    def __getattr__(self, attr):
        '''Check self.dict for the desired attribute and return it'''
        try:
            return self.__dict__['dict'][attr]
        except KeyError:
            raise AttributeError("'{}' object has no attribute '{}'".format(
                self.__class__.__name__, attr))

    def __str__(self):
        return self._make_help()

    def _make_help(self):
        help_str_dict = {}
        help_str_dict.update(self.dict)

        # format arguments
        sep = '\n\t'
        fmt = ['{:<15}{:<10}{}', '{:<25}{}']

        for key in ['pos', 'flag']:
            try:
                line_list = [
                    fmt[1].format(*[v for v in val if v is not None])
                    for val in help_str_dict[key]
                ]
                help_str_dict[key] = sep.join(line_list)
            except TypeError:
                pass  # no pos or flag found

        for key in ['kw']:
            try:
                line_list = [
                    fmt[0].format(*val)
                    for val in help_str_dict[key]
                ]
                help_str_dict[key] = sep.join(line_list)
            except TypeError:
                pass  # no pos or flag found
        help_str_dict['sep'] = sep
        return self.TEMPLATE.format(**help_str_dict)
