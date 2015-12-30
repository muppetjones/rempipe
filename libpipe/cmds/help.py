
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
        NAME
            {name}

        SYNOPSIS
            {synopsis}

        DESCRIPTION
            {description}

        POSITIONAL ARGUMENTS
            {pos}

        KEYWORD ARGUMENTS
            {kw}

        FLAGS
            {flag}

        EXAMPLES
            {examples}
    ''')

    def __init__(self,
                 name=None, synopsis=None, description=None,
                 args=None, examples=None,
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
                'pos':  [a for a in args if a[0] is None],
                'flag':  [a for a in args if a[1] is None],
                'kw':  [
                    a for a in args
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
        return self.TEMPLATE.format(**self.dict)
