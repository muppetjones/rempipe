
# import argparse
# import os.path

from abc import ABCMeta, abstractmethod

# import remsci.scripted.base as base
# from remsci.lib.parsers import gff


# from remsci.lib.parsers.readers import __scripted__

import logging
log = logging.getLogger(__name__)


# -----------------------------------------------------------------------------
# Subparser creation

class SubparserBase(metaclass=ABCMeta):

    def setup(self):
        self.subparser.set_defaults(func=self.run)

    @abstractmethod
    def run(self, args):
        pass
