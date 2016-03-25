'''A series of functions for setting up an argparse parser

Functions for setting up and adding to a single parser.

'''

import argparse

# Should we only use a single parser?
__parser = None


def parser():
    global __parser
    if not __parser:
        __parser = argparse.ArgumentParser()
    return __parser
