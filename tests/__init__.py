import sys
from os.path import dirname, abspath, join

print("rempipe test init -- add 'remsci' to path")
sys.path.insert(
    0, join(dirname(dirname(dirname(abspath(__file__)))), 'remsci'))

from remsci.lib.utility import customLogging
customLogging.config()
