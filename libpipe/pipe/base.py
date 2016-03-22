

from libpipe import cmd


class PipeBase(object):

    '''Pipe object for storing and handling a series of command-line objs

    PipeBase stores a series of command line objects and can execute those
    commands a number of fashions
        * Via a BASH script
        * Directly through python subprocesses
        * Through a resource manager system, such as Torque

    TODO(sjbush): Add Torque and GridEngine handling (both individually
        and through a pbs script).
    TODO(sjbush): Add direct execution via subprocesses (via cmd.run)

    Attributes:


    '''

    def __init__(self):

        cmd.CmdDummy()
