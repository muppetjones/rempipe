'''Define the interface for a command line object'''


import abc


class CmdInterface(metaclass=abc.ABCMeta):
    '''Define an interface for command-like objects

    Abstract methods:
        cmd: Return the bash-executable string of the command.
        link: Pair current 'output' with given command's 'input'.
            Should return the given command to allow chaining.
        output: Return a list of command output. Should be string paths
            in *most* cases.
    '''
    @abc.abstractmethod
    def cmd(self):
        pass

    @abc.abstractmethod
    def output(self):
        pass

    @abc.abstractmethod
    def link(self):
        pass
