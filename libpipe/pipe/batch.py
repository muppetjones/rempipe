

from libpipe.pipe import base


import logging
log = logging.getLogger(__name__)


class BatchPipe(base.Pipe):

    '''Pipe for running additional steps based on multiple individual pipes

    BatchPipe is meant to run commands that require the output from
    multiple, independent sources, such as running DESeq2 on several RNA-Seq
    samples.

    BatchPipe behaves like most other pipes with a few notable exceptions:
    *   BatchPipe should _not_ be nested in other pipes. There are possible
        exceptions, but for now, the link method will raise TypeError if
        called.
    *   BatchPipe accepts an additional output directory. The standard
        output_dir argument is passed to the sub-pipes and commands, but
        output for the batched commands may be redirected to a different
        folder.
    *   BatchPipe includes an additional method called 'batch' that behaves
        exactly like 'add', except that commands and pipes added through
        'add' will be run for each input and those added via 'batch' will
        be run _after_ each command has completed.
    *   Input must be given as a dict with a {jobname: filename(s)} format.
    '''

    def __init__(self, **kwargs):

        try:
            _input = kwargs['input']
            del kwargs['input']
        except KeyError:
            _input = {}
        finally:
            if not isinstance(_input, dict):
                msg = 'BatchPipe input must be a dict'
                raise ValueError(msg)

    #
    #   CmdInterface
    #

    def cmd(self):
        pass

    def link(self, *args, **kwargs):
        msg = 'BatchPipe objects are not linkable.'
        raise TypeError(msg)

    def output(self):
        pass
