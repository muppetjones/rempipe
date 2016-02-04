import os

from libpipe.pipes.base import PresetPipe
from libpipe.pipes.qc import TrimPipe
from libpipe.cmds import (
    Hisat2Cmd,
    SamtoolsSortCmd, SamtoolsIndexCmd,
    HtseqCountCmd,
)

from libpipe.cmds.assemble import (
    VelvetkCmd, VelvethCmd, VelvetgCmd,
)


import logging
log = logging.getLogger(__name__)


class AssemblePipe(PresetPipe):

    '''Assemble a bacterial genome

    NOTE: This pipe is intended for use as a subpipe.

    Velvet Pipeline:
        a. Velvetk
        b. Velveth
        c. Velvetg
    '''

    REQ_PARAM = ['input_list', 'genome']

    def _setup(self, *pipe_args, input_list=[], genome='', **pipe_kwargs):

        # 1. velvetk
        args = input_list[:]
        kwargs = {'--genome': genome}
        velvetk = VelvetkCmd(*args, wrap='k', **kwargs)

        # 2. velveth
        args = []
        kwargs = {'k': '${k}'}
        velveth = VelvethCmd(*args, **kwargs)

        # 3. velvetg
        args = []
        kwargs = {}
        velvetg = VelvetgCmd(*args, **kwargs)

        # add 'em to the pipe
        self.add(velvetk, velveth, velvetg)


class WgsPipe(PresetPipe):

    '''Assemble a bacterial genome, from sequencer to annotation

    Pipeline:
        1. TrimPipe
            a. FastQC (reads)
            b. Skewer (reads)
            c. FastQC (skewered reads)
        2. (optional) Filter reads through reference genome (usually human)
            a. Align to filter genome (will use unaligned output)
        3. Determine insertion length (skip)
            a. Remove unmapped reads
            b. Calculate insert length
        4. Run through Velvet
            a. estimate k (velvetk.pl)
            b. velveth
            c. velvetg
        5. Assembly statistics


    '''

    REQ_PARAM = ['input_list', 'odir']

    def _setup(self, *pipe_args, **pipe_kwargs):
        # NOTE: TrimPipe uses 'input_list' to set FastQC input.

        # check for filter genome input

        # 1. TrimPipe
        trim = TrimPipe(*pipe_args, **pipe_kwargs)

        # 2. (optional) Filter through reference genome
        # > input: link from previous (fastq)
        # > output: sam files to given directory (+ the unaligned fastq)
        if 'genome' in pipe_kwargs:
            args = []
            kwargs = {'-x': pipe_kwargs['genome'], '-p': os.cpu_count()}
            if kwargs['-p'] is None:
                del kwargs['-p']
            else:
                kwargs['-p'] = kwargs['-p'] - 1
            hisat = Hisat2Cmd(*args, **kwargs)

        # 3. Determine Insertion length
        # ** SKIP **

        # 4a. velvetk -- find hash length!
        args = []
        kwargs = {}

        # CRITICAL!! This will be a problem
        # we need to get the value returned from here and use it in
        # velveth -- as a value, not a file!!
        # two options:
        #   1. Add a special _prep_cmd in velveth to retrieve it
        #       -- BUT the data may not exist to calculate at that point
        #   2. Add a special *pre-cmd* cmd
        #       -- i.e., velveth.cmd will return something like:
        #           k="$(velvetk.pl ....)
        #           velveth $k
        #   3. Allow for types to be used in req args
        #       -- cmd.output will be able to return those types
        #       -- probably still need better linking...somehow
        #   NO! Must be 2, as a sub command or something...3 same issue as 1

        # 4b. velveth
        args = []
        kwargs = {}
        velveth = VelvethCmd(*args, **kwargs)

        # 4c. velvetg
        args = []
        kwargs = {}
        velvetg = VelvetgCmd(*args, **kwargs)

        # 5. Assembly statistics

        # Add 'em all up
        try:
            self.add(trim, hisat, velveth, velvetg)
        except UnboundLocalError:
            self.add(trim, velveth, velvetg)
