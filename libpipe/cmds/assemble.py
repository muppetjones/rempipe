
import os.path

import libpipe.scripts
from libpipe.cmds.base import BaseCmd, CmdAttributes


class VelvetkCmd(BaseCmd):

    '''Estimate the optimal k-mer for velvet assembly

    Script from the Victorian Bioinformatics Consortium
    bioinformatics.net.au/software.velvetk.shtml

    This script lives within this package @ libpipe.scripts.velvetk.

    Arguments:

    '''

    attr = CmdAttributes(
        name='velvetk',
        synopsis='Script to calculate optimal k-mer value for velvet assembly',
        invoke_str=os.path.join(
            os.path.dirname(libpipe.scripts.__file__),
            'velvetk.pl'),

        arguments=[
            (0, 'FILE', 'The input fast[aq] file.'),
            ('--size', 'CHAR', 'The estimated genome size, e.g., 4.4M'),
            ('--genome', 'FASTA', 'A genome sequence file.'),
            ('--best', None, 'Just print the best k-mer to stdout.'),
        ],
        defaults={},

        req_args=1,
        req_kwargs=[['--size', '--genome'], ],
        req_type=[
            [('--genome', ), ('.fa', '.fna')],
            [(0, ), ('.fa', '.fna', )],
        ],
    )


class VelvethCmd(BaseCmd):
    attr = CmdAttributes(
        name='velveth',
        synopsis='Script to calculate optimal k-mer value for velvet assembly',
        invoke_str=os.path.join(
            os.path.dirname(libpipe.scripts.__file__),
            'velvetk.pl'),

        arguments=[
            (0, 'DIR', 'The output directory'),
            (1, ''),
            (0, 'FILE', 'The input fast[aq] file.'),
            ('--size', 'CHAR', 'The estimated genome size, e.g., 4.4M'),
            ('--genome', 'FASTA', 'A genome sequence file.'),
            ('--best', None, 'Just print the best k-mer to stdout.'),
        ],
        defaults={},

        req_args=3,
        req_kwargs=[['--size', '--genome'], ],
        req_type=[
            [('--genome', ), ('.fa', '.fna')],
            [(0, ), ('.fa', '.fna', )],
        ],

    )


class VelvetgCmd(BaseCmd):
    pass
