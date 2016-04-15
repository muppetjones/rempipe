
import os.path

from libpipe.cmd import align  # HOTFIX: remove ASAP
from libpipe.cmd import attr
from libpipe.cmd import base
from libpipe.type import seq

import logging
log = logging.getLogger(__name__)


class SamtoolsIndexCmd(base.CmdBase):

    '''Samtools Index

    Command Usage:
        samtools index [-bc] [-m INT] <in.bam>

    TODO(sjbush): Better handling for [(-c, -m), -b]
    TODO(sjbush): Enable fall through of genome index.

    Attributes:
        fmt: The string flag used to denote the index format. Default: "-b".

    Known Issues:
        * If 'None' is used instead of '0' in attr.args, the required type
            does not match up input to the arg.
    '''

    attr = attr.CmdAttributes(
        name='samtools_index',
        invoke='samtools index',

        args=[
            (0, 'FILE', 'BAM file to index'),
            ('-b', None, 'Generate BAI-format index for BAM files [default]'),
            ('-c', None, 'Generate CAI-format index for BAM files'),
            ('-m', 'INT', 'Set min interval size for CSI indices to 2^INT'),
        ],
        defaults={},

        req_kwargs=[],
        req_args=1,
        req_types=[
            [(0, ), (seq.BamType, )]
        ],
    )

    def __init__(self, *args, fmt=None, **kwargs):
        super().__init__(*args, **kwargs)

        if not fmt and not self.flags:
            self.flags.append('-b')

    def output(self):
        '''Return given BAM file (unchanged) and genome index (if given)'''
        index_list = [
            i for i in self.input()
            if isinstance(i, align.Hisat2Index)  # index.IndexType)  # HOTFIX
        ]
        return [self.args[0]] + index_list


class SamtoolsSortCmd(base.CmdBase):

    '''Samtools Sort

    Creates command to generate a sorted BAM file from a SAM or BAM file.
    By default, the output is sent to a bam file and tagged with '.s' to
    help avoid name collision:
        e.g., data/sample.sam --> data/sample.s.bam

    TODO(sjbush): Enable fall through of genome index.

    Command Usage (v1.3+):
        samtools sort [options...] <in>

    Command Usage (v1.2 and older):
        NOT IMLEMENTED!
        samtools sort [options...] <in> <out.prefix>
    '''

    attr = attr.CmdAttributes(
        name='samtools_sort',
        invoke='samtools sort',

        args=[
            (None, 'FILE', 'Input SAM or BAM file'),
            ('-o', 'FILE', 'File to write final output to'),
            ('-O', 'FORMAT', 'Write output as FORMAT (sam/bam/cram)'),
            ('-T', 'PREFIX', 'Write temporary files to PREFIX.nnn.bam'),
        ],
        defaults={},

        req_kwargs=[],
        req_args=1,
        req_types=[
            [(0, ), ('.bam', '.sam')]
        ],
    )

    #
    #   Cmd Interface
    #

    def output(self):
        '''Return sorted BAM file and genome index (if given)'''
        index_list = [
            i for i in self.input()
            if isinstance(i, align.Hisat2Index)  # index.IndexType)  # HOTFIX
        ]
        return [self.kwargs['-o']] + index_list

    #
    #   Overrides
    #

    def _pre_cmd(self):

        # ensure the output is set
        if '-o' not in self.kwargs:
            self.kwargs['-o'] = '{}.s.bam'.format(
                os.path.splitext(self.args[0])[0])

        # ensure that -T is set if -o is (required)
        if '-o' in self.kwargs and '-T' not in self.kwargs:
            self.kwargs['-T'] = '{}.tmp'.format(
                os.path.splitext(self.kwargs['-o'])[0])
