
import os.path

from libpipe.cmd.attr import CmdAttributes
from libpipe.cmd.base import CmdBase


class SamtoolsSortCmd(CmdBase):

    '''Samtools Sort

    Creates command to generate a sorted BAM file from a SAM or BAM file.
    By default, the output is sent to a bam file and tagged with '.s' to
    help avoid name collision:
        e.g., data/sample.sam --> data/sample.s.bam

    Command Usage (v1.3+):
        samtools sort [options...] <in>

    Command Usage (v1.2 and older):
        NOT IMLEMENTED!
        samtools sort [options...] <in> <out.prefix>
    '''

    attr = CmdAttributes(
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
        return [self.kwargs['-o']]

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
