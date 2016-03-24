
from libpipe.cmd.attr import CmdAttributes
from libpipe.cmd.base import CmdBase


class SamtoolsSortCmd(CmdBase):

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

    def output(self):
        return []
