
import os.path
from libpipe.cmds.base import BaseCmd

import logging
log = logging.getLogger(__name__)


class FastqcCmd(BaseCmd):

    '''FastQC command

    Command Usage:
        fastqc [-o output dir] [--(no)extract] [-f fastq|bam|sam]
               [-c contaminant file] seqfile1 .. seqfileN
    '''

    NAME = 'fastqc'
    INVOKE_STR = 'fastqc'

    DEFAULTS = {}

    attributes = {
        '-o': "output directory",
    }

    REQ_KWARGS = []
    REQ_ARGS = 1

    def output(self):

        files = []
        for f in self.args:
            file_base = os.path.splitext(f)[0] + '_fastqc'
            files.extend([
                file_base + '.html',
                file_base + '.zip',
            ])

        # update the path if 'o' was given
        try:
            files = [
                os.path.join(self.kwargs['-o'], os.path.basename(f))
                for f in files
            ]
        except KeyError:
            pass

        return files

    def _prepcmd(self):
        # ensure the '-o' option is given
        # -- use the common prefix of the two input files
        # -- or the basename of the the first file otherwise
        if '-o' not in self.kwargs:
            prefix = os.path.commonprefix(self.args)
            if not prefix or prefix == self.args[0]:
                prefix = os.path.splitext(self.args[0])[0]
            self.kwargs['-o'] = prefix


class SamtoolsSortCmd(BaseCmd):

    '''Samtools Sort

    Creates command to generate a sorted BAM file from a SAM or BAM file.

    Command Usage:
        samtools sort [options...] <in>

    Command Usage (legacy):
        samtools sort [options...] <in> <out.prefix>
    '''

    NAME = 'samtools_sort'
    INVOKE_STR = 'samtools sort'

    ARGUMENTS = [
        ('-o', 'FILE', 'File to write final output to'),
        ('-O', 'FORMAT', 'Write output as FORMAT (sam/bam/cram)'),
        ('-T', 'PREFIX', 'Write temporary files to PREFIX.nnn.bam'),
    ]
    DEFAULTS = {}

    REQ_KWARGS = []
    REQ_ARGS = 1
    REQ_TYPE = [
        [(0, ), ('.bam', '.sam')]
    ]

    def _prepcmd(self):
        # if we were only given a single file, make the prefix
        if '-o' not in self.kwargs:
            if self.args[0].endswith('sam'):
                bam_file = os.path.splitext(self.args[0])[0] + '.bam'
            else:
                bam_file = os.path.splitext(self.args[0])[0] + '.s.bam'
            self.kwargs['-o'] = bam_file

        if '-T' not in self.kwargs:
            t_kw = os.path.splitext(self.kwargs['-o'])[0] + '.tmp'
            self.kwargs['-T'] = t_kw

    @property
    def output(self):
        return [self.kwargs['-o'], ]


class SamtoolsIndexCmd(BaseCmd):

    '''Samtools Index

    Command Usage:
        samtools index [-bc] [-m INT] <in.bam>
    '''

    NAME = 'samtools_index'
    INVOKE_STR = 'samtools index'

    ARGUMENTS = [
        (None, 'FILE', 'BAM file to index'),
        ('-b', None, 'Generate BAI-format index for BAM files [default]'),
        ('-c', None, 'Generate CAI-format index for BAM files'),
        ('-m', 'INT', 'Set min interval size for CSI indices to 2^INT'),
    ]
    DEFAULTS = {}

    REQ_KWARGS = []
    REQ_ARGS = 1
    REQ_TYPE = [
        [(0, ), ('.bam')]
    ]

    def __init__(self, *args, fmt='-b', **kwargs):
        super().__init__(*args, **kwargs)

        self.flags.extend([
            fmt,
        ])

    def output(self):
        return [self.args[0], ]  # return the given bam file
