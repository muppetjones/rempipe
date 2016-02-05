
import os.path
from libpipe.cmds.base import BaseCmd, CmdAttributes

import logging
log = logging.getLogger(__name__)


class FastqcCmd(BaseCmd):

    '''FastQC command

    Command Usage:
        fastqc [-o output dir] [--(no)extract] [-f fastq|bam|sam]
               [-c contaminant file] seqfile1 .. seqfileN
    '''

    attr = CmdAttributes(
        name='fastqc',
        invoke_str='fastqc',

        arguments=[
            (None, 'FILE', 'One or more FASTQ files'),
            ('-o', 'DIR', 'Output directory'),
        ],
        defaults={},

        req_kwargs=[],
        req_args=1,
        req_type=[
            [(0, 1), ('.fastq', '.fq')],
        ],
    )

    def output(self):
        '''FastQC output

        Returns:
            A list containing: (1-n) The original fastq files given,
            followed by ((n+2) * n) the html and zip file produced for
            each given fastq file.
        '''

        out_extn = ['.html', ]
        if '--extract' not in self.flags:
            out_extn.append('.zip')

        files = list(self.args)
        for f in self.args:
            file_base = os.path.splitext(f)[0] + '_fastqc'
            try:
                # update the path if '-o' was given
                file_base = os.path.join(
                    self.kwargs['-o'], os.path.basename(file_base))
            except KeyError:
                pass
            files.extend([
                file_base + extn for extn in out_extn
            ])

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

    attr = CmdAttributes(
        name='samtools_sort',
        invoke_str='samtools sort',

        arguments=[
            ('-o', 'FILE', 'File to write final output to'),
            ('-O', 'FORMAT', 'Write output as FORMAT (sam/bam/cram)'),
            ('-T', 'PREFIX', 'Write temporary files to PREFIX.nnn.bam'),
        ],
        defaults={},

        req_kwargs=[],
        req_args=1,
        req_type=[
            [(0, ), ('.bam', '.sam')]
        ],
    )

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

    def output(self):
        return [self.kwargs['-o'], ]


class SamtoolsIndexCmd(BaseCmd):

    '''Samtools Index

    Command Usage:
        samtools index [-bc] [-m INT] <in.bam>
    '''

    attr = CmdAttributes(
        name='samtools_index',
        invoke_str='samtools index',

        arguments=[
            (None, 'FILE', 'BAM file to index'),
            ('-b', None, 'Generate BAI-format index for BAM files [default]'),
            ('-c', None, 'Generate CAI-format index for BAM files'),
            ('-m', 'INT', 'Set min interval size for CSI indices to 2^INT'),
        ],
        defaults={},

        req_kwargs=[],
        req_args=1,
        req_type=[
            [(0, ), ('.bam')]
        ],
    )

    def __init__(self, *args, fmt='-b', **kwargs):
        super().__init__(*args, **kwargs)

        self.flags.extend([
            fmt,
        ])

    def output(self):
        return [self.args[0], ]  # return the given bam file


class SamtoolsViewCmd(BaseCmd):

    '''Samtools View

    Prints alignments in SAM format

    Command Usage:
        samtools view [options...] <in.bam|in.sam|in.cram> [region]

    Command Usage (legacy):
        samtools sort [options...] <in> <out.prefix>
    '''

    attr = CmdAttributes(
        name='samtools_view',
        invoke_str='samtools view',

        arguments=[
            (None, 'FILE', 'SAM|BAM|CRAM input file'),
            (None, 'REGION', 'Region(s) as 1-based RNAME[:STARTPOS[-ENDPOS]]'),
            ('-o', 'FILE', 'File to write final output to'),
            ('-f', 'INT',
             'Output alignments with INT bits set in the FLAG field'),
            ('-F', 'INT',
             'Output alignments with INT bits NOT set in the FLAG field'),
            ('-b', None, 'Output to BAM format'),
            ('-h', None, 'Include headers in output'),
        ],
        defaults={},

        req_kwargs=[],
        req_args=1,
        req_type=[
            [(0, ), ('.bam', '.sam')]
        ],
    )

    def __init__(self, *args, region=None, **kwargs):
        super().__init__(*args, **kwargs)

        self.region = region

    def _prepreq(self):
        # CRITICAL: Ensure the SAM file comes BEFORE the region(s)
        if len(self.args) > 1 and not self.args[0].endswith('am'):
            self.args = [arg for arg in self.args if arg.endswith(
                'am')] + [arg for arg in self.args if not arg.endswith('am')]

        # Append any regions given during init
        if self.region:
            self.args.append(self.region)

    def _prepcmd(self):
        if '-o' not in self.kwargs:
            # sam to bam
            if self.args[0].endswith('.sam') and '-b' in self.flags:
                ofile = self.args[0].replace('.sam', '.bam')
            # bam to sam
            elif self.args[0].endswith('.bam') and '-b' not in self.flags:
                ofile = self.args[0].replace('.bam', '.sam')
            # no type conversion (CRAM not implemented)
            else:
                ofile = self.args[0]

            # filter using -f and|or -F
            base, extn = os.path.splitext(ofile)
            for flag in ['-f', '-F']:
                try:
                    flag_str = '{}{}'.format(flag[1:], self.kwargs[flag])
                    base = '{}_{}'.format(base, flag_str)
                except KeyError:
                    pass
            ofile = base + extn

            # DO NOT OVERWRITE!! -- if it's the same name, just add SOMETHING
            if ofile == self.args[0]:
                ofile = base + '_filt' + extn

            self.kwargs['-o'] = ofile

    def output(self):
        return [self.kwargs['-o'], ]  # return the given bam file


class AbacasCmd(BaseCmd):

    '''Reorder contigs to match reference genome

    Usage:
        abacas.pl -r <reference file: single fasta> \
            -q <query sequence file: fasta> \
            -p <nucmer/promer>  [OPTIONS]

    Arguments:
        -r	reference sequence in a single fasta file
        -q	contigs in multi-fasta format
        -p	MUMmer program to use: 'nucmer' or 'promer'

    NOTE: If '-p' is not given, we will attempt to determine
        nucleotide or protein.

    NOTE: If a genome prefix is given, will attempt to identify correct
        sequence file.
    '''

    attr = CmdAttributes(
        name='abacas',
        synopsis='',
        invoke_str='',

        arguments=[
            ('-r', 'FILE', 'reference sequence (single file).'),
            ('-q', 'FILE', 'contigs (multi-fasta format).'),
            ('-p', 'nucmer|promer', 'MUMmer program to use.'),
            ('-o', 'CHAR', 'prefix for output files.'),
            ('-c', None, 'Reference sequence is circular. Default: True'),
            ('-m', None, 'Print ordered contigs to fasta file'),
            ('-b', None, 'Print contigs in bin to file'),
        ],

        req_kwargs=['-r', '-q', ],
        req_type=[
            [('-r', '-q'), ('.fa', '.fasta', '.fna'), ],
        ],
    )

    def __init__(
            self, *args,
            circular=True, multifasta=True, binfile=True, **kwargs):

        super().__init__(*args, **kwargs)

        if circular and '-c' not in self.flags:
            self.flags.append('-c')

        if multifasta and '-m' not in self.flags:
            self.flags.append('-m')

        if binfile and '-b' not in self.flags:
            self.flags.append('-b')

    def _prepreq(self):
        # ensure genome extension
        try:
            self.kwargs['-r'] = self._ensure_file_and_extension(
                self.kwargs['-r'], self.attr.get_types('-r'))
        except self.FileTypeError:
            # No extn given OR unable to find a file with an expected extn
            raise

    def _prepcmd(self):

        if '-o' not in self.kwargs:
            # assume the contigs are in <sample>/kXX/contigs.fa
            filepath = os.path.dirname(self.kwargs['-q'])
            filename = os.path.basename(os.path.dirname(filepath))
            genome = os.path.splitext(os.path.basename(self.kwargs['-r']))[0]

            self.kwargs['-o'] = os.path.join(
                filepath, '{}_{}'.format(filename, genome))

    def output(self):
        extns = ['.fasta', '.bin', '.tab']
        return [self.kwargs['-o'] + extn for extn in extns]
