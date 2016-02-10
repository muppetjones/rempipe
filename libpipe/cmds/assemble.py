
import os.path
import re

from remsci.lib.decorators import file_or_handle

import libpipe.scripts
from libpipe.cmds.base import BaseCmd, BashCmd, CmdAttributes


import logging
log = logging.getLogger(__name__)


class VelvetkCmd(BaseCmd):

    '''Estimate the optimal k-mer for velvet assembly

    Command usage:
        velvetk.pl [--size X | --genome F] [options] <readfile ...>

    Script from the Victorian Bioinformatics Consortium
    bioinformatics.net.au/software.velvetk.shtml

    This script lives within this package @ libpipe.scripts.

    Arguments:

    '''

    genome_sizes = {}  # cache genome sizes
    attr = CmdAttributes(
        name='velvetk',
        synopsis='Script to calculate optimal k-mer value for velvet assembly',
        invoke_str=os.path.join(
            os.path.dirname(libpipe.scripts.__file__),
            'velvetk.pl'),

        arguments=[
            (None, 'FILE', 'The input fast[aq] file.'),
            ('--size', 'CHAR', 'The estimated genome size, e.g., 4.4M'),
            ('--genome', 'FASTA', 'A genome sequence file.'),
            ('--best', None, 'Just print the best k-mer to stdout.'),
        ],
        defaults={},

        req_args=1,
        req_kwargs=[['--size', '--genome'], ],
        req_type=[
            [('--genome', ), ('.fa', '.fna')],
            [(0, 1), ('.fa', '.fna', '.fq', '.fastq')],
        ],
    )

    def __init__(self, *args, best=True, **kwargs):
        super().__init__(*args, **kwargs)

        # If we're wrapping, then we're only interested in the BEST k-mer
        # so MAKE SURE the '--best' flag is set.
        # Likewise, if we only want the best k-mer, close stderr
        if 'wrap' in kwargs:
            best = True
        if best:
            self.flags.append('--best')
            self.redirect = '2>&-'

    def output(self):
        # the ONLY args are the sequence files, which we should pass through
        return self.args

    def _prepreq(self):

        # Prefer '--size' over genome
        # We assume we'll be using the same genome multiple times,
        # so calculating and caching the genome length should significantly
        # speed up the actual runtime of the script.
        if '--genome' in self.kwargs:
            self.kwargs['--size'] = self._calc_genome_len(
                self.kwargs['--genome'])
            del self.kwargs['--genome']

    def _calc_genome_len(self, genome):

        # assume the genome does NOT have an extension
        genome_name = os.path.basename(genome)
        try:
            # get cached genome length
            genome_len = VelvetkCmd.genome_sizes[genome_name]

        except KeyError:
            # calculate and cache genome length
            # check all possible extensions
            try:
                seq_file = [
                    genome + extn for extn in ['.fa', '.fna']
                    if os.path.isfile(genome + extn)
                ][0]
            except IndexError:
                raise self.FileTypeError('missing')

            # count all base pairs & cache
            # NOTE: no content checking to ensure actual base pairs
            with open(seq_file, 'r') as fh:
                genome_len = sum(
                    len(line.rstrip()) for line in fh
                    if not line.startswith('>')
                )
            VelvetkCmd.genome_sizes[genome_name] = genome_len

        finally:
            return genome_len


class VelvethCmd(BaseCmd):

    '''Setup hashing for Velvet assembly

    Command usage:
        velveth directory hash_length \
            {[-file_format][-read_type][-separate|-interleaved] \
             filename1 [filename2 ...]} {...} [options]

    NOTE on run order:
        velvetk (or similar) should be run FIRST to estimate k-mer length.
        velvetg MUST be run AFTER.

    NOTE: velveth and velvetg SHOULD be run together but share
        data via a common directory. Because this directory will be
        created at runtime, preventing type checking, velveth.output will
        return [<outdir>/contigs.fa], which velvetg should match as *.fa
        to identify the output directory.

    NOTE: velveth requires the output directory and k-mer hash length
        to come BEFORE all other parameters. To accomodate this unorthodox
        requirement, plese use the 'output' and 'k' arguments during
        initialization.

    NOTE: velveth has a max k-mer length value that is defined during
        compliation, with a default of 31.

    Arguments:
        output: The output directory. Default: <Input file directory>/k<val>
        k: The hash length. Should be an odd integer or bash varaible
    '''

    attr = CmdAttributes(
        name='velveth',
        synopsis='Hashing for genome assembly',
        invoke_str='velveth',

        arguments=[
            (0, 'DIR', 'The output directory'),
            (1, 'INT', 'Hash length. May also pass in a BASH variable'),
            (None, 'FILE', 'The input fast[aq] file(s).'),
        ],
        defaults={},

        req_args=3,
        req_kwargs=[],
        req_type=[
            [('--genome', ), ('.fa', '.fna')],
            [(2, 3), ('.fa', '.fna', '.fasta', '.fq', '.fastq')],
        ],
        n_priority_args=2,  # output directory and k hash length FIRST!
        allow_bash_var=True,

    )

    def __init__(self, *args, output=None, k=None, **kwargs):
        super().__init__(*args, **kwargs)

        self.args = [output, k] + self.args

    def output(self):
        # velveth and velvetg use a common directory with static file names
        # return the output directory only
        return [self.args[0], ]

    def _prepreq(self):
        '''Prepare cmd for requirements check

        1. Ensure the output directory is set correctly
        '''

        # set directory if not given
        # NOTE: we should have the input files at this point
        # NOTE: If the m,M,s format is given for k, replace commas with '_'
        if not self.args[0]:
            self.args[0] = os.path.join(
                os.path.dirname(self.args[2]),
                'k' + str(self.args[1].replace(',', '_')),
            )

    def _additional_requirements(self):
        '''Ensure k is set'''

        if not self.args[1]:
            raise self.PositionalArgError('Missing required k-mer value')

    def _prepcmd(self):

        # assume single-end if only one file
        # assume paired-end if two files
        input_seq = self.args[2:]
        read_len = self._estimate_read_len(input_seq[0])
        flag = '-long' if read_len > 200 else '-short'
        if len(input_seq) > 1:
            flag = flag + 'Paired -separate'
        flag = flag.split()

        # add flag for type!
        extn = os.path.splitext(input_seq[0])[1]
        if extn in ['.fa', '.fasta']:
            flag.append('-fasta')
        elif extn in ['.fq', '.fastq']:
            flag.append('-fastq')
        else:
            flag.append('-fmtAuto')

        self.flags.extend(flag)

    @file_or_handle(mode='r')
    def _estimate_read_len(self, fh, n=10):
        '''Estimate the read length

        Arguments:
            seqfile: The file to estimate read length
            n: The number of sequences to use to estimate
        Return:
            The rounded average of the first n sequences
        '''

        rx_alphanum = re.compile(r'^\w+$')

        line_len = []
        for line in fh:
            line = line.rstrip()
            if line[0] in ['@', '>', '+']:
                continue  # fasta or fastq header
            if rx_alphanum.search(line) is None:
                continue  # not a sequence (propably a quality score)
            if len(line_len) > n:
                break

            line_len.append(len(line))

        return round(sum(line_len) / len(line_len))


class VelvetgCmd(BaseCmd):

    '''Construct a de Bruijn graph, etc., to assemble a bacterial genome

    TODO: More advanced handling of type matching to allow for directories
        (should check for at least partial path existence).

    HACK: Currently connects velveth.output directory to self.input directory
        via an empty extension match. See above TODO.

    '''

    attr = CmdAttributes(
        name='velvetg',
        synopsis='De Bruijn graph assembly for de novo genome assembly',
        invoke_str='velvetg',

        arguments=[
            (0, 'DIR', 'The output directory'),
            ('-min_contig_lgth', 'INT',
             'Minimum contig length exported to contigs.fa. Default: 500.'),
            ('-exp_cov', 'FLOAT|auto',
             'Expected coverage of unique regions (or infer). Default: auto'),
            ('-cov_cutoff', 'FLOAT|auto',
             'Coverage threshold for removing nodes. Default: auto'),
            ('-ins_length', 'INT', 'Expected insert length'),
            ('-ins_length_sd', 'INT', 'Expected insert length standard dev'),
            ('-unused_reads', 'yes|no',
             'export unused reads to UnusedReads.fa. Default: yes'),
            ('-clean', 'yes|no',
             'Remove some intermediary files. Default: yes'),
            ('-very_clean', 'yes|no',
             'Remove all intermediary files. Default: no'),
        ],
        defaults={
            '-min_contig_lgth': 500,
            '-exp_cov': 'auto',
            '-cov_cutoff': 'auto',
            '-unused_reads': 'yes',
            '-clean': 'yes',
            '-very_clean': 'no',
        },

        n_priority_args=1,
        req_args=1,
        req_kwargs=[],
        req_type=[
            [(0, ), ('', )],  # HACK!
        ],
        allow_bash_var=True,
    )

    def output(self):
        # velveth and velvetg use a common directory with static file names
        # return the output directory only
        output_files = [
            os.path.join(self.args[0], f)
            for f in ('contigs.fa', 'stats.txt')
        ]
        return [self.args[0]] + output_files


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
            [('-q', ), ('.fa', '.fasta', '.fna'), ],
            [('-r', ), ('.fa', '.fasta', '.fna'), ],
        ],
        allow_bash_var=True,
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

        log.debug(self.kwargs)

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


class ContigSummaryCmd(BashCmd):

    '''Write contig size summary

    Command:
        grep ">" "$f" |awk -F '_' '{print $4}' |sort -rn > "${contigsize}"

    NOTE: BASH one-liner. Must not do a full cleaning on the command, lest
        we mangle it beyond recognition.

    Arguments:
        filename: Contig file from velvet
        outfile: The file to write to. Default: <dirname>/contig_sizes.txt
            (on init only)
    Output:
        1) A file with all contig names and sizes (contig_sizes.txt)
        2) [fall_through] The input fasta file.
    '''

    attr = CmdAttributes(
        name='contig_summary',
        synopsis='De Bruijn graph assembly for de novo genome assembly',
        invoke_str='grep ">"',

        arguments=[
            (0, 'FILE', 'The contig file'),
        ],

        req_args=1,
        req_type=[
            [(0, ), ('.fa', '.fasta', )],
        ],
        allow_bash_var=True,
    )

    def __init__(self, *args, outfile=None, **kwargs):
        super().__init__(*args, **kwargs)

        if outfile:
            self.redirect = ('>', outfile)

        if len(self.args) > 1:
            raise self.PositionalArgError('unknown')

        # second part of one liner
        # -- get the fourth column and sort
        self.args.append("|awk -F '_' '{print $4}' |sort -rn")

    def _prepreq(self):

        # ensure file is first, rest of cmd is second
        if self.args[0].startswith('|awk'):
            self.args.reverse()
            log.debug(self.args)

    def _prepcmd(self):

        # ensure double quotes around input file
        # BREAKS EVERYTHING!
        # if not self.args[0].startswith('"'):
        #     self.args[0] = '"{}"'.format(self.args[0])

        if not self.redirect:
            outfile = '_'.join([
                os.path.splitext(self.args[0])[0],
                'contig_sizes.txt',
            ])
            self.redirect = ('>', outfile)

    def output(self):
        return [self.args[0], self.redirect[1]]
