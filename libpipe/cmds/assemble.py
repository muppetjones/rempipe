
import os.path
import re

from remsci.lib.decorators import file_or_handle

import libpipe.scripts
from libpipe.cmds.base import BaseCmd, CmdAttributes


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
            [(0, 1), ('.fa', '.fna', )],
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
            ('--size', 'CHAR', 'The estimated genome size, e.g., 4.4M'),
            ('--genome', 'FASTA', 'A genome sequence file.'),
            ('--best', None, 'Just print the best k-mer to stdout.'),
            (None, 'FILE', 'The input fast[aq] file(s).'),
        ],
        defaults={},

        req_args=3,
        req_kwargs=[['--size', '--genome'], ],
        req_type=[
            [('--genome', ), ('.fa', '.fna')],
            [(2, 3), ('.fa', '.fna', )],
        ],
        n_priority_args=2,  # output directory and k hash length FIRST!
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

        self.flags.extend(flag.split())

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
            if rx_alphanum(line) is None:
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
    )

    def output(self):
        # velveth and velvetg use a common directory with static file names
        # return the output directory only
        output_files = [
            os.path.join(self.args[0], f)
            for f in ('contigs.fa', 'stats.txt')
        ]
        return [self.args[0]] + output_files
