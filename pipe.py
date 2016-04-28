'''Run a genomics pipeline

Details here.
'''

# TODO(sjbush): Make subparsers.

import os.path

from libpipe import argp
from textwrap import dedent

from libpipe.pipe import align
from libpipe.util import path

import logging
log = logging.getLogger(__name__)


#
#   Setup Functions
#

def setup_logger():
    import libpipe.util.logging as pipe_log
    pipe_log.config()


def setup_parser(parser=None):

    if parser is None:
        parser = argp.core.parser()
    parser = argp.io.input_parser(parser, accept_dirs=True, accept_files=True)
    parser = argp.io.output_parser(parser, accept_dirs=False)
    parser = argp.pipe.pipe_parser(parser)

    return parser


#
#   Info Functions
#

def summarize_args(args):

    value_args = ['project', 'root', 'data']
    list_args = ['filter_list', 'file_list', 'genome_list']
    args_dict = {k: None for k in (value_args + list_args)}
    args_dict.update({k: getattr(args, k) for k in value_args})
    join_str = '\n' + ' ' * 8

    for arg in list_args:
        try:
            arg_str = join_str.join(getattr(args, arg))
        except TypeError:
            arg_str = None
        finally:
            args_dict[arg] = arg_str

    try:
        summary_list = join_str.join(
            '{}: {}'.format(k, v)
            for k, v in args.summary.items()
        )
    except AttributeError:
        summary_list = None
    finally:
        args_dict['summary'] = summary_list
    summary = dedent('''
        -----------------------------------------------------------------------
        Pipe args:
            Project: {project}
            Root dir: {root}
            Data dir: {data}
            Genome(s):
                {genome_list}
            Filter(s):
                {filter_list}
            Input (individual):
                {file_list}
            Input (summary):
                {summary}
        -----------------------------------------------------------------------
    ''').format(**args_dict)
    return(summary)


#
#   Pipe execution
#

def run_pipes(file_dict, genome_list, data=None, odir=None):

    for genome in genome_list:
        log.info('Running genome: {}'.format(genome))
        _genome = [genome]
        for name, file_list in sorted(file_dict.items()):
            log.info('   ...sample: {}'.format(name))

            try:
                _file_list = [os.path.join(data, f) for f in file_list]
                _input = _genome + _file_list
            except TypeError:
                _input = _genome + file_list

            # add name as a directory and ensure exists
            if odir:
                output_dir = os.path.join(odir, name)
                path.makedirs(output_dir)
            else:
                output_dir = None

            pipe = align.AlignPipe(
                input=_input, job_name=name, odir=output_dir)
            pipe.write()  # '~/dev/tempus/data/test_script.pbs')

    return


#
#   Main
#

def main(args):
    '''Execute pipeline using given args'''

    log.info(summarize_args(args))

    if args.summary is not None:
        # create a file_dict with {name: file} format
        file_dict = args.summary
    else:
        # convert file_list to to dict (akin to args.summary)
        # NOTE: will NOT handle pe in current form
        file_dict = {
            os.path.basename(
                os.path.splitext(f)[0]
            ): [f] for f in args.file_list
        }

    # assume the output directory based on given root and project
    try:
        odir = os.path.join(args.root, args.project, 'samples')
    except (AttributeError, TypeError):
        # arg not set or set to None
        odir = None  # output will go directly to input folder

    run_pipes(file_dict, args.genome_list, data=args.data, odir=odir)

    pass


if __name__ == '__main__':

    setup_logger()
    parser = setup_parser()
    args = parser.parse_args()
    args.find_files(args)
    args.build_args(args)
    main(args)
