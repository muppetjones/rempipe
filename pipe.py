'''Run a genomics pipeline

Details here.
'''

# TODO(sjbush): Make subparsers.

import os.path

from libpipe import argp
from textwrap import dedent

from libpipe.pipe import align

import logging
log = logging.getLogger(__name__)


#
#   Setup Functions
#

def setup_logger():
    import libpipe.util.logging as pipe_log
    pipe_log.config()
    return log


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

    value_args = ['project', 'root', 'data', 'genome']
    list_args = ['filter_list', 'file_list']
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
            Genome: {genome}
            Filters:
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

def run_pipes(file_dict):

    for name, file_list in sorted(file_dict.items()):
        log.debug(name)
        log.debug(file_list)
        _input = file_list
        pipe = align.AlignPipe(input=file_list)
        pipe.write('~/dev/tempus/data/test_script.pbs')


#
#   Main
#

def main(args):

    log.info(summarize_args(args))

    if args.summary is not None:
        file_dict = args.summary
    else:
        # convert file_list to to dict (akin to args.summary)
        # NOTE: will NOT handle pe in current form
        file_dict = {
            os.path.basename(
                os.path.splitext(f)[0]
            ): [f] for f in args.file_list
        }
    run_pipes(file_dict)

    pass


if __name__ == '__main__':

    log = setup_logger()
    parser = setup_parser()
    args = parser.parse_args()
    args.find_files(args)
    main(args)
