'''Run a genomics pipeline

Details here.
'''

# TODO(sjbush): Make subparsers.

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


def setup_parser():

    parser = argp.core.parser()
    _ = argp.io.input_parser(parser, accept_dirs=True, accept_files=True)
    _ = argp.io.output_parser(parser, accept_dirs=False)
    _ = argp.pipe.pipe_parser(parser)
    return _


#
#   Info Functions
#

def summarize_args(args):

    value_args = ['project', 'root', 'data', 'genome', 'summary']
    list_args = ['filter_list', 'file_list']
    args_dict = {k: getattr(args, k) for k in value_args}
    args_dict.update({
        k: ('\n{}'.format(' ' * 8).join(getattr(args, k))
            if getattr(args, k) is not None else None)
        for k in list_args
    })
    log.debug(args_dict)
    summary = dedent('''
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
    ''').format(**args_dict)
    return(summary)


#
#
#


#
#   Main
#

def main(args):

    log.info(summarize_args(args))

    #

    pipe = align.AlignPipe()
    pass


if __name__ == '__main__':

    log = setup_logger()
    parser = setup_parser()
    args = parser.parse_args()
    args.find_files(args)
    main(args)
