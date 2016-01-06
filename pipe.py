#!/usr/bin/env python

import os.path
import subprocess
import sys
import time
sys.path.append("../remsci/")

import remsci.scripted.base as base
from remsci.lib.utility import path
from libpipe.pipes.genomics import NestedGenomicsPipe
from libpipe.parsers.fastq import FastqScripted
# from libpipe.cmds import (
#     SkewerCmd, HisatCmd, Bowtie2Cmd, FastqcCmd,
#     SamtoolsSortCmd, SamtoolsIndexCmd, BedtoolsMulticovCmd,
# )


ROOT_DIR = path.protect('~/work/projects')
PIPE_PBS_TEMPLATE = ''


def add_subparsers():
    '''Use this function to add subparsers from modules'''
    subparser = base.get_subparser()
    fastq_parser = FastqScripted(subparser)
    fastq_parser.setup()


def setup_logger():
    # setup logger
    import logging
    from remsci.lib.utility import customLogging
    customLogging.config()
    log = logging.getLogger(__name__)
    return log


def parse_args():
    # setup the parser
    add_subparsers()
    parser = base.get_parser()
    subpars = base.get_subparser()
    subpars.required = True

    args = parser.parse_args()

    # parse directory file list
    # (requires input_directory_parser to be included)
    try:
        args.find_files(args)
    except AttributeError:
        pass

    # protect file names
    protect = ['summary', 'root_dir', 'data_dir']
    for p in protect:
        try:
            setattr(args, p, path.protect(getattr(args, p)))
        except TypeError:
            pass  # arg might not be set (is None)

    return args


def read_summary(args):
    '''Reads a summary file and returns a list of names and files

    Expected format:
            <name>  <file>  [<file> ...]
        The first two columns should be the name of the sample
        and the file. Paired-end should be listed on the same line.
    '''

    with open(args.summary, 'r') as fh:
        rows = [line.rstrip().split() for line in fh]
        data_dir = (args.data_dir
                    if args.data_dir else os.path.dirname(fh.name))

    # add and protect path to second (and third) columns in rows
    for row in rows:
        row[1:] = [
            path.protect(os.path.join(data_dir, col))
            for col in row[1:]
        ]

    return rows


def run_pipe(summary, genome, project_dir, force=False):

    # initialize all pipes
    pipes = []
    for row in summary:

        job_name = row[0]
        files = row[1:]
        log.info('Processing "{}"'.format(job_name))

        sample_dir = os.path.join(project_dir, job_name)
        path.makedirs(sample_dir)

        pipe = NestedGenomicsPipe(
            job_name=job_name,
            odir=sample_dir,
            input_list=files,
            force=force,
            genome=genome,
        )
        pipe.write_script(directory=sample_dir)
        pipes.append(pipe)

    # execute each pipe in turn
    for pipe in pipes:
        relpath = os.path.relpath(pipe.pbs_file, ROOT_DIR)
        log.info('Running pbs script: "{}"'.format(relpath))
        pipe.run()


def main():
    global ROOT_DIR
    args = parse_args()
    summary = read_summary(args)
    force = args.force

    if args.root_dir:
        ROOT_DIR = path.protect(args.root_dir)

    genome = args.genome
    project_name = args.project if args.project else 'new_project'
    project_dir = os.path.join(ROOT_DIR, project_name, 'samples')

    run_pipe(summary, genome, project_dir, force=force)


if __name__ == '__main__':

    log = setup_logger()
    main()
