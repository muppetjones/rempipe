#!/usr/bin/env python

import os.path
import subprocess
import sys
import time
sys.path.append("../remsci/")

import remsci.scripted.base as base
from remsci.lib.utility import path
from libpipe.parsers.fastq import FastqScripted
from libpipe.cmds import (SkewerCmd, HisatCmd, Bowtie2Cmd)


ROOT_DIR = path.protect('~/data/rempipe')


def add_subparsers():
    '''Use this function to add subparsers from modules'''
    subparser = base.get_subparser()
    fastq_parser = FastqScripted(subparser)
    fastq_parser.setup()


def pipe(file_list, genome, project_dir):

    with open('./template.pbs', 'r') as ih:
        pbs = ih.read()

    timestamp = time.strftime("%y%m%d-%H%M%S")

    for f in file_list:
        name = f[0]
        files = f[1:]
        out_dir = os.path.join(project_dir, name)
        log.debug(out_dir)
        path.makedirs(out_dir)

        # trimming
        out_prefix = os.path.join(out_dir, name)
        skewer = SkewerCmd(*files, o=out_prefix)

        # setup alignment
        # NOTE: need to check for encoding
        hisat_kwargs = {
            'x': genome,
            'S': '{}_{}.sam'.format(
                out_prefix,
                os.path.basename(genome),
            ),
        }
        trimmed_fastq = skewer.output()
        if len(trimmed_fastq) == 1:
            hisat_kwargs['U'] = trimmed_fastq[0]
        else:
            hisat_kwargs['1'], hisat_kwargs['2'] = trimmed_fastq
        hisat = Bowtie2Cmd(timestamp=timestamp, **hisat_kwargs)

        # write pbs file
        pbs_file = '{}_pipe_{}.pbs'.format(out_prefix, timestamp)
        log_file = pbs_file.replace('pbs', 'log')
        with open(pbs_file, 'w') as oh:
            file_content = '\n'.join([
                pbs, skewer.cmd(), hisat.cmd(),
            ]) + '\n'

            oh.write(file_content)

        # run pbs file
        qid_file = '{}_pipe_{}.qid'.format(out_prefix, timestamp)
        qsub_command = ['qsub', '-N', name, '-o', log_file, pbs_file,]
        with open(qid_file, 'w') as qid:
            log.debug(' '.join(qsub_command))
            retcode = subprocess.check_call(qsub_command, stdout=qid)
            log.debug('retcode: {}'.format(retcode))
            

def main():
    pass


if __name__ == '__main__':

    # setup logger
    import logging
    from remsci.lib.utility import customLogging
    customLogging.config()
    log = logging.getLogger(__name__)

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

    # call the default function
    try:
        file_list = args.func(args)
    except AttributeError:
        raise
        pass

    # update root directory
    if args.root:
        ROOT_DIR = path.protect(args.root)

    project_dir = os.path.join(ROOT_DIR, 'samples', args.project)
    pipe(file_list, args.genome, project_dir)
