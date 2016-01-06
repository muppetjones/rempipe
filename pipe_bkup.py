#!/usr/bin/env python

import os.path
import subprocess
import sys
import time
sys.path.append("../remsci/")

import remsci.scripted.base as base
from remsci.lib.utility import path
from libpipe.pipes.base import BasePipe
from libpipe.parsers.fastq import FastqScripted
from libpipe.cmds import (
    SkewerCmd, HisatCmd, Bowtie2Cmd, FastqcCmd,
    SamtoolsSortCmd, SamtoolsIndexCmd, BedtoolsMulticovCmd,
)


ROOT_DIR = path.protect('~/data/rempipe')
PIPE_PBS_TEMPLATE = ''
DO_RUN = lambda x: not _has_output(x)


def add_subparsers():
    '''Use this function to add subparsers from modules'''
    subparser = base.get_subparser()
    fastq_parser = FastqScripted(subparser)
    fastq_parser.setup()


def pipe(file_list, genome, project_dir, force=False):

    timestamp = time.strftime("%y%m%d-%H%M%S")

    for f in file_list:
        name = f[0]
        files = f[1:]
        out_dir = os.path.join(project_dir, name)
        path.makedirs(out_dir)

        # # 1st fast qc
        # fastqc_1 = FastqcCmd(*files, o=out_dir)
        #
        # # trimming
        # out_prefix = os.path.join(out_dir, name)
        # trim = SkewerCmd(*files, o=out_prefix)
        # trimmed_fastq = trim.output()
        #
        # # 2nd fastqc
        # fastqc_2 = FastqcCmd(*trimmed_fastq, o=out_dir)
        #
        # # setup alignment
        # # NOTE: need to check for encoding
        # align_kwargs = {
        #     '-x': genome,
        #     '-S': '{}_{}.sam'.format(
        #         out_prefix,
        #         os.path.basename(genome),
        #     ),
        #     '-p': 3,  # set for local (should use pbs paramters on qsub)
        # }
        # if len(trimmed_fastq) == 1:
        #     align_kwargs['U'] = trimmed_fastq[0]
        # else:
        #     align_kwargs['1'], align_kwargs['2'] = trimmed_fastq
        # align = HisatCmd(timestamp=timestamp, **align_kwargs)
        # # human_kw = [m for m in ['human', 'sapien', 'G37RCh'] if m in genome]
        # # if human_kw:
        # #     align = HisatCmd(timestamp=timestamp, **align_kwargs)
        # # else:
        # #     align = Bowtie2Cmd(timestamp=timestamp, **align_kwargs)
        #
        # # samtools
        # sam_sort = SamtoolsSortCmd(*(align.output()))
        # sam_index = SamtoolsIndexCmd(*(sam_sort.output()))

        # UPDATED
        fastqc_1 = FastqcCmd(*files, o=out_dir)

        # trimming
        out_prefix = os.path.join(out_dir, name)
        trim = SkewerCmd(*files, o=out_prefix)

        # 2nd fastqc
        fastqc_2 = FastqcCmd(o=out_dir)

        # setup alignment
        # NOTE: need to check for encoding
        align_kwargs = {
            '-x': genome,
            '-S': '{}_{}.sam'.format(
                out_prefix,
                os.path.basename(genome),
            ),
            '-p': 3,  # set for local (should use pbs paramters on qsub)
        }
        align = HisatCmd(timestamp=timestamp, **align_kwargs)

        # samtools
        sam_sort = SamtoolsSortCmd()
        sam_index = SamtoolsIndexCmd()

        # count
        kwargs = {'-bed': genome}
        bedtools_multicov = BedtoolsMulticovCmd(**kwargs)

        # Setup pipe
        # NOTE: This is the alpha test of the pipe class.
        job_name = name + '_' + os.path.basename(genome)
        pipe = BasePipe(job_name=job_name, force=force)
        pipe.add(
            fastqc_1, trim, fastqc_2, align, sam_sort, sam_index,
            bedtools_multicov,
        )

        # write pbs file & run
        # pbs_file = '{}  , timestamp, os.path.basename(genome))
        pipe.write_script()
        pipe.run()
        # cmds = [fastqc_1, trim, fastqc_2, align, sam_sort, sam_index]
        # _write_pbs(pbs_file, cmds)

        # run pbs in qsub
        # qid = _run_pbs(name, pbs_file)
        # log.debug("Running as '{}'".format(qid))


def _write_pbs(pbs_file, cmds):
    global PIPE_PBS_TEMPLATE
    if not PIPE_PBS_TEMPLATE:
        with open('./template.pbs', 'r') as ih:
            PIPE_PBS_TEMPLATE = ih.read()

    # write pbs file
    with open(pbs_file, 'w') as oh:
        oh.write(PIPE_PBS_TEMPLATE + "\n")

        omit_msg = '# Output found. Skip.\n# '
        comment_str = 'echo "Running {}..."\n'

        for cmd in cmds:
            cmd_str = str(cmd) + "\n\n"
            if DO_RUN(cmd):
                oh.write(comment_str.format(cmd.name))
                oh.write(cmd_str)
            else:
                oh.write(omit_msg + cmd_str)


def _run_pbs(job_name, pbs_file, log_file=''):

    qid_file = pbs_file.replace('pbs', 'qid')
    if not log_file:
        log_file = pbs_file.replace('pbs', 'log')

    # run pbs file
    qsub_command = ['qsub', '-N', job_name, '-o', log_file, pbs_file, ]
    with open(qid_file, 'w') as qid_fh:
        log.debug(' '.join(qsub_command))
        retcode = subprocess.check_call(qsub_command, stdout=qid_fh)
        log.debug('retcode: {}'.format(retcode))

    # get the qsub job id
    with open(qid_file, 'r') as qid_fh:
        qid = qid_fh.read()

    # link the log file (last run)
    link_base = os.path.join(os.path.dirname(pbs_file), job_name)
    _renew_link(link_base + '.log', log_file)
    _renew_link(link_base + '.pbs', pbs_file)
    _renew_link(link_base + '.qid', qid_file)

    return qid


def _renew_link(link_name, src_file):
    if os.path.isfile(link_name):
        os.unlink(link_name)
    os.symlink(src_file, link_name)


def _has_output(cmd):
    '''Checks for expected output from a command

    Arguments:
        cmd     Command object.
    Returns:
        True if ALL output is found.
        False otherwise
    '''

    for f in cmd.output:
        if not os.path.isfile(f):
            return False
    return True


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
    if args.root_dir:
        ROOT_DIR = path.protect(args.root_dir)

    if args.force:
        DO_RUN = lambda x: True

    project_dir = os.path.join(ROOT_DIR, args.project, 'samples')
    pipe(file_list, args.genome, project_dir, force=args.force)
