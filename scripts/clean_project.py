#!/usr/bin/env python

'''Clean up project files

Removes known intermediate files to preserve disk space.
> Compress FASTQ files
> Convert SAM to BAM; delete SAM

Arguments:
    source directory: Directory to clean
    destination directory: [optional] Directory to move FASTQ files to after
        zipping.
'''

import argparse
import os
import subprocess
import sys
import tarfile

from textwrap import dedent


# print("libpipe init -- add 'remsci' to path")
# sys.path.insert(
#     0, os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(
#         os.path.abspath(__file__)))), 'remsci'))

# from remsci.lib.utility import path

from libpipe.cmds.utility import SamtoolsSortCmd, SamtoolsViewCmd

import logging
log = logging.getLogger(__name__)


def setup_logger():
    # setup logger
    import logging
    from libpipe.utility import logging as pipelog
    pipelog.config()
    log = logging.getLogger(__name__)
    return log


def setup_parser():

    parser = argparse.ArgumentParser(description=dedent('''
        Clean a genomics project.
        > Compress FASTQ files (and optionally move to destination).
        > Convert SAM to BAM.
        > Optionally delete originals.
    '''))
    parser.add_argument(
        'src_dir', help="Source directory (to clean)")
    parser.add_argument(
        'dest_dir', nargs='?', default=None,
        help="[optional] Destination directory for FASTQ files")
    parser.add_argument(
        '--unlink', action='store_true', default=False,
        help='Delete original files.')
    parser.add_argument(
        '--force', action='store_true', default=False,
        help='Overwrite existing output.')
    return parser


def find_fastq_files(srcdir):
    return path.walk_file(srcdir, pattern=r'\.f(?:ast)?q$')


def find_sam_files(srcdir):
    return path.walk_file(srcdir, pattern=r'\.bam$')


def compress_files(file_list, unlink=False, force=False, destdir=None):

    # using the dsrc lossless compression algorithm
    # -- best mode with CRC32 checksum
    # http://sun.aei.polsl.pl/
    cmd = 'dsrc c -m2 -c {raw_file} {compressed_file}'
    ext = '.dsrc'

    for f in file_list:
        log.debug('Compressing {}'.format(os.path.basename(f)))
        if destdir:
            tar_name = os.path.join(destdir, os.path.basename(f) + ext)
        else:
            tar_name = f + ext
        tar_cmd = cmd.format(
            raw_file=f,
            compressed_file=tar_name,
        )

        log.debug(tar_cmd)
        # don't recompress if we don't need to
        if not force and os.path.isfile(tar_name):
            continue

        try:
            subprocess.check_call(tar_cmd, shell=True)
        except Exception as e:
            log.error(str(e))
        else:
            os.unlink(f)
    return


def move_files(file_list, destdir):
    for f in file_list:
        shutil.move(f, destdir)
    return


def convert_to_bam(file_list, unlink=False, force=False):

    for f in file_list:

        log.debug('Converting {}'.format(f))
        sam_sort = SamtoolsSortCmd(f)
        sam_view = SamtoolsViewCmd(f, '-b')

        sort_cmd = sam_sort.cmd()
        view_cmd = sam_view.cmd()

        if not force and (sam_sort._has_output() or sam_view._has_output()):
            continue  # already have output, just skip

        # should run using the cmd itself once implemented
        try:
            subprocess.check_call(sort_cmd.split())
        except Exception as e:
            log.error(str(e))
        else:
            os.unlink(f)


if __name__ == '__main__':

    log = setup_logger()
    parser = setup_parser()
    args = parser.parse_args()

    args.src_dir = args.src_dir.rstrip('/\\')
    try:
        args.dest_dir = args.dest_dir.rstrip('/\\')
    except AttributeError:
        pass  # dest_dir not given

    log.info('Source directory: {}'.format(args.src_dir))
    log.info('Destination directory: {}'.format(args.dest_dir))

    fastq_list = find_fastq_files(args.src_dir)
    sam_list = find_sam_files(args.src_dir)

    log.debug(fastq_list)

    if args.dest_dir:
        common_prefix = os.path.commonprefix(fastq_list)
        args.dest_dir = os.path.join(
            args.dest_dir,
            os.path.basename(args.src_dir),
            os.path.basename(os.path.dirname(common_prefix)),
        )
        try:
            os.makedirs(args.dest_dir)
        except OSError as e:
            log.debug(str(e))
            pass  # we don't care if the directory already exists
        log.info('Moving FASTQ files to {}'.format(args.dest_dir))

    compress_files(
        fastq_list, unlink=args.unlink, force=args.force,
        destdir=args.dest_dir)
    # convert_to_bam(
    #     sam_list, unlink=args.unlink, force=args.force)
