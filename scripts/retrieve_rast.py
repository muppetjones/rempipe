
import os.path
import pathlib
import re
import subprocess
import sys
import time

from itertools import cycle
from textwrap import dedent

print("libpipe init -- add 'remsci' to path")
rempath = pathlib.PurePath('~/dev/remsci')
rempipe = pathlib.PurePath('~/dev/pipe')
sys.path.insert(0, os.path.expanduser(str(rempath)))
sys.path.insert(0, os.path.expanduser(str(rempipe)))

# join(dirname(dirname(dirname(abspath(__file__)))), 'remsci'))

from remsci.lib.utility import path
from remsci.lib.decorators import file_or_handle

import logging
log = logging.getLogger(__name__)


FORMAT = 'rast_tarball'
FORMAT_EXT = {
    'genbank': '.gb',
    'embl': '.embl',
    'gff3': '.gff3',
    'rast_tarball': '.tgz',
}
ROOT_DIR = '/Users/biiremployee/work/projects/lzd_babies_wgs/samples'


def usage():
    msg = dedent('''
        Usage:
            python retrieve_rast.py <root dir> <format> [pattern]

        Synopsis:
            Identify files containing RAST job ids, parse for the id, and
            retrieve the job in the requested format.

        Arguments:
            root dir: The root directory where the RAST job ids may be found.
            format: Retrieval format (genbank|embl|gff3|rast_tarball[default]).
                Default: tarball.
            pattern: Optional pattern for files containing rast id.
                Default: "rast.txt".
    ''')
    return msg


def get_args(idx, default=None):
    '''REPLACE THIS WITH ARGPARSE!!'''

    try:
        return sys.argv[idx]
    except IndexError:
        if default:
            return default
        else:
            print(usage())
            raise ValueError('Bad arguments')


def get_rast_ids(root_dir, pattern='rast.txt'):

    file_list = path.walk_file(root_dir, pattern=pattern)
    rast_ids = list(filter(None, [parse_rast_id(f) for f in file_list]))
    bad_ids = [f for f in file_list if parse_rast_id(f) is None]

    if bad_ids:
        log.warning('Bad ids:\n  {}'.format('\n  '.join(bad_ids)))

    return rast_ids, file_list


@file_or_handle(mode='r')
def parse_rast_id(fh):

    txt = fh.read()
    m = re.search(r"Job '(?P<id>\d+)'", txt)
    try:
        return m.group('id')
    except:
        return None


def retrieve_jobs(rast_dict, fmt, root_dir, rast_config):

    root_path = pathlib.Path(root_dir)

    id_loop = cycle(sorted(rast_dict.keys()))
    while rast_dict:

        rast_id = next(id_loop)
        download_file = parse_output_filename(
            rast_id, rast_dict, root_path, fmt)
        if os.path.isfile(str(download_file)):
            print('Output found: {}. Skipping'.format(download_file))
            del rast_dict[rast_id]
            id_loop = cycle(rast_dict)
            continue

        if check_job_status(rast_id, rast_config):
            try:
                retrieve_job(rast_id, fmt, download_file, rast_config)
            except Exception as e:
                log.debug(str(e))
                raise
                pass  # just continue on to the next job
            else:
                del rast_dict[rast_id]
                id_loop = cycle(sorted(rast_dict.keys()))


def parse_output_filename(rast_id, rast_dict, root_path, fmt):
    # parse the download filename
    # -- sample_annotation.format
    # -- assume sample name is the first dict after root
    infile = rast_dict[rast_id]
    file_path = pathlib.Path(infile)

    file_relpath = file_path.relative_to(root_path)
    assumed_sample_name = file_relpath.parts[0]
    print('Checking sample {}'.format(assumed_sample_name))

    download_file = file_path.parent.joinpath(
        '{}_annotation{}'.format(
            assumed_sample_name, FORMAT_EXT[fmt]
        ))
    return download_file


def check_job_status(rast_id, rast_config):
    print('Checking job status [{}]'.format(rast_id))

    check_status_call = "svr_status_of_RAST_job {uname} {passwd} {job_id}"

    cmd_dict = {
        'job_id': rast_id,
    }
    cmd_dict.update(rast_config)

    cmd_chk = check_status_call.format(**cmd_dict).split()

    try:
        status = subprocess.check_output(cmd_chk, universal_newlines=True)
    except:
        raise

    if 'complete' in status:
        return True

    return False


def retrieve_job(rast_id, fmt, download_file, rast_config):
    print('Retrieving [{}] to {}'.format(rast_id, download_file))
    log_file = download_file.with_suffix('.log')

    retrieve_call = "svr_retrieve_RAST_job {uname} {passwd} {job_id} {format}"

    cmd_dict = {
        'job_id': rast_id,
        'format': fmt,
    }
    cmd_dict.update(rast_config)

    cmd_ret = retrieve_call.format(**cmd_dict).split()
    with open(str(download_file), 'w') as oh, open(str(log_file), 'w') as lh:
        lh.write(' '.join(cmd_ret))
        try:
            subprocess.check_call(cmd_ret, stdout=oh, stderr=lh)
        except Exception as e:
            lh.write(str(e))
            raise
        else:
            return True


def main(rast_config):

    root_dir = get_args(1, None)
    fmt = get_args(2, 'rast_tarball')
    pattern = get_args(3, 'rast.txt')

    rast_ids, file_list = get_rast_ids(root_dir, pattern)
    rast_lookup = {k: v for k, v in zip(rast_ids, file_list)}

    # next, retrieve the rast ids
    # -- rm the id if successful, keep trying otherwise
    # -- need some way to mark if done (for future script runs)
    retrieve_jobs(rast_lookup, fmt, root_dir, rast_config)

if __name__ == '__main__':

    from remsci.lib.utility import customLogging
    customLogging.config()

    # load config info
    import json
    import libpipe.config
    with open(os.path.join(os.path.dirname(
            libpipe.config.__file__), 'rast.config')) as fh:
        rast_config = json.load(fh)

    main(rast_config)
