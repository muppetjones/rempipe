
import argparse
import os.path
import re
import sys

print("libpipe init -- add 'remsci' to path")
sys.path.insert(
    0, os.path.join(os.path.dirname(os.path.dirname(
        os.path.dirname(os.path.abspath(__file__)))), 'remsci'))

import remsci.scripted.base as remparse
from remsci.lib.decorators import file_or_handle


import remsci.scripted.base as base
from remsci.scripted.interface import SubparserBase

import logging
log = logging.getLogger(__name__)


class CoverageScripted(SubparserBase):

    def __init__(self, subparser=None):

        if not subparser:
            subparser = base.get_subparser()

        # create the subparser
        self.subparser = subparser.add_parser(
            'coverage',
            parents=[
                base.input_file_parser(),
                base.input_directory_parser(),
                base.output_file_parser(),
            ],
            formatter_class=argparse.ArgumentDefaultsHelpFormatter,
            help='Calculate coverage summaries',
        )

        self.setup()

    def setup(self):
        super().setup()

        self.subparser.add_argument(
            '--min-count', dest='min_count',
            default=1,
            help=('Minimum counts required per feature'),
        )

        return

    def run(self, args):

        # parse directory file list
        # (requires input_directory_parser to be included)
        try:
            dir_path = args.dir_path
            args.find_files(args)
        except argparse.ArgumentError as e:
            print(e.message)
            exit()
        except AttributeError:
            pass

        log.debug('\n'.join(
            os.path.relpath(f, dir_path)
            for f in args.file_list
        ))

        count_list = self.count_files(args.file_list, args.min_count)
        if not args.outfile:
            common = os.path.commonprefix(args.file_list)
            dir_str = os.path.dirname(common)
            args.outfile = os.path.join(dir_str, 'coverage.txt')

        self.write_coverage(args.outfile, count_list)
        return args.outfile

    def count_files(self, file_list, min_count=1):

        count_list = []
        for f in sorted(file_list):
            name = os.path.splitext(os.path.basename(f))[0]
            count_list.append(
                [name, ] + self.count_file(f, min_count)
            )
        return count_list

    @file_or_handle(mode='r')
    def count_file(self, fh, min_count=1):
        '''Count total features in htseq-count output with counts > min

        Arguments:
            fh          FILE    htseq-count output
            min_count   INT     minimum count threshold

        Return:
            [total count, count > min, ratio]
        '''
        total_above = 0
        total_non_feat = 0
        for i, line in enumerate(fh):
            if line.startswith('__'):
                total_non_feat = total_non_feat + 1
                continue
            if int(line.split()[1]) > min_count:
                total_above = total_above + 1
        i = i - total_non_feat + 1
        return [str(i), str(total_above), '{:.2f}'.format(total_above / i)]

    @file_or_handle(mode='w')
    def write_coverage(self, fh, cov_mat):
        fh.write('# Created by:\n# {}\n'.format(' '.join(sys.argv)))
        row_str = '{}\n'
        for row in cov_mat:
            fh.write(row_str.format('\t'.join(row)))
        return


def setup_logger():
    # setup logger
    import logging
    from remsci.lib.utility import customLogging
    customLogging.config()
    log = logging.getLogger(__name__)
    return log

if __name__ == '__main__':

    log = setup_logger()

    parser = remparse.get_parser()
    subpar = remparse.get_subparser()
    CoverageScripted(subpar)

    args = parser.parse_args()
    args.func(args)
