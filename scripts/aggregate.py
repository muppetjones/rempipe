'''Summarize alignment and count info for samples



'''

import abc
import os.path
import re
import sys
from itertools import cycle


print("libpipe init -- add 'remsci' to path")
sys.path.insert(
    0, os.path.join(os.path.dirname(os.path.dirname(
        os.path.dirname(os.path.abspath(__file__)))), 'remsci'))

import remsci.scripted.base as remparse
from remsci.lib.decorators import file_or_handle

import libpipe.templates
from libpipe.parsers.aggregate import AggregateScripted


class Aggregator(metaclass=abc.ABCMeta):
    TEMPLATE = None
    TEMPLATE_FILE = 'aggregator.html'
    SUFFIX = None

    def __init__(self, agg_files):

        if not self.TEMPLATE:
            self.TEMPLATE = os.path.join(
                os.path.dirname(libpipe.templates.__file__),
                self.TEMPLATE_FILE
            )

        self.agg_files = sorted(agg_files)
        self.write_order = self.WRITE_ORDER
        self.suffix = self.SUFFIX

    @abc.abstractmethod
    def aggregate(self):

        # first, we need to find out how may times the header will be split
        # --usually the number of underscores, but we remove 'unal'
        #   and timestamps and everything following, so...
        # HACK: this is ugly and inefficient (runs each header twice!)
        self._split_headers()
        self._order_cols()

    def combine(self, agg):
        new_summary = [
            [list(rowA[0]), list(rowA[1]) + list(rowB[1])]
            for rowA, rowB in zip(self.summary, agg.summary)
            if rowA[0] == rowB[0]
        ]

        self.summary = new_summary
        self.write_order = self.write_order + agg.write_order

        self.suffix = ''.join(s.capitalize()
                              for s in [self.suffix, agg.suffix])

    def write_html(self):

        html_dir, html_name, html_path = self._parse_html_name()

        # sets self.col_sep and self.template
        base_html = self._parse_template(html_name)

        # write order should be set in aggregate
        try:
            rows = [self._create_row(None, self.write_order), ]
        except AttributeError:
            raise NotImplementedError('Aggregator MUST set WRITE_ORDER')

        previous = ''
        color_cycle = cycle(['negative', ''])
        for row in self.summary:
            sample_group = row[0][0]
            if sample_group != previous:
                color = next(color_cycle)
                previous = sample_group

            row = self._create_row(
                row[0], row[1], class_name=color
            )
            rows.append(row)

        base_html = base_html.replace('{{ template }}', '\n'.join(rows))

        with open(html_path, 'w') as fh:
            fh.write(base_html)

    def _parse_html_name(self):
        suffix = '_{}.html'.format(self.suffix)

        html_dir = os.path.dirname(os.path.commonprefix(self.agg_files))
        html_name = os.path.basename(os.path.commonprefix(self.agg_files))
        html_path = os.path.join(html_dir, html_name + suffix)

        return html_dir, html_name, html_path

    def _parse_template(self, name):
        with open(self.TEMPLATE, 'r') as fh:
            html = fh.read()

        # Parse out the template
        # > first and last line will be template markers
        # > remove leading and trailing white space
        m = re.search(
            '(.*)BEGIN TEMPLATE.+JOIN="(.+?)"(.+?)END TEMPLATE(.*)', html,
            re.DOTALL)
        if not m:
            raise ValueError('template not found in html')
        prefix_html = m.group(1).split('\n')[:-1]
        suffix_html = m.group(4).split('\n')[1:]
        join_str = m.group(2)
        template = '\n'.join([
            line.rstrip().lstrip()
            for line in m.group(3).split('\n')[1:-1]
        ])

        # add command to suffix
        suffix_html = "\n".join(suffix_html)
        suffix_html = suffix_html.replace(
            '{{ command }}',
            ' '.join(sys.argv),
        ).lstrip().rstrip()

        base_html = (
            "\n".join(prefix_html) +
            "\n{{ template }}\n" +
            suffix_html
        )
        base_html = base_html.replace('{{ name }}', name)

        self.template, self.col_sep = template, join_str
        return base_html

    def _filter_header(self, header):
        header_list = [
            h for h in re.split(r'[_\-]', header)
            if h != 'unal'
        ]
        # keep up to timestamp
        last_index = 0
        for i, h in enumerate(header_list):
            last_index = i + 1
            try:
                int(h)
            except ValueError:
                continue
            else:
                last_index = last_index - 1
                break
        return header_list[0:last_index]

    def _split_headers(self):
        new_summary = [
            [self._filter_header(k), v]
            for k, v in sorted(self.summary.items())
        ]
        hmax = [len(row[0]) for row in new_summary]
        self.summary = new_summary
        self.header_len = max(hmax)
        return

    def _order_cols(self):
        # assumes header already split
        for row in self.summary:
            row[1] = [row[1][o] for o in self.write_order]

    def _create_row(self, header, row, class_name=''):

        # HACK
        header_list = [''] * self.header_len
        if header:
            header_list = header + header_list

        header = '</th><th>'.join(header_list[:self.header_len])

        return self.template.replace(
            '{{ class_name }}', class_name
        ).replace(
            '{{ sample }}', header
        ).replace(
            '{{ results }}', self.col_sep.join(row)
        )

    @abc.abstractmethod
    def _parse_summary(self, fh):
        pass


class AlignmentAggregator(Aggregator):
    SUFFIX = 'align'
    WRITE_ORDER = [
        'total_reads', 'total_aligned', 'total align %', 'aligned_once', 'align once %']

    def aggregate(self):
        summary = {}
        for agg in self.agg_files:
            base_name = os.path.splitext(os.path.basename(agg))[0]
            base_name = base_name.replace('-trimmed_', '_')
            if agg.endswith('.log'):
                summary[base_name] = self._parse_summary(agg)

        self.summary = summary
        super().aggregate()

    @file_or_handle(mode='r')
    def _parse_summary(self, fh):
        '''Read and parse tuxedo family alignment summary

        Example:
            575434 reads; of these:
              575434 (100.00%) were unpaired; of these:
                566497 (98.45%) aligned 0 times
                8088 (1.41%) aligned exactly 1 time
                849 (0.15%) aligned >1 times
            1.55% overall alignment rate
        '''

        summary = [
            line.lstrip().rstrip().split()
            for line in fh
        ]

        aligned_once = summary[3][0]
        aligned_mult = summary[4][0]
        total_reads = summary[0][0]
        total_aligned = str(int(aligned_once) + int(aligned_mult))
        perc_total = summary[-1][0]
        perc_once = summary[3][1][1:-1]
        summary_details = {
            'total_reads': total_reads,
            'total_aligned': total_aligned,
            'total align %': perc_total,
            'aligned_once': aligned_once,
            'align once %': perc_once,
        }

        return summary_details


class CountAggregator(Aggregator):

    SUFFIX = 'count'
    WRITE_ORDER = [
        'not_aligned',
        'no_feature',
        'ambiguous',
        'too_low_aQual',
        'alignment_not_unique',
    ]

    def aggregate(self):
        summary = {}
        for agg in self.agg_files:
            base_name = os.path.splitext(os.path.basename(agg))[0]
            base_name = base_name.replace('-trimmed_', '_')
            if agg.endswith('.count'):
                summary[base_name] = self._parse_summary(agg)

        self.summary = summary
        super().aggregate()

    @file_or_handle(mode='r')
    def _parse_summary(self, fh):
        '''Read and parse tuxedo family alignment summary

        Example:
            ...
            __no_feature	211341
            __ambiguous	2709
            __too_low_aQual	3916
            __not_aligned	304930
            __alignment_not_unique	69314
        '''

        summary = [
            line.lstrip().rstrip().split()
            for line in fh if line.startswith('__')
        ]

        summary_details = {
            k[2:]: v
            for k, v in summary
        }

        return summary_details


def setup_logger():
    # setup logger
    import logging
    from remsci.lib.utility import customLogging
    customLogging.config()
    log = logging.getLogger(__name__)
    return log


def add_subparsers(subparser):
    '''Use this function to add subparsers from modules'''
    AggregateScripted(subparser)


if __name__ == '__main__':
    import sys
    from os.path import dirname, abspath, join

    print("libpipe init -- add 'remsci' to path")
    sys.path.insert(
        0, join(dirname(dirname(dirname(abspath(__file__)))), 'remsci'))

    from remsci.lib.utility import path

    log = setup_logger()

    parser = remparse.get_parser()
    subpar = remparse.get_subparser()
    add_subparsers(subpar)

    args = parser.parse_args()
    file_list = args.func(args)

    # read_dir = '/Users/biiremployee/work/projects/lzd_babies_rnaseq/samples'
    # count_files = path.walk_file(read_dir, pattern=r'\.count[^\.]')
    # align_files = path.walk_file(read_dir, pattern=r'hisat2?\.log')
    #

    if args.alignment_summary:
        aggA = AlignmentAggregator(file_list)
        aggA.aggregate()
        aggA.write_html()

    if args.count_summary:
        aggC = CountAggregator(file_list)
        aggC.aggregate()
        aggC.write_html()

    if args.alignment_summary and args.count_summary:
        aggA.combine(aggC)
        aggA.write_html()
