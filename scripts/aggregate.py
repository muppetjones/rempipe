'''Summarize alignment and count info for samples
'''

import abc
import numpy
import os.path
import re
import sys
from itertools import cycle

from libpipe.decorators import file_or_handle

import libpipe.templates
from libpipe.parsers.aggregate import AggregateScripted
from libpipe.parsers import base as remparse


def setup_logger():
    # setup logger
    import logging
    from libpipe.utility import logging as pipelog
    pipelog.config()
    log = logging.getLogger(__name__)
    return log

log = setup_logger()
log.debug('hello world')


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

    def combine(self, *args):
        for agg in args:
            self._combine(agg)
        return
        new_summary = self.summary
        new_summary = [
            [list(rowA[0]), list(rowA[1]) + list(rowB[1])]
            for rowA, rowB in zip(new_summary, agg.summary)
            if rowA[0] == rowB[0]
        ]

        self.summary = new_summary
        self.write_order = self.write_order + agg.write_order

        self.suffix = ''.join(s.capitalize()
                              for s in [self.suffix, agg.suffix])

    def _combine(self, agg):

        for this, that in zip(self.summary, agg.summary):
            if this[0] == that[0]:
                this[1].extend(that[1])
        self.write_order.extend(agg.write_order)
        self.suffix = '{}{}'.format(
            self.suffix, agg.suffix.capitalize())
        return

    def _calc_stats(self):

        log.debug([' '.join(row[0][1:])
                   for row in self.summary])

        unique_grp = sorted(list(set(
            ' '.join(row[0][1:])
            for row in self.summary
        )))

        stats = [numpy.mean, numpy.median]
        stat_rows = []
        for stat in stats:
            stat_rows.extend(self._calc_stat(stat, unique_grp))
        return stat_rows

    def _calc_stat(self, stat, unique_grp):
        stats_rows = []
        log.debug(stat)
        log.debug(unique_grp)
        for grp in unique_grp:

            # get relevant data
            # also convert string to number
            rows = [
                row[1]
                for row in self.summary
                if ' '.join(row[0][1:]) == grp]
            for i, row in enumerate(rows):
                row = [
                    float(cell) if not cell.endswith('%') else float(cell[:-1])
                    for cell in row
                ]
                rows[i] = row
            mat = numpy.array(rows)

            stat_list = stat(mat, axis=0)

            grp_tup = [stat.__name__] + grp.split()
            stats_rows.append([
                grp_tup, ['{:.2f}'.format(v) for v in stat_list]
            ])
        return stats_rows

    def write_html(self, omit=[]):

        html_dir, html_name, html_path = self._parse_html_name()

        # sets self.col_sep and self.template
        base_html = self._parse_template(html_name)

        # update header_len
        no_omit = [elem for elem in self.summary[1][0] if elem not in omit]
        self.header_len = len(no_omit)

        for row in self.summary:
            row[0] = [x for x in row[0] if x not in omit]

        # write order should be set in aggregate
        try:
            rows = [self._create_row(None, self.write_order, omit=omit), ]
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
                row[0], row[1], class_name=color, omit=omit
            )
            rows.append(row)

        # log.debug('\n'.join(rows))

        stats_rows = self._calc_stats()
        for row in stats_rows:
            sample_group = row[0][0]
            if sample_group != previous:
                color = next(color_cycle)
                previous = sample_group
            row = self._create_row(
                row[0], row[1], class_name=color, omit=omit,
            )
            rows.append(row)
        # log.debug('\n'.join(rows))

        base_html = base_html.replace('{{ template }}', '\n'.join(rows))

        with open(html_path, 'w') as fh:
            fh.write(base_html)
        return html_path

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

        try:
            index = header.index('trimmed')
        except ValueError:
            index = 0

        prefix = header[:index].replace('-', '_').rstrip('-_ ')
        header_list = [prefix] + [
            h for h in re.split(r'[_\-]', header[index:])
            if h != 'unal'
            and h not in ['trimmed', 'pair', 'hisat2']
        ]
        # header_list = [
        #     v for v in header_list
        # ]

        # keep up to timestamp
        last_index = len(header_list)
        for i, h in enumerate(header_list):
            try:
                int(h)
                if len(h) == 6:
                    last_index = i
                    break
            except ValueError:
                pass
        return header_list[:last_index]

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

    def _create_row(self, header, row, class_name='', omit=[]):

        # HACK
        header_list = [''] * self.header_len
        if header:
            header_list = header + header_list

        try:
            header = '</th><th>'.join(header_list[:self.header_len])
        except TypeError:
            raise

        try:
            row = [
                '{:,}'.format(int(cell)) if '.' not in cell else
                '{:,}'.format(float(cell)) if not cell.endswith('%')
                else '{:,}'.format(float(cell[:-1]))
                for cell in row
            ]
        except ValueError:
            pass

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
        'total_features',
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

        counts = [line for line in fh]
        summary = [
            line.lstrip().rstrip().split()
            for line in counts if line.startswith('__')
        ]

        summary_details = {
            k[2:]: v
            for k, v in summary
        }

        summary_details['total_features'] = str(len(counts) - 5)
        del counts

        return summary_details


class CoverageAggregator(Aggregator):

    SUFFIX = 'cov'
    WRITE_ORDER = [
        'total_features',
        'features > 1',
        'ratio'
    ]

    def aggregate(self):
        if len(self.agg_files) > 1:
            raise ValueError('too many files. only one needed')
        self.summary = self._parse_summary(self.agg_files[0])
        super().aggregate()
        return

    @file_or_handle(mode='r')
    def _parse_summary(self, fh):
        '''Parse output of 'coverage.py' to include in aggregate'''

        # lines = [line for line in ]
        coverage = [
            line.lstrip().rstrip().split()
            for line in fh if not line.startswith('#')
        ]

        summary_details = {
            line[0].replace('-trimmed_', '_'): {
                k: v
                for k, v in zip(self.WRITE_ORDER, line[1:])
            }
            for line in coverage
        }

        return summary_details


def add_subparsers(subparser):
    '''Use this function to add subparsers from modules'''
    AggregateScripted(subparser)


if __name__ == '__main__':

    log = setup_logger()
    log.debug('hello world')

    parser = remparse.get_parser()
    subpar = remparse.get_subparser()
    add_subparsers(subpar)

    args = parser.parse_args()
    file_list = args.func(args)

    agg_list = []

    omit = ['trimmed', 'pair']

    def aggregate_with(file_list, agg_cls):
        agg = agg_cls(file_list)
        agg.aggregate()
        agg.write_html(omit=omit)
        return agg

    log.debug('alignment')
    if args.alignment_summary:
        agg_list.append(aggregate_with(file_list, AlignmentAggregator))

    log.debug('count')
    if args.count_summary:
        agg_list.append(aggregate_with(file_list, CountAggregator))

    log.debug('finished alignment and count')
    common_dir = os.path.dirname(os.path.commonprefix(file_list))
    cov_file = os.path.join(common_dir, 'coverage.txt')
    if os.path.isfile(cov_file):
        agg_list.append(aggregate_with([cov_file], CoverageAggregator))
    log.debug('finished cov')

    for agg in agg_list[1:]:
        agg_list[0].combine(agg)
    log.debug('finished combinging')
    ofile = agg_list[0].write_html(omit=omit)