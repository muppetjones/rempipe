
'''Aggregate FastQC summary data into single page

Creates a table of FastQC data pulled from FastQC output.
Each row gives the sample name and the pass/fail status for each test.
Clicking on each row will open a section below the row, displaying the
FastQC generated graphs for the sample.

Required files:
    The FastQC zipped output (for the pass/fail status)
    The FastQC html file (assumed in the same directory) for the
    test results from each section.

TODO: De-couple from Semantic-UI (html template uses stylesheet).
TODO: Recreate for Django. Will greatly simplify the Python code
        and the template dependance.
TODO: Add coloring to cell based on P/F. (use :before & :after for PE)

Caveats:
    Will not function in Chrome without '--allow-file-access-from-files'

'''

import os.path
import re
from zipfile import ZipFile

import libpipe.templates


class FastQCAggregator(object):

    SUMMARY_COL_INDEX = {
        'result': 0,
        'category': 1,
        'file': 2,
    }

    TEMPLATE = None

    def __init__(self, zipped_files):

        if not self.TEMPLATE:
            self.TEMPLATE = os.path.join(
                os.path.dirname(libpipe.templates.__file__),
                'fastqc_aggregate.html'
            )

        self.zipped = zipped_files

    def __str__(self):
        row_fmt = '{:<5}' * len(self.categories)
        row_list = []
        for sample, results in zip(self.samples, self.results):
            row_str = row_fmt.format(*results)
            row_list.append('{:<12}{}'.format(sample, row_str))
        return '\n'.join(row_list)

    def aggregate(self):
        log.debug('aggregate')
        categories = self._read_categories(self.zipped[0])
        results = [
            self._read_summary(zip_file)
            for zip_file in self.zipped
        ]

        prefix = os.path.dirname(os.path.commonprefix(self.zipped))
        samples = [
            os.path.dirname(os.path.relpath(sample, prefix))
            for sample in self.zipped
        ]

        # collapse results for paired end
        if samples[0] == samples[1]:
            samples, results = self._collapse(samples, results)

        self.categories = categories
        self.samples = samples
        self.results = results

    def write_html(self):

        log.debug('write')
        with open(self.TEMPLATE, 'r') as fh:
            template = fh.read()

        html_dir = os.path.dirname(os.path.commonprefix(self.zipped))
        html_name = os.path.basename(os.path.commonprefix(self.zipped))
        html_path = os.path.join(html_dir, html_name + '_fastqc.html')

        template = template.replace('{{ name }}', html_name)

        # Parse out the template
        # > first and last line will be template markers
        # > remove leading and trailing white space
        m = re.search(
            '(.*)BEGIN TEMPLATE.+JOIN="(.+?)"(.+?)END TEMPLATE(.*)', template,
            re.DOTALL)
        if not m:
            raise ValueError('template not found in html')
        prefix_html = m.group(1).split('\n')[:-1]
        postfix_html = m.group(4).split('\n')[1:]
        join_str = m.group(2)
        subtemp = '\n'.join([
            line.rstrip().lstrip()
            for line in m.group(3).split('\n')[1:-1]
        ])

        # Parse join information from the template
        sample_html = [subtemp.replace(
            '{{ test }}',
            ''
        ).replace(
            '{{ results }}',
            join_str.join(self.categories)
        ).replace(
            '{{ data }}',
            '''<p>This table dynamically pulls the FastQC information for the
            samples listed below.<br />Each row lists the test results.
            Click each row to display the FastQC graphs.</p>'''
        ).replace(
            '{{ n_results }}',
            str(len(self.categories) + 1)
        )]
        log.debug(sample_html)
        for sample, results, zipped in zip(self.samples, self.results, self.zipped):
            sample_html.append(subtemp.replace(
                '{{ test }}',
                sample
            ).replace(
                '{{ results }}',
                join_str.join(results)
            ).replace(
                '{{ data }}',
                '<div class="loader" data-html="' +
                # '<div class="ui embed" data-url="' +
                os.path.relpath(zipped.replace('.zip', '.html'), html_dir) +
                '"></div>'
            ).replace(
                '{{ n_results }}',
                str(len(self.categories) + 1)
            ))

        log.debug('Writing to {}'.format(html_path))
        with open(html_path, 'w') as fh:
            fh.write('\n'.join(prefix_html))
            fh.write('\n'.join(sample_html))
            fh.write('\n'.join(postfix_html).replace(  # add command to suffix
                '{{ command }}',
                ' '.join(sys.argv),
            ))

        return

    @classmethod
    def _read_categories(cls, zip_file, target_file='summary.txt'):
        return cls.__read_column(
            zip_file, target_file, cls.SUMMARY_COL_INDEX['category'])

    @classmethod
    def _read_summary(cls, zip_file, target_file='summary.txt'):
        return cls.__read_column(
            zip_file, target_file, cls.SUMMARY_COL_INDEX['result'])

    @classmethod
    def __read_column(cls, zip_file, target_file, col):

        with ZipFile(zip_file, mode='r') as zf:

            try:
                zf.getinfo(target_file)
            except KeyError:
                target_file = os.path.join(
                    os.path.splitext(os.path.basename(zip_file))[0],
                    target_file,
                )

            with zf.open(target_file, 'r') as sf:
                categories = [
                    line.decode('utf8').lstrip().rstrip().split('\t')[col]
                    for line in sf
                ]
                return categories

    @classmethod
    def _collapse(cls, samples, results):
        '''Collapse every two rows into single row -- paired-end'''

        first_set = results[0::2]
        second_set = results[1::2]
        samples = samples[0::2]

        combined = [
            [
                '{}/{}'.format(x[0], y[0])
                for x, y in zip(first_set[i], second_set[i])
            ]
            for i in range(len(first_set))
        ]
        return samples, combined


def setup_logger():
    # setup logger
    import logging
    from remsci.lib.utility import customLogging
    customLogging.config()
    log = logging.getLogger(__name__)
    return log

if __name__ == '__main__':
    import sys
    from os.path import dirname, abspath, join

    print("libpipe init -- add 'remsci' to path")
    sys.path.insert(
        0, join(dirname(dirname(dirname(abspath(__file__)))), 'remsci'))

    from remsci.lib.utility import path

    log = setup_logger()

    read_dir = '/Users/biiremployee/work/projects/laura_wgs/samples'
    zipped_files = path.walk_file(read_dir, pattern=r'trimmed.*\.zip')

    log.debug('foo')
    fastqc_agg = FastQCAggregator(zipped_files)
    log.debug('bar')
    fastqc_agg.aggregate()

    fastqc_agg.write_html()
