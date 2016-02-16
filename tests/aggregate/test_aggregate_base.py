
import io
import os.path
import unittest
import re
from textwrap import dedent
from unittest.mock import patch

from remsci.lib.tests.base import BaseTestCase


from libpipe.aggregate.base import Aggregator


import logging
log = logging.getLogger(__name__)

sample_aggregate = dedent('''
    row1_a    1
    row2_b    2
    row3_c    3
''').rstrip().lstrip()


class ExampleAggregator(Aggregator):

    def aggregate(self):

        data = sample_aggregate.split('\n')

        self.rownames = [os.path.basename(f) for f in self.file_list]
        self.colnames = [line.split()[0] for line in data]

        rows = []
        for f in self.file_list:
            rows.append([line.split()[1] for line in data])
        self.data = rows


class TestExampleAggregator(BaseTestCase):

    '''Test that the example is working as expected'''

    def test_aggregate_sets_rownames(self):
        file_list = list('akdnehs')
        obj = ExampleAggregator(file_list)
        obj.aggregate()

        expected_rownames = sorted(file_list)
        self.assertEqual(obj.rownames, expected_rownames)

    def test_aggregate_sets_colnames(self):
        file_list = list('akdnehs')
        obj = ExampleAggregator(file_list)
        obj.aggregate()

        expected_data = [row.split()[0]
                         for row in sample_aggregate.split('\n')]
        self.assertEqual(obj.colnames, expected_data)

    def test_aggregate_sets_data(self):
        file_list = list('akdnehs')
        obj = ExampleAggregator(file_list)
        obj.aggregate()

        expected_data = [[
            row.split()[1]
            for row in sample_aggregate.split('\n')
        ], ] * len(file_list)
        self.assertEqual(obj.data, expected_data)


class TestAggregateBase(BaseTestCase):

    def mock_CsvReader(self):
        patcher = patch('libpipe.aggregate.base.CsvReader')
        m = patcher.start()
        self.addCleanup(patcher.stop)

        m.keys = lambda: [line.split()[0] for line in sample_aggregate]
        self.csv = m

    def setup_aggregator(self):
        file_list = list('abcd')
        obj = ExampleAggregator(file_list=file_list)
        obj.aggregate()
        return obj

    # misc

    def test_cannot_init_Aggregator(self):
        with self.assertRaises(TypeError):
            Aggregator()  # abstract methods not implemented!

    # template

    def test_init_loads_valid_template_path(self):
        '''Test that the template path is set and a valid file

        NOTE: Do NOT mock isfile--the default template should exist!
        '''

        obj = ExampleAggregator()
        self.assertTrue(
            os.path.isfile(obj.template_path),
            'No template file found'
        )

    def test_init_given_template_file_changes_instance_only(self):
        '''Test the class template is not modified'''

        orig = ExampleAggregator.template
        obj = ExampleAggregator(template='foo/bar.html')
        self.assertEqual(obj.template_path, 'foo/bar.html')
        self.assertEqual(Aggregator.template, orig)

    # file list

    def test_init_deep_copies_and_sorts_given_files(self):
        file_list = list('akdnehs')
        obj = ExampleAggregator(file_list)

        self.assertEqual(obj.file_list, sorted(file_list), 'Files not sorted')
        file_list[3] = 'foo'
        self.assertNotIn('foo', obj.file_list, 'File list not deep copied')

    def test_init_deep_copies_given_files_without_sort(self):
        file_list = list('akdnehs')
        obj = ExampleAggregator(file_list, sort=False)

        self.assertEqual(obj.file_list, file_list, 'Files are sorted')
        file_list[3] = 'foo'
        self.assertNotIn('foo', obj.file_list, 'File list not deep copied')

    # aggregate

    @unittest.skip('This should be determined by children!')
    def test_aggregate_calls_open_for_each_file_by_default(self):
        m = self.setup_mock_read(sample_aggregate)

        file_list = list('abcd')
        obj = ExampleAggregator(file_list=file_list)
        obj.aggregate()

        self.assertEqual(m.call_count, len(file_list), 'Bad call count')

    # write_html

    def test_write_html_reads_from_set_template(self):
        m = self.setup_mock_read('template data')
        obj = self.setup_aggregator()

        try:
            obj.write_html()
        except ValueError:
            pass  # bad template error--expected
        m.assert_any_call(obj.template_path, 'r')

    def test_write_html_creates_filename_from_commonprefix(self):
        m = self.setup_mock_write()
        obj = self.setup_aggregator()

        obj.file_list = ['foo/bar_' + f for f in obj.file_list]

        with patch.object(obj, '_parse_template'):
            obj.write_html()

        m.assert_any_call('foo/bar_aggregate.html', 'w')

    def test_write_html_raises_ValueError_if_no_template_found(self):
        self.setup_mock_read('template data')

        obj = self.setup_aggregator()

        with self.assertRaises(ValueError):
            obj.write_html()

    # template loading

    def test_write_html_sets_base_html(self):
        obj = self.setup_aggregator()

        obj.write_html()
        self.assertIn('<html>', obj.base_html)

    def test_write_html_sets_template_html(self):
        obj = self.setup_aggregator()

        obj.write_html()
        self.assertIn('</tr>', obj.template_html)

    def test_write_html_sets_join_html(self):
        obj = self.setup_aggregator()

        obj.write_html()
        self.assertIn('</td>', obj.join_html['data'])

    def test_write_html_buffers_colnames_with_empty_cell(self):
        obj = self.setup_aggregator()

        expected_data = '<th></th><td>' + '</td><td>'.join(obj.colnames)
        with io.StringIO() as h:
            obj.write_html(html_file=h)
            written = h.getvalue()

        written = re.sub(r'\s', ' ', written)  # shrink all whitespace
        written = re.sub(r'>\s*<', '><', written)  # remove btwn tag spaces
        self.assertIn(expected_data, written)
