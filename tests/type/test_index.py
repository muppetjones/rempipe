
# import unittest
# from unittest import mock


# from libpipe.util import path
from libpipe.type import index
from tests.type import test_base

import logging
log = logging.getLogger(__name__)


class TestTypeIndexBase(test_base.TypeBaseTestCase):

    def setUp(self):
        self.TYPE = index.IndexType

    def test_factory_sets_expected_extns_attr_to_None(self):
        index_obj = index.IndexType.factory()
        self.assertIsNone(index_obj.expected_extns)

    def test_factory_sets_expected_counts_attr_to_None(self):
        index_obj = index.IndexType.factory()
        self.assertIsNone(index_obj.expected_counts)

    def test_factory_sets_expected_counts_to_None(self):
        _index = index.IndexType.factory()
        self.assertIsNone(_index.expected_counts)

    def test_factory_list_sets_extns_in_alpha_order(self):
        expected = ['.rev.bt2', '.bt2', ]
        _index = index.IndexType.factory(expected)
        self.assertEqual(_index.expected_extns, sorted(expected))

    def test_factory_multi_sets_extns_in_alpha_order(self):
        expected = ['.rev.bt2', '.bt2', ]
        _index = index.IndexType.factory(*expected)
        self.assertEqual(_index.expected_extns, sorted(expected))

    def test_factory_sets_expected_attr_from_kwargs(self):
        extn_dict = {'.bt2': 4, '.rev.bt2': 2, '.foo': None, '.bar': None}

        # repeat check multiple times to ensure consistent
        for i in range(50):
            _index = index.IndexType.factory(**extn_dict)

            # extns
            expected = sorted(list(extn_dict.keys()))
            self.assertEqual(_index.expected_extns, expected)

            # counts
            expected = [extn_dict[k] for k in sorted(list(extn_dict.keys()))]
            self.assertEqual(_index.expected_counts, expected)

    def test_factory_dict_sets_extns_in_alpha_order(self):
        extn_dict = {'.bt2': 4, '.rev.bt2': 2, '.foo': None, '.bar': None}

        # repeat check multiple times to ensure consistent
        for i in range(50):
            _index = index.IndexType.factory(extn_dict)
            expected = sorted(list(extn_dict.keys()))
            self.assertEqual(_index.expected_extns, expected)

    def test_factory_dict_sets_counts_in_extn_order(self):
        extn_dict = {'.bt2': 4, '.rev.bt2': 2, '.foo': None, '.bar': None}

        # repeat check multiple times to ensure consistent
        for i in range(50):
            _index = index.IndexType.factory(extn_dict)
            expected = [extn_dict[k] for k in sorted(list(extn_dict.keys()))]
            self.assertEqual(_index.expected_counts, expected)
