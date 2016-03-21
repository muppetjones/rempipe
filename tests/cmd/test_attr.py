'''Test the CmdAttributes class

'''

import unittest

from libpipe.cmd.attr import CmdAttributes

sample_attributes = {
    'name': 'test command',
    'synopsis': 'Use to test the CmdAttributes class',
    'description': 'More info about the CmdAttributes class',
    'invoke': 'test_command',

    'args': [
        (0, 'INPUT', 'Something we need'),
        ('-f', 'FILE', 'Better to be explicit'),
        ('-o', 'FILE', 'Output file'),
        ('-n', int, 'A number'),
        ('--foo', 'FILE', 'verbose arg'),
        ('-v', None, 'A random flag'),
        ('-x', None, 'A random flag'),
    ],
    'defaults': {
        '-n': 5,
    },

    'req_args': 1,
    'req_kwargs': ['-f', ],
    'req_types': [
        [(0, ), ('.txt', '.csv', ), ],
    ],
}


class TestCmdAttributes(unittest.TestCase):

    def setUp(self):
        self.kwargs = dict.copy(sample_attributes)

    def test_init_raises_ValueError_if_missing_required_args(self):
        '''Ensure we can't create the obj w/o the basics'''
        required = ['invoke', 'args']
        for req in required:
            setup_dict = dict.copy(self.kwargs)
            del setup_dict[req]
            with self.subTest(missing_req=req):
                with self.assertRaises(ValueError):
                    CmdAttributes(**setup_dict)

    def test_init_deep_copies_all_attributes(self):
        '''Make sure that changes to the input dict don't affect obj'''
        ca = CmdAttributes(**self.kwargs)
        for k, v in self.kwargs.items():
            with self.subTest(attr=k):
                if isinstance(self.kwargs[k], list):
                    self.kwargs[k].append('new_value')  # mutable!
                else:
                    self.kwargs[k] = 'Something completely different'
                self.assertNotEqual(getattr(ca, k), self.kwargs[k])

    def test_init_ignores_unknown_attributes_if_strict(self):
        '''Ensure strict (default=T) ignores unknown attributes'''
        self.kwargs['new_key'] = 'Should not be set'
        ca = CmdAttributes(**self.kwargs)
        with self.assertRaises(AttributeError):
            ca.new_key

    def test_init_sets_unknown_attributes_if_not_strict(self):
        '''Test that obj has unknown attributes if strict=False'''
        self.kwargs['new_key'] = 'Should not be set'
        ca = CmdAttributes(strict=False, **self.kwargs)
        ca.new_key  # should not raise

    def test_init_raises_ValueError_if_defaults_has_unknown_arg(self):
        '''Ensure default values match known args.

        NOTE: We really want to make sure that positional args cannot be
            default, BUT...this is over restrictive.
        NOTE: CmdAttributes does NOT store actual values!
        '''
        self.kwargs['defaults'] = {'--unknown': 'bad_arg'}
        with self.assertRaises(ValueError):
            _ = CmdAttributes(**self.kwargs)

    def test_init_sets_name_from_invoke_if_not_given(self):
        del self.kwargs['name']
        ca = CmdAttributes(**self.kwargs)
        self.assertEqual(ca.invoke, ca.name)

    @unittest.skip('nyi')
    def test_duplicate_creates_new_object(self):
        ca = CmdAttributes(**self.kwargs)
        self.assertNotEqual(ca, ca.duplicate())

    @unittest.skip('nyi')
    def test_duplicate_deep_copies_object(self):
        ca = CmdAttributes(**self.kwargs)
        dup = ca.duplicate()

        for k, v in ca.__dict__.items():
            with self.subTest(attr=k):
                self.assertEqual(getattr(ca, k), getattr(dup, k))
                setattr(ca, k, 'not equal')
                self.assertNotEqual(getattr(ca, k), getattr(dup, k))

    @unittest.skip('nyi')
    def test_allows_custom_attr_if_not_strict(self):
        self.kwargs['strict'] = False
        self.kwargs['custom_attr'] = 5
        ca = CmdAttributes(**self.kwargs)

        self.assertEqual(ca.custom_attr, 5)

    @unittest.skip('nyi')
    def test_disallows_custom_attr_if_strict(self):
        self.kwargs['strict'] = True
        self.kwargs['custom_attr'] = 5
        ca = CmdAttributes(**self.kwargs)

        with self.assertRaises(AttributeError):
            ca.custom_attr
