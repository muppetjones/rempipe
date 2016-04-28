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

    def test_init_raises_TypeError_if_unknown_args(self):
        with self.assertRaises(TypeError):
            CmdAttributes('huh')

    def test_init_accepts_dict_as_single_arg_param(self):
        ca = CmdAttributes(self.kwargs)  # should not raise
        self.assertEqual(ca.invoke, self.kwargs['invoke'])

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
            CmdAttributes(**self.kwargs)

    def test_init_sets_name_from_invoke_if_not_given(self):
        del self.kwargs['name']
        ca = CmdAttributes(**self.kwargs)
        self.assertEqual(ca.invoke, ca.name)

    def test_init_sets_defaults_to_empty_dict_if_not_given(self):
        del self.kwargs['defaults']
        ca = CmdAttributes(**self.kwargs)
        self.assertIsNotNone(ca.defaults)
        self.assertIsInstance(ca.defaults, dict)

    def test_init_defaults_output_args_to_empty_list(self):
        ca = CmdAttributes(**self.kwargs)
        self.assertIsNotNone(ca.output_args, [])
        self.assertIsInstance(ca.output_args, list)

    # def test_init_throws_ValueError_if_multiple_IndexTypes_required(self):
    #     '''Ensure that each cmd can require only one index type
    #
    #     Ensure that a command cannot be initialized if the attribute
    #     requirements list more than one index.
    #     NOTE: This may be limiting, but, e.g., hisat2 should only
    #         expect hisat2 indices.
    #     '''
    #
    #     import libpipe.type.index as index
    #
    #     extns = list('abc')
    #     counts = [2] * len(extns)
    #     i1 = index.factory(extns=extns, counts=counts)
    #     i2 = index.factory(extns=extns, counts=counts)
    #     self.kwargs['req_types'] = [
    #         [('-x', ), (i1, i2)],
    #     ]
    #     with self.assertRaisesRegex(ValueError, 'multiple.*type'):
    #         CmdAttributes(**self.kwargs)
