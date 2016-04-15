'''Test the CmdBase class ability to use custom types'''

import copy
from unittest import mock

from libpipe.type import file as _file
from libpipe.type import index

from tests.cmd import test_base


class TestBaseCmd_CustomTypes(test_base.BaseTestCase):

    '''Ensure type registry is not affected'''

    @classmethod
    def setUpClass(cls):
        cls.old_registry = copy.deepcopy(cls.TYPE.registry)
        cls.old_children = copy.deepcopy(cls.TYPE.children)

    @classmethod
    def tearDownClass(cls):
        # Reload saved registries after each test
        cls.TYPE.registry = copy.deepcopy(cls.old_registry)
        cls.TYPE.children = copy.deepcopy(cls.old_children)

    def tearDown(self):
        # Reload saved registries after each test
        self.TYPE.registry = copy.deepcopy(self.old_registry)
        self.TYPE.children = copy.deepcopy(self.old_children)


class TestBaseCmd_FileType(TestBaseCmd_CustomTypes):

    TYPE = _file.FileType

    def test_match_file_type(self):
        FType = _file.factory(
            name='Cmd_FileType_Interaction',
            extns=['.hello', '.my', '.baby'],
        )
        self.CMD.attr.req_types = [
            [('-U', ), (FType, )],
        ]
        cmd = self.CMD()
        cmd.input = lambda: ['oh.my', 'somthin.else']

        cmd._match_input_with_args()
        self.assertEqual(cmd.kwargs['-U'], 'oh.my')


class TestBaseCmd_IndexType(TestBaseCmd_CustomTypes):

    TYPE = index.IndexType

    def test_match_index(self):
        '''Test that a custom type is matched appropriately'''
        # TODO(sjbush): replace index type with the base type

        SampleIndex = index.factory(
            name='SampleIndex', extns=['.bt2'], counts=[2])
        # log.debug(type(SampleIndex))
        # log.debug(isinstance(SampleIndex, cls.TYPE))
        self.CMD.attr.req_types = [
            [('-x', ), (SampleIndex, )],
        ]
        cmd = self.CMD()
        cmd.input = lambda: ['file1.txt', 'path/to/idx']

        # CRITICAL: Mocking 'walk_file' does NOT restrict checking
        #   based on the given argument. This will cause issues as
        #   any str type WILL be identified as an index.
        #   Therefore, we mock 'walk_safe' to allow selection based
        #   on the given value. Unfortunately, this is partially testing
        #   walk_file.
        with mock.patch('libpipe.util.path.walk_safe', ) as mock_walk:
            mock_walk.return_value = {
                'file': ['idx.1.bt2', 'idx.2.bt2', 'idx.txt']}
            cmd._match_input_with_args()
        self.assertEqual(cmd.kwargs['-x'], 'path/to/idx')
        self.assertIsInstance(cmd.kwargs['-x'], SampleIndex)

    def test_match_converts_index_path_str_to_index_type(self):
        SampleIndex = index.factory(
            name='SampleIndex', extns=['.bt2'], counts=[2])
        self.CMD.attr.req_types = [
            [('-x', ), (SampleIndex, )],
        ]
        cmd = self.CMD()
        cmd.input = lambda: ['file1.txt', 'path/to/idx']

        with mock.patch('libpipe.util.path.walk_safe', ) as mock_walk:
            mock_walk.return_value = {
                'file': ['idx.1.bt2', 'idx.2.bt2', 'idx.txt']}
            cmd._match_input_with_args()
        self.assertTrue(
            issubclass(cmd.kwargs['-x'].__class__, SampleIndex),
            "Given path not converted to IndexType",
        )
        self.assertTrue(hasattr(cmd.kwargs['-x'], 'extns'))
