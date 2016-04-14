import os.path

from functools import partial
from unittest import mock

# from libpipe.util import path
from libpipe.type import base as _type
from libpipe.type import file as _file
from tests import base
from tests.type import test_base

import logging
log = logging.getLogger(__name__)


class FileTypeTestCase(base.LibpipeTestCase):

    '''Define default setup and helper functions'''

    def setUp(self):
        super().setUp()
        self.PARENT = _file.FileType
        self.FACTORY = _file.factory

        n_foo, n_bar = (4, 2)
        self.CHILD = self.FACTORY(
            name='TestFile', extns=['.bar', '.foo'], counts=[n_bar, n_foo])


class TestFileMeta(base.LibpipeTestCase):

    '''Tests against the meta directly'''

    def test_FileMeta_inherits_from_TypeMeta(self):
        self.assertTrue(
            issubclass(_file.FileMeta, _type.TypeMeta),
            "FileMeta is NOT subclassed from TypeMeta",
        )


class TestFileType__TypeBase(test_base.TestTypeBase):

    '''Test against FileType directly (inherits from TypeBase tests)'''

    def setUp(self):
        super().setUp()

        # define variables used by TypeBase tests
        # -- note the partial factory declaration!
        self.PARENT = _file.FileType
        self.FACTORY = partial(_file.factory, extns=['.foo'])
        self.CHILD = self.FACTORY('FileFooType')  # includes extns=['.foo']

        # mock out walk_file so we can test check_extns
        patcher = mock.patch('libpipe.util.path.walk_file')
        m = patcher.start()
        self.addCleanup(patcher.stop)
        m.return_value = ['file.foo']
