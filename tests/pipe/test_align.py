
import os.path
from unittest import mock

from libpipe.pipe import align
from tests.base import LibpipeTestCase  # includes read and write mock


import logging
log = logging.getLogger(__name__)

DEFAULT_INPUT = ['data/seq_1.fq', 'data/seq_2.fq', 'genomes/GRCh38p5']


class TestAlignPipe_MockedCmds(LibpipeTestCase):

    def setUp(self):

        self.used_cmds = [
            'align.Hisat2Cmd',
            'samtools.SamtoolsSortCmd',
            'samtools.SamtoolsIndexCmd',
            'count.HtseqCountCmd',
        ]

        # create mock tuples [mod.NameCmd] -> (name, cmd)
        patchers = [
            mock.patch('libpipe.cmd.{}'.format(cmd), mock.MagicMock())
            for cmd in self.used_cmds
        ]
        self.mock_cmds = [patcher.start() for patcher in patchers]
        for patcher in patchers:
            self.addCleanup(patcher.stop)

    def test_init_creates_default_cmd_objects(self):
        '''Ensure hisat2, samtools sort & index, & htseq-count cmds created'''

        align.AlignPipe()
        for name, mock_cmd in zip(self.used_cmds, self.mock_cmds):
            with self.subTest(cmd=name):
                self.assertEqual(
                    mock_cmd.call_count, 1, '{} not called'.format(name))

    def test_setup_adds_all_cmd(self):
        pipe = align.AlignPipe()

        # must check each instead of call to pipe.add b/c of HACK
        cmds = [cmd() for cmd in self.mock_cmds]
        for i, cmd in enumerate(cmds[:-1]):
            cmd.link.assert_called_once_with(cmds[i + 1])
        self.assertEqual(pipe.cmds, cmds)
        # mock_add.assert_called_once_with(
        #     *[cmd() for name, cmd in self.mock_cmds])


class TestAlignPipe_Linking(LibpipeTestCase):

    def setUp(self):

        # override walk_file method
        # -- used when checking for index
        patcher = mock.patch('libpipe.util.path.walk_file')
        self.mock_walk = patcher.start()
        self.mock_walk.return_value = ['x.{}.ht2'.format(i) for i in range(8)]
        self.addCleanup(patcher.stop)

        self.input = DEFAULT_INPUT
        self.pipe = align.AlignPipe(input=self.input)
        self.pipe.cmd()  # cause matching and checking

        self.cmds = {
            k: v for k, v in zip(
                ['align', 'sort', 'index', 'count'], self.pipe.cmds)
        }

    def test_hisat2_input_matched_with_initial_fastq(self):
        check_keys = ['-1', '-2', '-x']
        found = {
            k: v for k, v in self.cmds['align'].kwargs.items()
            if k in check_keys
        }
        expected = {k: v for k, v in zip(['-1', '-2', '-x'], self.input)}
        self.assertEqual(found, expected)

    def test_samtoolssort_input_matched_with_hisat2_sam(self):
        prefix = os.path.commonprefix(self.input[:2]).rstrip('_')
        self.assertEqual(self.cmds['sort'].args, [prefix + '.sam'])

    def test_samtoolsindex_input_matched_with_sort_bam(self):
        prefix = os.path.commonprefix(self.input[:2]).rstrip('_')
        self.assertEqual(self.cmds['index'].args, [prefix + '.s.bam'])

    def test_htseqcount_input_matched_with_index_bam_and_genome(self):
        prefix = os.path.commonprefix(self.input[:2]).rstrip('_')
        expected = [prefix + '.s.bam', self.input[-1] + '.gtf']
        self.assertEqual(self.cmds['count'].args, expected)

    def test_output_contains_count_file(self):
        pipe = align.AlignPipe(input=DEFAULT_INPUT)
        pipe.cmd()
        prefix = os.path.commonprefix(self.input[:2]).rstrip('_')
        self.assertIn(prefix + '.s.count', pipe.output())

    def test_output_contains_unaligned_seq_fq(self):
        '''Test hisat2 unal (default) returned in output'''
        pipe = align.AlignPipe(input=DEFAULT_INPUT)
        pipe.cmd()
        prefix = os.path.commonprefix(self.input[:2]).rstrip('_')
        expected = [
            prefix + suffix
            for suffix in ['.s.count', '_unal.1.fastq', '_unal.2.fastq']
        ]
        self.assertEqual(pipe.output(), expected)
