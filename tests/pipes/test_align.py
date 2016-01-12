from unittest.mock import patch

from remsci.lib.tests.base import BaseTestCase

import libpipe
from libpipe.pipes.align import _AlignPipe, AlignPipe, NestedAlignPipe

import logging
log = logging.getLogger(__name__)


class Test_AlignPipe(BaseTestCase):

    def test_output_returns_count_file_only(self):
        self.mock_method(
            libpipe.cmds.align.HisatCmd, '_additional_requirements')

        with patch('os.path.isfile', return_value=True):
            ap = _AlignPipe(input_list=['test.fq'], genome='hsapiens')
            ap.cmd()

        self.assertEqual(ap.output(), ('test_hsapiens.count.bed', ))

    def test_soft_output_includes_same_number_of_fastq_files_as_input(self):
        self.mock_method(
            libpipe.cmds.align.HisatCmd, '_additional_requirements')

        input_list = ['test.fq']
        for i in range(2):
            input_list = ['test' + str(j) + '.fq' for j in range(i + 1)]
            with patch('os.path.isfile', return_value=True):
                ap = _AlignPipe(input_list=input_list, genome='hsapiens')
                ap.cmd()

            ap.soft_output = True
            olist = ap.output()
            iolist = [f for f in olist if f.endswith('.fastq')]
            self.assertEqual(len(iolist), len(input_list))

    def test_raises_ValueError_if_multiple_genomes_given_during_init(self):
        with self.assertRaises(ValueError):
            _AlignPipe(genome=['hsapiens', 'mtuberculosis'])

    def test_accepts_genome_list_with_single_genome_given_during_init(self):
        _AlignPipe(genome=['hsapiens', ])  # no exception!


class TestAlignPipe(BaseTestCase):

    def setUp(self):
        self.kwargs = {
            'input_list': ['file_1.fq', 'file_2.fq'],
            'genome': 'hsapiens',
            'odir': './'
        }

    def test_init_raises_ValueError_if_no_input_list_given(self):
        with self.assertRaises(ValueError):
            AlignPipe(**{'genome': 'hsapiens', })

    def test_init_raises_ValueError_if_no_genome_given(self):
        with self.assertRaises(ValueError):
            AlignPipe(**{'input_list': ['file_1.fq', 'file_2.fq'], })

    def test_fastqc_args_set_to_input_list_during_init(self):
        align = AlignPipe(**self.kwargs)
        trim = align.cmds[0]
        fastqc = trim.cmds[0]
        self.assertEqual(fastqc.args, self.kwargs['input_list'])

    def test_align_input_set_to_trim_output(self):
        self.mock_method(
            libpipe.cmds.align.HisatCmd, '_additional_requirements')
        self.mock_class('os.path.isfile', return_value=True)
        align = AlignPipe(**self.kwargs)
        align.cmd()
        self.assertEqual(align.cmds[0].output(), align.cmds[1].input())

    def test_setup_calls_TrimPipe_and_AlignPipe(self):
        sub1 = self.mock_class('libpipe.pipes.align.TrimPipe')
        sub2 = self.mock_class('libpipe.pipes.align._AlignPipe')

        AlignPipe(**self.kwargs)
        self.assertEqual(sub1.call_count, 1)
        self.assertEqual(sub2.call_count, 1)


class TestNestedAlignPipe(BaseTestCase):

    def setUp(self):
        self.kwargs = {
            'input_list': ['file_1.fq', 'file_2.fq'],
            'genome': ['hsapiens', 'mtuberculosis'],
            'odir': './'
        }

    def test_raises_ValueError_if_single_genome_given_during_init(self):
        self.kwargs['genome'] = 'hsapiens'
        with self.assertRaises(ValueError):
            NestedAlignPipe(**self.kwargs)

    def test_subsequent_align_pipes_use_unaligned_output_from_prev_pipe(self):
        self.mock_method(
            libpipe.cmds.align.HisatCmd, '_additional_requirements')
        self.mock_class('os.path.isfile', return_value=True)

        pipe = NestedAlignPipe(**self.kwargs)
        pipe.cmd()  # force input update

        subpipe = pipe.cmds[-1]
        hisat = subpipe.cmds[0]

        self.assertIn('hsapiens_unal', hisat.kwargs['-1'])
