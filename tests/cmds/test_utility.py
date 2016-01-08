import os.path
import unittest
from unittest.mock import patch

import libpipe
from libpipe.cmds.utility import FastqcCmd

import logging
log = logging.getLogger(__name__)


class TestFastqcCmd(unittest.TestCase):

    def test_output_list_starts_with_input(self):

        args = [v + '.fq' for v in list('abcd')]
        cmd = FastqcCmd(*args)
        output = cmd.output()

        self.assertEqual(output[:len(args)], args)

    def test_output_list_includes_html_for_each_input(self):

        args = [v + '.fq' for v in list('abcd')]
        cmd = FastqcCmd(*args)
        output_html = [o for o in cmd.output() if o.endswith('.html')]

        # REMEMBER: fastqc adds it's own name!
        expected = [a.replace('.fq', '_fastqc.html') for a in args]
        self.assertEqual(output_html, expected)

    def test_output_list_includes_zip_for_each_input(self):

        args = [v + '.fq' for v in list('abcd')]
        cmd = FastqcCmd(*args)
        output = [o for o in cmd.output() if o.endswith('.zip')]

        # ensure fastqc will not try to extract the data
        try:
            cmd.flags.remove(['--extract'])
        except ValueError:
            pass  # already gone, no problem

        # REMEMBER: fastqc adds it's own name!
        expected = [a.replace('.fq', '_fastqc.zip') for a in args]
        self.assertEqual(output, expected)

    def test_output_list_does_not_include_zip_if_extract_flag_set(self):

        args = [v + '.fq' for v in list('abcd')]
        cmd = FastqcCmd(*args)
        cmd.flags.append('--extract')

        # REMEMBER: fastqc adds it's own name!
        output = [o for o in cmd.output() if o.endswith('.zip')]
        self.assertEqual(output, [])

    def test_output_noninput_items_use_input_directory(self):
        args = ['~/tmp/' + v + '.fq' for v in list('abcd')]
        cmd = FastqcCmd(*args)
        output_html = [o for o in cmd.output() if o.endswith('.html')]

        self.assertEqual(
            [os.path.dirname(o) for o in output_html],
            ['~/tmp'] * len(args)
        )

    def test_output_noninput_items_use_o_flag_directory_when_given(self):
        args = ['~/tmp/' + v + '.fq' for v in list('abcd')]
        kwargs = {'-o': '~/work'}
        cmd = FastqcCmd(*args, **kwargs)
        output_html = [o for o in cmd.output() if o.endswith('.html')]

        self.assertEqual(
            [os.path.dirname(o) for o in output_html],
            [kwargs['-o']] * len(args)
        )
