import unittest

from libpipe.cmds.help import HelpCmd


class TestHelpCmd(unittest.TestCase):

    def test_sections_default_to_None(self):
        hc = HelpCmd()

        for k, v in hc.dict.items():
            with self.subTest(attr=k):
                self.assertIsNone(v)

    def test_str_fills_template_with_non_None_sections(self):
        sect = {'name': 'my name', 'synopsis': 'my synopsis'}
        hc = HelpCmd(**sect)

        help_txt = str(hc)  # implicit string conversion test
        for k, v in sect.items():
            with self.subTest(attr=k):
                self.assertIn(k, help_txt)

    def test_automatically_creates_section_attributes(self):
        sect = {'name': 'my name', 'synopsis': 'my synopsis'}
        hc = HelpCmd(**sect)

        self.assertEqual(hc.name, sect['name'])

    def test_still_raises_attribute_error_for_non_attributes(self):

        hc = HelpCmd()
        with self.assertRaises(AttributeError):
            hc.missing_attr

    def test_str_returns_template_string(self):
        hc = HelpCmd()

        hc_str = str(hc)
        self.assertIn('NAME', hc_str)  # test for expected 'NAME' in template
        self.assertNotIn('name', hc_str)  # ignore case test
