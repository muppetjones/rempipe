
import abc
import os.path
import re

from remsci.lib.decorators import file_or_handle

import libpipe.templates


import logging
log = logging.getLogger(__name__)


class Aggregator(metaclass=abc.ABCMeta):

    '''Aggregate data from file series into an HTML table


    '''

    template = 'aggregator.html'
    name = None

    def __init__(self, file_list=[], template='', sort=True):

        if template:
            self.template_path = template
        else:
            self.template_path = os.path.join(
                os.path.dirname(libpipe.templates.__file__),
                self.template,
            )

        if sort:
            self.file_list = sorted(file_list)
        else:
            self.file_list = file_list[:]

        # expected template parameters
        # NOTE: these are reset each write_template
        self.base_html = ''
        self.template_html = ''
        self.join_html = {'data': '', 'header': ''}

    @abc.abstractmethod
    def aggregate(self):
        raise NotImplementedError('abstract method')

    def write_html(self, html_file=None, template=None):

        if not template:
            template = self.template_path

        if not html_file:
            common = os.path.commonprefix(self.file_list)
            common = re.sub('[^a-zA-Z0-9]*$', '', common)  # rm trailing punc
            filename = '_'.join(filter(None, [
                os.path.basename(common), self.name, 'aggregate.html']
            ))
            html_file = os.path.join(
                os.path.dirname(common),
                filename,
            )

        # sets base_html, template_html, and join_html
        try:
            self._parse_template(template)
        except ValueError:
            raise

        # write to file
        self._write(html_file)

    @file_or_handle(mode='r')
    def _parse_template(self, fh):
        self.__compile_regex()

        base_html = fh.read()

        # parse out the template
        m = self.rx['extract_template'].search(
            base_html,
        )

        if not m:
            raise ValueError('Could not find template')

        prefix_html = m.group('prefix')
        template = m.group('template')
        suffix_html = m.group('suffix')

        # parse the relevant tag names
        matches = self.rx['extract_tag'].findall(template)
        joins = {k: '</{0}><{0}>'.format(v) for k, v in matches}

        # remove script lines
        template_str = '\n'.join(template.split('\n')[1:-1])

        self.base_html = prefix_html + '__template__' + suffix_html
        self.template_html = template_str
        self.join_html = joins

        return

    @file_or_handle(mode='w')
    def _write(self, fh):

        # write header
        # -- include extra blanks if rowname is a list
        header_row = self.template_html.replace('{{ rowname }}', '')
        header_row = header_row.replace(
            '{{ data }}',
            self.join_html['data'].join(self.colnames)
        )
        fh.write(header_row)

    # regular expressions

    def __compile_regex(self):
        if 'rx' not in self.__dict__:
            self.rx = {}
            self.__compile_extract_template_regex()
            self.__compile_extract_tag_regex()

    def __compile_extract_tag_regex(self):
        self.rx['extract_tag'] = re.compile(
            r'(?P<name>\w+)-tag=[\'\"](?P<tag>.+?)[\'\"]',
        )

    def __compile_extract_template_regex(self):
        self.rx['extract_template'] = re.compile(
            r'(?P<prefix>.*)' +
            r'(?P<template><script type="text/.*?template.*?".+?</script>)' +
            r'(?P<suffix>.*)',
            re.DOTALL,
        )
