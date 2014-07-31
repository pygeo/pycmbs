# -*- coding: utf-8 -*-
"""
This file is part of pyCMBS. (c) 2012-2014
For COPYING and LICENSE details, please refer to the file
COPYRIGHT.md
"""

from behave import *
from pycmbs.benchmarking import config2
import tempfile

@given(u'config file is empty')
def step_impl(context):
    tfile = tempfile.TemporaryFile()
    tfile.write('')
    tfile.seek(0)
    context.config_contents = tfile.read()

@then(u'read yaml config file')
def step_impl(context):
    config_contents = context.config_contents
    yaml_config_contents = config2.parse_config(config_contents, fmt='yaml')
    assert yaml_config_contents is None

@then(u'read json config file')
def step_impl(context):
    config_contents = context.config_contents
    json_config_contents = config2.parse_config(config_contents, fmt='json')
    assert json_config_contents is None
