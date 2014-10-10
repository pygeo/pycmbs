# -*- coding: utf-8 -*-
"""
This file is part of pyCMBS. (c) 2012-2014
For COPYING and LICENSE details, please refer to the file
COPYRIGHT.md
"""

from pycmbs.benchmarking.external_analysis import ExternalAnalysis

from unittest import TestCase
import unittest

import os
import numpy as np
import datetime

from nose.tools import assert_raises

class TestData(unittest.TestCase):
    def setUp(self):
        self.test_file = 'test_template.txt'
        if os.path.exists(self.test_file):
            os.remove(self.test_file)
        o = open(self.test_file, 'w')
        o.write('#This is a testfile\n')
        o.write('echo <mytag>\n')
        o.write('touch <outputfile>')
        o.close()

    def tearDown(self):
        if os.path.exists(self.test_file):
            os.remove(self.test_file)

    def test_external_invalid_executable(self):
        with self.assertRaises(ValueError):
            E = ExternalAnalysis('python ', self.test_file, {'nothing':'first'})

    def test_external_output_dir(self):
        E = ExternalAnalysis('python <INPUTFILE>', self.test_file, {'nothing':'second'}, output_directory='./test/nothing')
        self.assertEqual(E.output_directory, './test/nothing/')
        self.assertTrue(os.path.exists('./test/nothing'))
        os.system('rm -rf ./test')

    def test_external_missing_template(self):
        with self.assertRaises(ValueError):
            E = ExternalAnalysis('python <INPUTFILE>', 'notemplate.txt', {'nothing':'first'}, output_directory='./test/nothing')

    def test_external_create_script(self):
        E = ExternalAnalysis('python <INPUTFILE>', self.test_file, {'mytag':'this_is_the_output', 'outputfile':'pycmbstestfile.txt'}, output_directory='./test/nothing')
        E._create_script()

        self.assertTrue(os.path.exists('./test/nothing/' + self.test_file))
        o = open('./test/nothing/' + self.test_file, 'r')
        l = o.readline()
        self.assertEqual(l, '#This is a testfile\n')
        l = o.readline()
        self.assertEqual(l, 'echo this_is_the_output\n')
        l = o.readline()
        self.assertEqual(l, 'touch pycmbstestfile.txt')

    def test_external_run(self):
        E = ExternalAnalysis('bash <INPUTFILE>', self.test_file, {'mytag':'this_is_the_output', 'outputfile':'pycmbstestfile.txt'}, output_directory='./test/nothing')
        E.run()
        self.assertTrue(os.path.exists('./test/nothing/pycmbstestfile.txt'))
        os.remove('./test/nothing/pycmbstestfile.txt')
        os.system('rm -rf ./test')
