import unittest
from nose.tools import assert_raises
from pycmbs.benchmarking.report import Report
import os


class TestData(unittest.TestCase):

    def setUp(self):
        self.R = Report('testfile', 'myreport', 'Alex Loew', outdir='/tmp/')

    def test_ReportInit(self):
        self.assertEqual(self.R.filename, '/tmp/testfile.tex')
        self.assertEqual(self.R.format, 'png')
        self.assertEqual(self.R.author, 'Alex Loew')

    def test_open_report(self):
        self.R.open()
        self.assertTrue(os.path.exists(self.R.filename))

    def test_report_features(self):
        self.R._write_separator()
        self.R._write_header()
