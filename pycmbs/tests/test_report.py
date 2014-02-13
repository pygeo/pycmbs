import unittest
from nose.tools import assert_raises
from pycmbs.benchmarking.report import Report
import os
import numpy as np
import matplotlib.pyplot as plt


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
        if os.path.exists(self.R.filename):
            os.remove(self.R.filename)

        f = plt.figure()
        ax = f.add_subplot(111)
        ax.plot(np.random.random(200))

        self.R.open()
        self.R.section('This is a section')
        self.R.subsection('This is a subsection')
        self.R.subsubsection('This is a subsubsection')

        self.R.newpage()
        self.R.clearpage()

        self.R.barrier()
        self.R.open_table()
        self.R.close_table(caption='This is my caption')
        #self.input('filename_to_input.tex')

        self.R.newpage()
        self.R.figure(f, caption='My figure caption')

        self.R.close()
        self.R.compile()

        if os.path.exists(self.R.filename):
            os.remove(self.R.filename)
        if os.path.exists(self.R.filename[:-3]+'.pdf'):
            os.remove(self.R.filename[:-3]+'.pdf')

    def test_input(self):
        self.R.input('testname')




