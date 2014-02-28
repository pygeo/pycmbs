# -*- coding: UTF-8 -*-
import unittest
from nose.tools import assert_raises
from pycmbs.benchmarking.report import Report
import os
import numpy as np
import matplotlib.pyplot as plt


class TestData(unittest.TestCase):

    def setUp(self):
        if not os.path.exists('./tmp'):
            os.makedirs('./tmp')
        self.R = Report('testfile', 'myreport', 'Alex Loew', outdir='./tmp/')

    def tearDown(self):
        if os.path.exists('./tmp/'):
            os.system('rm -rf ./tmp/')

    def test_ReportInit(self):
        self.assertEqual(self.R.filename, './tmp/testfile.tex')
        self.assertEqual(self.R.format, 'png')
        self.assertEqual(self.R.author, 'Alex Loew')

    def test_open_report(self):
        self.R.open()
        self.assertTrue(os.path.exists(self.R.filename))

    def test_open_report_MissingDirectory(self):
        self.R = Report('testfile', 'myreport', 'Alex Loew', outdir='./tmp/nixdir')
        self.R.open()
        self.assertTrue(os.path.exists(self.R.filename))
        self.assertTrue(os.path.exists('./tmp/nixdir'))

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

    def test_report_InvalidFigure(self):
        f = None
        r = self.R.figure(f, caption='My figure caption')
        self.assertEqual(r, None)

    #~ def test_report_CaptureFigures(self):
        #~ if os.path.exists(self.R.outdir + 'fig_00001.png'):
            #~ os.remove(self.R.outdir + 'fig_00001.png')
        #~ if os.path.exists(self.R.outdir + 'fig_00002.png'):
            #~ os.remove(self.R.outdir + 'fig_00002.png')
        #~ if os.path.exists(self.R.outdir + 'fig_00003.png'):
            #~ os.remove(self.R.outdir + 'fig_00003.png')
#~
        #~ f1 = plt.figure()
        #~ ax1 = f1.add_subplot(111)
        #~ ax1.plot(np.random.random(100))
        #~ f2 = plt.figure()
        #~ ax2 = f2.add_subplot(111)
        #~ ax2.plot(np.random.random(200))
        #~ f3 = plt.figure()
        #~ ax3 = f3.add_subplot(111)
        #~ ax3.plot(np.random.random(300))
        #~
        #~ self.R.capture_figures()
#~
        #~ self.assertTrue(os.path.exists(self.R.outdir + 'fig_00001.png'))
        #~ self.assertTrue(os.path.exists(self.R.outdir + 'fig_00002.png'))
        #~ self.assertTrue(os.path.exists(self.R.outdir + 'fig_00003.png'))

    def test_input(self):
        self.R.input('testname')




