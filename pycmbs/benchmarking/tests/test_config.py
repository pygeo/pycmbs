# -*- coding: utf-8 -*-
"""
This file is part of pyCMBS.
(c) 2012- Alexander Loew
For COPYING and LICENSE details, please refer to the LICENSE file
"""

import unittest
from pycmbs.benchmarking import config
import tempfile
import os

class TestPycmbsBenchmarkingConfig(unittest.TestCase):

    def setUp(self):
        self.odir = tempfile.mkdtemp()
        self.test_cfg = self.odir + os.sep + 'test.cfg'
        self.temp_dir = tempfile.mkdtemp()

    def tearDown(self):
        pass

    def test_read_write(self):
        # write some dummy CFG file and try to read it again
        models = []
        models.append({'id' : 'MPI-ESM-LR', 'type' : 'CMIP5', 'experiment' : 'amip', 'path' : 'testpath'})
        models.append({'id' : 'MPI-ESM-MR', 'type' : 'CMIP5', 'experiment' : 'amip', 'path' : 'testpath'})

        CW = config.CFGWriter(self.test_cfg)
        CW.save(temp_dir=self.temp_dir, vars='default', start_date='YYYY-MM-DD',
             stop_date='YYYY-MM-DD', models=models)

        CR = config.ConfigFile(self.test_cfg)
        self.assertTrue(os.path.exists(self.test_cfg))



if __name__ == "__main__":
    unittest.main()

