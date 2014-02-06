# -*- coding: utf-8 -*-


import unittest
from pycmbs.icon import Icon

class TestPycmbsIcon(unittest.TestCase):

    def setUp(self):
        # requires local installation of ICON sample files!
        self.gridfile ='../..//example_data/icon/r2b4_amip.nc'
        self.datafile = '../../example_data/icon/rms0006_atm_phy_DOM01_ML_0001.nc'

    def test_DummyTest(self):
        pass

    def test_IconInit(self):
        x = Icon(None, None, 'None')

    def test_IconInitMissingFile(self):
        x = Icon('no.nc', 'nothing.nc', 'novar')
        with self.assertRaises(ValueError):
            x.read()

    def test_IconInitMissingGridFile(self):
        x = Icon(self.datafile, 'nothing.nc', 'novar')
        with self.assertRaises(ValueError):
            x.read()


    #~ def test_IconReadOK(self):
        #~ x = Icon(self.datafile, self.datafile, 'rsns')
        #~ x.read()

if __name__ == "__main__":
    unittest.main()
