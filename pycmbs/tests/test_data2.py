#
# TESTS with some test data
#
"""
This file is part of pyCMBS. (c) 2012-2014
For COPYING and LICENSE details, please refer to the files
LICENSE.md and COPYRIGHT.md
"""


from unittest import TestCase
import unittest

from pycmbs.data import Data

import os
import pylab as pl
import numpy as np

from nose.tools import assert_raises


class TestData(unittest.TestCase):
    def setUp(self):
        pass
        #~ self.file = 'air.mon.mean.nc' # testfilename
        #~ if not os.path.exists(self.file):
            #~ self._download_sample()
        #~ if os.path.exists(self.file):
            #~ # check if download was sucessfull
            #~ self.D = Data(self.file, 'air', read=True)
            #~ self.areafile = self.file[:-3] + '_cell_area.nc'
        #~ else:
            #~ # download failed for some reason
            #~ self.file = None

    def _download_sample(self):
        if not os.path.exists(self.file):
            try:
                print 'Downloading sample data ... this might take a few minutes (only needed at first run)'
                os.system('wget --ftp-user=anonymous --ftp-password=nix ftp://ftp.cdc.noaa.gov/Datasets/ncep.reanalysis.derived/surface/air.mon.mean.nc')
            except:
                raise ValueError('Can not download sampled data file from NCEP server. Please try manually')

            if not os.path.exists(self.file):
                raise ValueError('Something with the file download went wrong!')

    def test_cell_area(self):
        pass
        #~ if self.file is not None:
            #~ pass
            #~ del self.D.cell_area
            #~ os.remove(self.areafile)
            #~ self.D._set_cell_area()
#~
            #~ # check if cell area file is existing
            #~ self.assertTrue(os.path.exists(self.areafile))
#~
            #~ # now remove the cell area file to force calculation
            #~ os.remove(self.areafile)
            #~ del self.D.cell_area
            #~ self.D._set_cell_area()
            #~ self.assertTrue(os.path.exists(self.areafile))

    def test_cell_area_InvalidLatLon(self):
        pass
        #~ if self.file is not None:
            #~ self.D.lon = None
            #~ self.D.lat = None
            #~ self.D.cell_area = None
#~
            #~ os.remove(self.areafile)
            #~ self.D._set_cell_area()
            #~ self.assertTrue(np.all(self.D.cell_area == 1.))

    def test_zonal_mean(self):
        pass
        #~ if self.file is not None:
            #~ zm = self.D.get_zonal_mean()
            #~ ZM = self.D.get_zonal_mean(return_object=True)
            #~ self.assertTrue(np.all(np.abs(1.-zm / ZM.data.T) < 1.E-6 ))

    def test_get_deseasonalized_anomaly(self):
        pass
        #~ if self.file is not None:
            #~ c = self.D.get_deseasonalized_anomaly(base='current')
            #~ c = self.D.get_deseasonalized_anomaly(base='all')

if __name__ == '__main__':
    unittest.main()



