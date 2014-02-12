# -*- coding: utf-8 -*-


import unittest
import os
import pycmbs.benchmarking.utils as utils


class TestUtils(unittest.TestCase):

    def test_get_data_pool_OK(self):
        os.environ['SEP'] = os.getcwd()
        p = utils.get_data_pool_directory()
        self.assertEqual(p, os.getcwd())

    #~ def test_get_data_pool_Default(self):
        #~ if 'SEP' in os.environ.keys():
            #~ xx = os.environ.pop('SEP')
#~
        #~ p = utils.get_data_pool_directory()
        #~ self.assertEqual(p, '/pool/SEP')

    @unittest.skip('TODO: acquire files using wget script!')
    def test_get_generic_landseamask_DEFAULT(self):
        #XXX TODO analyze also correctness of results
        ls = utils.get_generic_landseamask(False)
        ls = utils.get_generic_landseamask(True, area='ocean')
        ls = utils.get_generic_landseamask(True, area='global')
        ls = utils.get_generic_landseamask(False, mask_antarctica=False)
        fname = ['land_sea_fractions_remapnn_t63grid.nc', 'land_sea_fractions_remapnn_t63grid_cell_area.nc']
        self.assertTrue(os.path.exists(fname[0]))
        for f in fname:
            if os.path.exists(f):
                os.remove(f)


if __name__ == "__main__":
    unittest.main()
