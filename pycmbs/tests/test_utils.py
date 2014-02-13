# -*- coding: utf-8 -*-


import unittest
import os
import pycmbs.benchmarking.utils as utils


class TestUtils(unittest.TestCase):

    def setUp(self):
        pass

    def test_get_data_pool_OK(self):
        os.environ['SEP'] = os.getcwd()
        p = utils.get_data_pool_directory()
        self.assertEqual(p, os.getcwd())

    def test_get_data_pool_Default(self):
        if 'SEP' in os.environ.keys():
            xx = os.environ.pop('SEP')
        p = utils.get_data_pool_directory()
        self.assertEqual(p, '/pool/SEP/')
        os.environ.update({'SEP': xx})

    def test_get_data_pool_MultipleDirs(self):
        r = os.environ['SEP']
        os.environ.update({'SEP' : '/some/first/path:' + r })
        p = utils.get_data_pool_directory()
        self.assertEqual(p, r)
        os.environ.update({'SEP': r})

    def test_get_data_pool_MultipleDirs2(self):
        r = os.environ['SEP']
        os.environ.update({'SEP' : r + ':/some/second/dir' })
        p = utils.get_data_pool_directory()
        self.assertEqual(p, r)
        os.environ.update({'SEP': r})

    def test_cdo_tempdir_fromENV(self):
        d = '/some/directory/path'
        os.environ.update({'CDOTEMPDIR': d})
        r = utils.get_temporary_directory()
        self.assertEqual(r, d + '/')

    def test_cdo_tempdir_fromENV1(self):
        d = '/some/directory/path/'
        os.environ.update({'CDOTEMPDIR': d})
        r = utils.get_temporary_directory()
        self.assertEqual(r, d)

    def test_cdo_tempdir_DefaultNoEnv(self):
        if 'CDOTEMPDIR' in os.environ.keys():
            d = os.environ.pop('CDOTEMPDIR')
        r = utils.get_temporary_directory()
        self.assertEqual(r, './')


    #~ def test_get_generic_landseamask_DEFAULT(self):
        #~ #XXX TODO analyze also correctness of results
        #~ ls = utils.get_generic_landseamask(False)
        #~ ls = utils.get_generic_landseamask(True, area='ocean')
        #~ ls = utils.get_generic_landseamask(True, area='global')
        #~ ls = utils.get_generic_landseamask(False, mask_antarctica=False)
        #~ fname = ['land_sea_fractions_remapnn_t63grid.nc', 'land_sea_fractions_remapnn_t63grid_cell_area.nc']
        #~ self.assertTrue(os.path.exists(fname[0]))
        #~ for f in fname:
            #~ if os.path.exists(f):
                #~ os.remove(f)


if __name__ == "__main__":
    unittest.main()
