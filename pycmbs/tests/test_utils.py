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

    #~ def test_get_data_pool_Default(self):
        #~ if 'SEP' in os.environ.keys():
            #~ xx = os.environ.pop('SEP')
#~
        #~ p = utils.get_data_pool_directory()
        #~ self.assertEqual(p, '/pool/SEP')

    def test_get_generic_landseamask_DEFAULT(self):
        #XXX TODO analyze also correctness of results
        ls = utils.get_generic_landseamask(False)
        ls = utils.get_generic_landseamask(True, area='ocean')
        ls = utils.get_generic_landseamask(True, area='global')
        ls = utils.get_generic_landseamask(False, mask_antarctica=False)



if __name__ == "__main__":
    unittest.main()
