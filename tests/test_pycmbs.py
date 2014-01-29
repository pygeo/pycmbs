import unittest
import pycmbs

class testPycmbs(unittest.TestCase):
    def setUp(self):
        pass

    def testPycmbsName(self):
        import pycmbs as module
        ref_name = 'pycmbs'
        test_name = module.__name__
        self.assertEqual(ref_name, test_name)
    
    def testPycmbsConstants(self):
        test_earth_radius = pycmbs.EarthRadius
        ref_earth_radius = 6371000. 
        self.assertEqual(test_earth_radius, ref_earth_radius)



if __name__ == "__main__":
    unittest.main()
