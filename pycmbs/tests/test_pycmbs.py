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
    

if __name__ == "__main__":
    unittest.main()
