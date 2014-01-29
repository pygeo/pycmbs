import unittest

class testPycmbs(unittest.TestCase):
    def setUp(self):
        pass

    def testImportPycmbs(self):
        import pycmbs as module
        ref_name = 'pycmbs'
        test_name = module.__name__
        self.assertEqual(ref_name, test_name)

if __name__ == "__main__":
    unittest.main()
