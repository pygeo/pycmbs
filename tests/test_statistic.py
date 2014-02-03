import unittest

class TestStatistic(unittest.TestCase):
    def test_get_significance(self):
        from pycmbs.statistic import get_significance
        self.assertAlmostEqual(get_significance(0.5,100.),1.18049202704e-07,delta=1.e-6)

if __name__ == "__main__":
    unittest.main()
