# -*- coding: utf-8 -*-

import unittest
from pycmbs import region
import tempfile
import os

class TestRegion(unittest.TestCase):

    def setUp(self):
        pass

    def test_DummyTest(self):
        pass

    def test_shape_write(self):
        # example taken from pyshp examples: https://code.google.com/p/pyshp/
        import shapefile

        # Make a point shapefile
        w = shapefile.Writer(shapefile.POINT)
        w.point(90.3, 30)
        w.point(92, 40)
        w.point(-122.4, 30)
        w.point(-90, 35.1)
        w.field('FIRST_FLD')
        w.field('SECOND_FLD','C','40')
        w.record('First','Point')
        w.record('Second','Point')
        w.record('Third','Point')
        w.record('Fourth','Point')
        w.save('shapefiles/test/point')
        # Create a polygon shapefile
        w = shapefile.Writer(shapefile.POLYGON)
        w.poly(parts=[[[1,5],[5,5],[5,1],[3,3],[1,1]]])
        w.field('FIRST_FLD','C','40')
        w.field('SECOND_FLD','C','40')
        w.record('First','Polygon')
        tfile = tempfile.mkdtemp() + os.sep + 'mytestshape'
        w.save(tfile)

        # try to read some information
        sf = shapefile.Reader(tfile)
        shapes = sf.shapes()
        print shapes[0].bbox






if __name__ == "__main__":
    unittest.main()

