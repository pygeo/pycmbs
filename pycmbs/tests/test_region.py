# -*- coding: utf-8 -*-

import unittest
from pycmbs import region
import tempfile
import os
import shapefile
from pycmbs.region import RegionShape

class TestRegion(unittest.TestCase):

    def _generate_line_sample_shapefile(self):
        w = shapefile.Writer(shapefile.POLYLINE)
        #~ w.line(parts=[[[1,5],[5,5],[5,1],[3,3],[1,1]]])
        w.poly(parts=[[[1,3],[5,3],[3,5], [1,3]   ]], shapeType=shapefile.POLYLINE)
        w.field('FIRST_FLD','C','40')
        w.field('SECOND_FLD','C','40')
        w.record('First','Line')
        w.record('Second','Line')
        tfile = tempfile.mkdtemp() + os.sep + 'mypolylinefile'
        w.save(tfile)
        return tfile

    def _generate_poly_sample_shapefile(self):
        # https://code.google.com/p/pyshp/wiki/PyShpDocs
        w = shapefile.Writer(shapefile.POLYGON)
        w.poly(parts=[[[1,5],[5,5],[5,1],[3,3],[1,1]]])
        w.field('FIRST_FLD','C','40')
        w.field('SECOND_FLD','C','40')
        w.record('First','Polygon')
        tfile = tempfile.mkdtemp() + os.sep + 'mypolyfile'
        w.save(tfile)
        return tfile

    def _generate_point_sample_shapefile(self):
        # example taken from pyshp examples: https://code.google.com/p/pyshp/

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

        # Create a polygon shapefile
        w = shapefile.Writer(shapefile.POLYGON)
        w.poly(parts=[[[1,5],[5,5],[5,1],[3,3],[1,1]]])
        w.field('FIRST_FLD','C','40')
        w.field('SECOND_FLD','C','40')
        w.record('First','Polygon')
        tfile = tempfile.mkdtemp() + os.sep + 'mytestshape'
        w.save(tfile)
        return tfile

    def setUp(self):
        pass

    def test_shape_write(self):
        shp_file = self._generate_point_sample_shapefile()
        # try to read some information
        sf = shapefile.Reader(shp_file)
        shapes = sf.shapes()
        print shapes[0].bbox

    def test_RegionShape(self):
        # test to read a shapefile
        shp_file = self._generate_poly_sample_shapefile()
        RS = RegionShape(shp_file)
        for k in RS.regions.keys():
            r = RS.regions[k]
            print r.lon
            print r.lat





if __name__ == "__main__":
    unittest.main()

