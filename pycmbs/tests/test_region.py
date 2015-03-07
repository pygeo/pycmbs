# -*- coding: utf-8 -*-
"""
This file is part of pyCMBS.
(c) 2012- Alexander Loew
For COPYING and LICENSE details, please refer to the LICENSE file
"""

import unittest
from pycmbs import region

from pycmbs.region import RegionParser
from pycmbs.region import RegionPolygon
from pycmbs.region import RegionBboxLatLon
from pycmbs.region import RegionIndex
from pycmbs.region import RegionGeneric
from nose.tools import assert_raises
import os
import numpy as np


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
        self.file = tempfile.mkdtemp() + os.sep + 'test.reg'
        o=open(self.file, 'w')
        o.write('[COMMENT]\n')
        o.write('description=IPCC regions as defined in SREX Appendix-3A, Table 3.A.1\n')
        o.write('label=IPCC regions\n')
        o.write('\n')
        o.write('[ALA]\n')
        o.write('coordinates=[(60.000, -105.000),(60.000, -168.022),(72.554, -168.022),(72.554, -105.000)]\n')
        o.write('id=1\n')
        o.write('\n')
        o.write('[AMZ]\n')
        o.write('coordinates=[(-20.000, -66.377),(-1.239, -79.729),(11.439, -68.800),(11.439, -50.000),(-20.000, -50.000)]\n')
        o.write('id=7\n')
        o.write('\n')
        o.write('[CAM]\n')
        o.write('coordinates=[(11.439, -68.800),(-1.239, -79.729),(28.566, -118.323),(28.566, -90.315)]\n')
        o.write('id=10\n')
        o.close()

        # rectangular box file
        self.boxfile = tempfile.mkdtemp() + os.sep + 'test.box'
        o=open(self.boxfile, 'w')
        o.write('# This file specifies a region by rectangular lat/lon coordinates\n')
        o.write('# comments can be always introduced using the hash\n')
        o.write('Tropics,1,-180.,180.,-30.,30.\n')
        o.write('Npolar,2,-180.,180.,63.,90.\n')
        o.write(' # some dummy data\n')  # note that the initial space is by purpose!
        o.write('Spolar,3,-180.,180.,-90.,-63.\n')
        o.close()


    def tearDown(self):
        if os.path.exists(self.file):
            os.remove(self.file)

    def test_regionparser_filenotexisting(self):
        with self.assertRaises(ValueError):
            R = RegionParser('invalid_file.txt')

    def test_regionparser_invalidformat(self):
        with self.assertRaises(ValueError):
            R = RegionParser(self.file, format='noformat')

    def test_regionparser_read_ini(self):
        R = RegionParser(self.file, format='ini')
        self.assertEqual(R.regions['ALA'].lon[0], -105.)
        self.assertEqual(R.regions['ALA'].lon[1], -168.022)
        self.assertEqual(R.regions['ALA'].lon[2], -168.022)
        self.assertEqual(R.regions['ALA'].lon[3], -105.)

        self.assertEqual(R.regions['AMZ'].lat[0], -20.)
        self.assertEqual(R.regions['AMZ'].lat[1], -1.239)
        self.assertEqual(R.regions['AMZ'].lat[2], 11.439)
        self.assertEqual(R.regions['AMZ'].lat[3], 11.439)
        self.assertEqual(R.regions['AMZ'].lat[4], -20.)

    def test_regionparser_read_box(self):
        R = RegionParser(self.boxfile, format='box')
        self.assertTrue(len(R.regions) == 3)
        for k in R.regions.keys():
            self.assertTrue(k in ['Tropics','Spolar','Npolar'])
            self.assertTrue(isinstance(R.regions[k], RegionBboxLatLon))
        self.assertEqual(R.regions['Tropics'].id,1)
        self.assertEqual(R.regions['Npolar'].id,2)
        self.assertEqual(R.regions['Spolar'].id,3)

    def test_RegionPolygon(self):
        lon = np.random.random(10)
        lat = np.random.random(10)
        R = RegionPolygon(111, lon, lat, label='Testregion')
        for i in xrange(len(lon)):
            self.assertEqual(R.lon[i], lon[i])
            self.assertEqual(R.lat[i], lat[i])
        self.assertEqual(111, R.id)

    def test_Region_invalid_label(self):
        lon = np.random.random(10)
        lat = np.random.random(10)
        with self.assertRaises(ValueError):
            R = RegionPolygon(55, lon, lat, label=None)

    def test_Region_Bbox(self):
        R = RegionBboxLatLon(55, -20., 20., -50., 30., label='TestBbox')
        self.assertEqual(R.lon[0], -20.)
        self.assertEqual(R.lon[1], -20.)
        self.assertEqual(R.lon[2], 20.)
        self.assertEqual(R.lon[3], 20.)

        self.assertEqual(R.lat[0], -50.)
        self.assertEqual(R.lat[1], 30.)
        self.assertEqual(R.lat[2], 30.)
        self.assertEqual(R.lat[3], -50.)

    def test_Region_Index(self):
        R = RegionIndex(77, 20., 30., 50., 60., label='TestIndex')
        self.assertTrue(R.lon is None)
        self.assertTrue(R.lat is None)

    def test_Region_InvalidBbox(self):
        R = RegionGeneric(5, label='test')
        with self.assertRaises(ValueError):
            R._check_bbox_validity(5, 4, 1, 2)
        with self.assertRaises(ValueError):
            R._check_bbox_validity(3, 4, 2, 1)

    def test_Region_label(self):
        R = RegionGeneric(5, label='test xx')
        self.assertEqual(R._get_label(), 'testxx')

    def test_Region_Bbox(self):
        R = RegionIndex(555, 2, 5, 3, 7, label='test')
        x1,x2,y1,y2 = R.get_Bbox()
        self.assertEqual(x1, (2,3))
        self.assertEqual(x2, (2,7))
        self.assertEqual(y1, (5,7))
        self.assertEqual(y2, (5,3))

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


