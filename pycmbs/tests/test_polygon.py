# -*- coding: utf-8 -*-
"""
This file is part of pyCMBS. (c) 2012-2014
For COPYING and LICENSE details, please refer to the file
COPYRIGHT.md
"""

from unittest import TestCase
import unittest

from pycmbs.data import Data
import scipy as sc
import matplotlib.pylab as pl
import numpy as np
from scipy import stats
import datetime
import matplotlib.pylab as plt

from nose.tools import assert_raises

from pycmbs.polygon import Polygon, Raster


class TestData(unittest.TestCase):

    def setUp(self):
        self.D = Data(None, None)
        self.D._init_sample_object(nt=1000, ny=1, nx=1)

    def test_is_closed(self):
        poly = [(150.,20.), (-160.,30.), (-170.,10.), (170.,10.)]
        P = Polygon(3, poly)
        self.assertFalse(P.is_closed())

        poly1 = [(150.,20.), (-160.,30.), (-170.,10.), (170.,10.), (150.,20.)]
        P1 = Polygon(3, poly1)
        self.assertTrue(P1.is_closed())

    def test_convertOGR(self):
        poly = [(150.,20.), (-160.,30.), (-170.,10.), (170.,10.)]
        P = Polygon(3, poly)
        A = P.convertToOGRPolygon()
        B = P.convertToOGRPolygon(ensure_positive=True)

    def test_shift(self):
        poly3 = [(150.,20.), (-160.,30.), (-170.,10.), (170.,10.)]
        P3 = Polygon(3, poly3)
        P3._shift_coordinates()  # shift longitudes by 200 degree
        self.assertEqual(P3.poly[0][0], 150.)
        self.assertEqual(P3.poly[1][0], 200.)
        self.assertEqual(P3.poly[2][0], 190.)
        self.assertEqual(P3.poly[3][0], 170.)

    def test_point_in_polygon(self):
        x = 1
        y = 1
        poly = [(0,0), (2,0), (2,2), (0,2)]
        P = Polygon(1, poly)
        self.assertTrue(P.point_in_poly(x,y))
        x = 4
        y = 4
        self.assertFalse(P.point_in_poly(x,y))

    def test_point_in_polygon_latlon(self):
        # test for point in polygon across dateline
        x1 = -175.
        y1 = 50.
        poly1= [(150.,60.), (-160.,60.), (-170.,45.), (170.,45.)]
        P1 = Polygon(1, poly1)
        self.assertTrue(P1.point_in_poly_latlon(x1,y1))


    def test_polygon_min_max(self):
        x = 1
        y = 1
        poly = [(-5,0), (2,0), (2,2), (0,6)]
        P = Polygon(1, poly)

        bbox = P.bbox()
        self.assertEqual(bbox[0], -5.)
        self.assertEqual(bbox[1], 2.)
        self.assertEqual(bbox[2], 0.)
        self.assertEqual(bbox[3], 6.)

    def test_raster_wrong_geometry(self):
        lon = np.random.random((10,20))
        lat = np.random.random((11,20))
        with self.assertRaises(ValueError):
            R = Raster(lon, lat)

    def test_raster_wrong_latlon(self):
        lon = np.random.random(10)
        lat = np.random.random(10)
        print lon.ndim
        with self.assertRaises(ValueError):
            R = Raster(lon, lat)

    def test_raster_no_Polygon(self):
        lon = np.random.random((10,20))
        lat = np.random.random((10,20))
        R = Raster(lon, lat)
        with self.assertRaises(ValueError):
            P = np.arange(10)
            R._rasterize_single_polygon(P)

    def test_raster_single_polygon(self):
        lon = np.linspace(-180., 180., 361)
        lat = np.linspace(-90., 90., 181)
        LON,LAT=np.meshgrid(lon, lat)

        # test a single polygon
        poly = [(-10.,-10.), (-10.,20), (15.,0.), (0.,-25.)]
        P = Polygon(5, poly)
        R=Raster(LON,LAT)
        R.mask = np.zeros(LON.shape)*np.nan
        R._rasterize_single_polygon(P)
        R.mask = np.ma.array(R.mask, mask=np.isnan(R.mask))

        u = np.unique(R.mask[~R.mask.mask])
        self.assertTrue(len(u)==1)
        self.assertTrue(5. in u)

    #~ def test_raster_single_polygon_fast(self):
        #~ lon = np.linspace(-180., 180., 361)
        #~ lat = np.linspace(-90., 90., 181)
        #~ LON,LAT=np.meshgrid(lon, lat)
#~
        #~ # test a single polygon
        #~ poly = [(-10.,-10.), (-10.,20), (15.,0.), (0.,-25.)]
        #~ P = Polygon(5, poly)
        #~ R=Raster(LON,LAT)
        #~ R.mask = np.zeros(LON.shape)*np.nan
        #~ R._rasterize_single_polygon(P, method='fast')
        #~ R.mask = np.ma.array(R.mask, mask=np.isnan(R.mask))
#~
        #~ u = np.unique(R.mask[~R.mask.mask])
        #~ self.assertTrue(len(u) == 1)
        #~ self.assertTrue(5. in u)


    def test_raster_multiple_polygon(self):  # this is quite slow!
        lon = np.linspace(-180., 180., 361)
        lat = np.linspace(-90., 90., 181)
        LON,LAT=np.meshgrid(lon, lat)

        # test a single polygon
        poly=[]
        poly1 = [(-10.,-10.), (-10.,20), (15.,0.), (0.,-15.)]
        poly.append(Polygon(1, poly1))

        poly2 = [(-50.,-80.), (-50.,-70.), (-40.,-70.), (-40.,-75.)]
        poly.append(Polygon(2, poly2))

        R=Raster(LON,LAT)
        R.rasterize_polygons(poly)

        u = np.unique(R.mask[~R.mask.mask])
        self.assertTrue(len(u)==2)
        self.assertTrue(1 in u)
        self.assertTrue(2 in u)


    #~ def test_raster_multiple_polygon_fast(self):
        #~ lon = np.linspace(-180., 180., 361)
        #~ lat = np.linspace(-90., 90., 181)
        #~ LON,LAT=np.meshgrid(lon, lat)
#~
        #~ # test a single polygon
        #~ poly=[]
        #~ poly1 = [(-10.,-10.), (-10.,20), (15.,0.), (0.,-15.)]
        #~ poly.append(Polygon(1, poly1))
#~
        #~ poly2 = [(-50.,-80.), (-50.,-70.), (-40.,-70.), (-40.,-75.)]
        #~ poly.append(Polygon(2, poly2))
#~
        #~ R=Raster(LON,LAT)
        #~ R.rasterize_polygons(poly, method='fast')
#~
        #~ u = np.unique(R.mask[~R.mask.mask])
        #~ self.assertTrue(len(u)==2)
        #~ self.assertTrue(1 in u)
        #~ self.assertTrue(2 in u)


