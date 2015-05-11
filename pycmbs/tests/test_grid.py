# -*- coding: utf-8 -*-
"""
This file is part of pyCMBS.
(c) 2012- Alexander Loew
For COPYING and LICENSE details, please refer to the LICENSE file
"""

import unittest
from pycmbs import grid

from nose.tools import assert_raises
import numpy as np


class TestPycmbsGrid(unittest.TestCase):

    def setUp(self):
        self.lon = np.linspace(-120., 130.)
        self.lat = np.linspace(-30., 30.)
        self.grid  = grid.Grid(np.deg2rad(self.lat), np.deg2rad(self.lon), sphere_radius=6000000.)

    def test_DummyTest(self):
        pass

    def test_grid_init(self):
        with self.assertRaises(ValueError):
            G = grid.Grid(np.deg2rad(self.lat), np.deg2rad(self.lon))

        G = grid.Grid(np.deg2rad(self.lat), np.deg2rad(self.lon), sphere_radius=6000000.)

    def test_grid_cellarea(self):
        with self.assertRaises(ValueError):
            self.grid.calc_cell_area()

    def test_grid_plot(self):
        self.grid.plot()

    def test_grid_plot_voronoi(self):
        self.grid.plot_voronoi()

    def test_grid_plot_delaunay(self):
        self.grid.plot_delaunay_grid()

    def test_grid_draw_edge(self):
        self.grid.draw_edges()






if __name__ == "__main__":
    unittest.main()
