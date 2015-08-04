# -*- coding: utf-8 -*-
"""
This file is part of pyCMBS.
(c) 2012- Alexander Loew
For COPYING and LICENSE details, please refer to the LICENSE file
"""

import numpy as np


def get_albedo_colortable():
    """
    colors(*,i)=[0,     0, 050] & boundary[i]=0.000 & i=i+1
    colors(*,i)=[0,     0, 200] & boundary[i]=0.020 & i=i+1  ; 0.020
    colors(*,i)=[0,     0, 255] & boundary[i]=0.040 & i=i+1  ; 0.040
    colors(*,i)=[255,  24,   0] & boundary[i]=0.060 & i=i+1  ; 0.060
    colors(*,i)=[220,  40,   4] & boundary[i]=0.080 & i=i+1  ; 0.080
    colors(*,i)=[192,  65,   7] & boundary[i]=0.100 & i=i+1  ; 0.100
    colors(*,i)=[129,  25,  14] & boundary[i]=0.120 & i=i+1  ; 0.120
    colors(*,i)=[ 74, 134,   0] & boundary[i]=0.140 & i=i+1  ; 0.140
    colors(*,i)=[152, 186,   0] & boundary[i]=0.160 & i=i+1  ; 0.160
    colors(*,i)=[153, 147,   0] & boundary[i]=0.180 & i=i+1  ; 0.180
    colors(*,i)=[139, 123,   0] & boundary[i]=0.200 & i=i+1  ; 0.200
    colors(*,i)=[125,  99,   0] & boundary[i]=0.220 & i=i+1  ; 0.220
    colors(*,i)=[111,  75,   0] & boundary[i]=0.240 & i=i+1  ; 0.240
    colors(*,i)=[126,  91,  14] & boundary[i]=0.260 & i=i+1  ; 0.260
    colors(*,i)=[141, 108,  28] & boundary[i]=0.280 & i=i+1  ; 0.280
    colors(*,i)=[156, 125,  42] & boundary[i]=0.300 & i=i+1  ; 0.300
    colors(*,i)=[171, 142,  56] & boundary[i]=0.325 & i=i+1  ; 0.325
    colors(*,i)=[186, 159,  71] & boundary[i]=0.350 & i=i+1  ; 0.350
    colors(*,i)=[201, 176,  85] & boundary[i]=0.375 & i=i+1  ; 0.375
    colors(*,i)=[216, 193,  99] & boundary[i]=0.400 & i=i+1  ; 0.400
    colors(*,i)=[231, 210, 113] & boundary[i]=0.450 & i=i+1  ; 0.450
    colors(*,i)=[240, 220, 120] & boundary[i]=0.500 & i=i+1  ; 0.500
    colors(*,i)=[246, 225, 135] & boundary[i]=0.550 & i=i+1  ; 0.550
    colors(*,i)=[246, 235, 155] & boundary[i]=0.600 & i=i+1  ; 0.600
    colors(*,i)=[240, 240, 180] & boundary[i]=0.650 & i=i+1  ; 0.650
    colors(*,i)=[250, 250, 210] & boundary[i]=0.700 & i=i+1  ; 0.750
    colors(*,i)=[230, 253, 200] & boundary[i]=0.750 & i=i+1  ; 0.700

    which means for instance that the interval 0.18 - 0.20 is coded with the RGB value [139,123,0]. Missing values (255) are coded in white. If you multiply these intervals by 254 you have the equivalent intervals directly in the way the albedo product is coded.
    """
    ct = [[0, 0, 050],
          [0, 0, 200],
          [0, 0, 255],
          [255, 24, 0],
          [220, 40, 4],
          [30, 70, 0],
          [50, 100, 0],
          [74, 134, 0],
          [152, 186, 0],
          [153, 147, 0],
          [139, 123, 0],
          [125, 99, 0],
          [111, 75, 0],
          [126, 91, 14],
          [141, 108, 28],
          [156, 125, 42],
          [171, 142, 56],
          [186, 159, 71],
          [201, 176, 85],
          [216, 193, 99],
          [231, 210, 113],
          [240, 220, 120],
          [246, 225, 135],
          [246, 235, 155],
          [240, 240, 180],
          [250, 250, 210],
          [230, 253, 200]]

    ct = np.asarray(ct)
    ct = ct / 255.

    # define  boundaries
    lbounds = [0.000,
               0.020,
               0.040,
               0.060,
               0.080,
               0.100,
               0.120,
               0.140,
               0.160,
               0.180,
               0.200,
               0.220,
               0.240,
               0.260,
               0.280,
               0.300,
               0.325,
               0.350,
               0.375,
               0.400,
               0.450,
               0.500,
               0.550,
               0.600,
               0.650,
               0.700,
               0.750]
    return lbounds, ct


def get_albedo_colortable1():
    """
    colors(*,i)=[0,     0, 050] & boundary[i]=0.000 & i=i+1
    colors(*,i)=[0,     0, 200] & boundary[i]=0.020 & i=i+1  ; 0.020
    colors(*,i)=[0,     0, 255] & boundary[i]=0.040 & i=i+1  ; 0.040
    colors(*,i)=[255,  24,   0] & boundary[i]=0.060 & i=i+1  ; 0.060
    colors(*,i)=[220,  40,   4] & boundary[i]=0.080 & i=i+1  ; 0.080
    colors(*,i)=[192,  65,   7] & boundary[i]=0.100 & i=i+1  ; 0.100
    colors(*,i)=[129,  25,  14] & boundary[i]=0.120 & i=i+1  ; 0.120
    colors(*,i)=[ 74, 134,   0] & boundary[i]=0.140 & i=i+1  ; 0.140
    colors(*,i)=[152, 186,   0] & boundary[i]=0.160 & i=i+1  ; 0.160
    colors(*,i)=[153, 147,   0] & boundary[i]=0.180 & i=i+1  ; 0.180
    colors(*,i)=[139, 123,   0] & boundary[i]=0.200 & i=i+1  ; 0.200
    colors(*,i)=[125,  99,   0] & boundary[i]=0.220 & i=i+1  ; 0.220
    colors(*,i)=[111,  75,   0] & boundary[i]=0.240 & i=i+1  ; 0.240
    colors(*,i)=[126,  91,  14] & boundary[i]=0.260 & i=i+1  ; 0.260
    colors(*,i)=[141, 108,  28] & boundary[i]=0.280 & i=i+1  ; 0.280
    colors(*,i)=[156, 125,  42] & boundary[i]=0.300 & i=i+1  ; 0.300
    colors(*,i)=[171, 142,  56] & boundary[i]=0.325 & i=i+1  ; 0.325
    colors(*,i)=[186, 159,  71] & boundary[i]=0.350 & i=i+1  ; 0.350
    colors(*,i)=[201, 176,  85] & boundary[i]=0.375 & i=i+1  ; 0.375
    colors(*,i)=[216, 193,  99] & boundary[i]=0.400 & i=i+1  ; 0.400
    colors(*,i)=[231, 210, 113] & boundary[i]=0.450 & i=i+1  ; 0.450
    colors(*,i)=[240, 220, 120] & boundary[i]=0.500 & i=i+1  ; 0.500
    colors(*,i)=[246, 225, 135] & boundary[i]=0.550 & i=i+1  ; 0.550
    colors(*,i)=[246, 235, 155] & boundary[i]=0.600 & i=i+1  ; 0.600
    colors(*,i)=[240, 240, 180] & boundary[i]=0.650 & i=i+1  ; 0.650
    colors(*,i)=[250, 250, 210] & boundary[i]=0.700 & i=i+1  ; 0.750
    colors(*,i)=[230, 253, 200] & boundary[i]=0.750 & i=i+1  ; 0.700

    which means for instance that the interval 0.18 - 0.20 is coded with the RGB value [139,123,0]. Missing values (255) are coded in white. If you multiply these intervals by 254 you have the equivalent intervals directly in the way the albedo product is coded.
    """
    ct = [[0, 0, 050],
          [0, 0, 200],
          [0, 0, 255],
          [255, 24, 0],
          [220, 40, 4],
          [192, 65, 7],
          [129, 25, 14],
          [74, 134, 0],
          [152, 186, 0],
          [153, 147, 0],
          [139, 123, 0],
          [125, 99, 0],
          [111, 75, 0],
          [126, 91, 14],
          [141, 108, 28],
          [156, 125, 42],
          [171, 142, 56],
          [186, 159, 71],
          [201, 176, 85],
          [216, 193, 99],
          [231, 210, 113],
          [240, 220, 120],
          [246, 225, 135],
          [246, 235, 155],
          [240, 240, 180],
          [250, 250, 210],
          [230, 253, 200]]

    ct = np.asarray(ct)
    ct = ct / 255.

    # define  boundaries
    lbounds = [0.000,
               0.020,
               0.040,
               0.060,
               0.080,
               0.100,
               0.120,
               0.140,
               0.160,
               0.180,
               0.200,
               0.220,
               0.240,
               0.260,
               0.280,
               0.300,
               0.325,
               0.350,
               0.375,
               0.400,
               0.450,
               0.500,
               0.550,
               0.600,
               0.650,
               0.700,
               0.750]
    return lbounds, ct


class ColorMapGenerator(object):

    """
    Generate colormaps from RGB value lists
    """

    def __init__(self):
        pass

    def albedo(self):
        lb, ct = get_albedo_colortable()
        return self.rgb_to_cmap(lb, ct, name='albedo')

    def rgb_to_cmap(self, lbound, rgb, name='mymap'):
        """
        generate a colormap based on a list of lower boundaries
        and an RGB list

        inspired by http://faculty.washington.edu/rjl/clawpack/trunk/python/pyclaw/plotters/colormaps.py

        Parameters
        ----------
        lbound : array
            array which specifies the lower boundaries
        rgb : array
            a list of [n,3] dimension which specifies the RGB values

        Example
        -------
        > import matplotlib.pylab as plt
        > x = plt.randn(100,100) + 1.
        > lb, ct = get_albedo_colortable()
        > C = ColorMapGenerator()
        > cma = C.rgb_to_cmap(lb, ct)
        > plt.imshow(x, cmap=cma, vmin=0., vmax=1.)
        > plt.colorbar()
        > plt.show()
        """

        from matplotlib.colors import LinearSegmentedColormap, ColorConverter
        import numpy as np

        if len(lbound) != len(rgb):
            raise ValueError(
                'Inconsistent geometries for boundaries and RGB table')

        lbound = np.asarray(lbound)
        # check that boundaries in ascending order
        if np.any(np.diff(lbound) < 0.):
            raise ValueError('Boundaries are not in ascending order!')

        n = len(lbound)
        bmin = lbound.min()
        bmax = lbound.max()
        CC = ColorConverter()
        R = []
        G = []
        B = []
        x = []
        for i in xrange(n):
            R.append(rgb[i, 0])
            G.append(rgb[i, 1])
            B.append(rgb[i, 2])
            x.append((lbound[i] - bmin) / (bmax - bmin))
        x = np.asarray(x)

        cmap_dict = {}
        red_list = []
        green_list = []
        blue_list = []
        # in case of homogeneous colors, generate a tuple with 0 and 1 (see:
        # http://stackoverflow.com/questions/16267143/matplotlib-single-colored-colormap-with-saturation)
        if bmax == bmin:
            for xx in [0., 1.]:
                red_list.append((xx, R[i], R[i]))
                green_list.append((xx, G[i], G[i]))
                blue_list.append((xx, B[i], B[i]))
        else:
            for i in xrange(n):
                red_list.append((x[i], R[i], R[i]))
                green_list.append((x[i], G[i], G[i]))
                blue_list.append((x[i], B[i], B[i]))

        cmap_dict['red'] = red_list
        cmap_dict['green'] = green_list
        cmap_dict['blue'] = blue_list

        cmap = LinearSegmentedColormap(name, cmap_dict)
        return cmap
