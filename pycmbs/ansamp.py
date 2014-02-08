#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Alexander Loew"
__version__ = "0.1.4"
__date__ = "2012/10/29"
__email__ = "alexander.loew@mpimet.mpg.de"

'''
# Copyright (C) 2012 Alexander Loew, alexander.loew@mpimet.mpg.de
# See COPYING file for copying and redistribution conditions.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; version 2 of the License.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
'''


"""
file    twoway_interaction.py
author  Ernesto P. Adorio, Ph.D.
        ernesto.adorio@gmail.com
        UPDEPP at Clarkfield
desc    Performs an anova with interaction
        on replicated input data.
        Each block must have the same number
        of input values.
version 0.0.2 Sep 12, 2011

"""

import numpy as np


def twoway_interaction(groups, format="html"):
    b = len(groups[0][0])
    a = len(groups)
    c = len(groups[0])
    groupsums = [0.0] * c

    print "blocks, a, c=", b, a, c

    print "Input groups:"
    v = 0.0  # total variation
    vs = 0.0  # subtotal variation
    vr = 0.0  # variation between rows
    GT = 0
    for i in range(a):
       vsx = 0.0
       vrx = 0.0
       for j in range(c):
         vsx = sum(groups[i][j])
         groupsums[j] += vsx
         print "debug vsx", vsx
         vrx += vsx
         vs += vsx * vsx
         for k in range(b):
            x = groups[i][j][k]
            v += x * x
            GT += x
       vr += vrx * vrx

    print "groupsums=", groupsums, vs
    print 'Total variation: ', v

    totadjustment = GT * GT / (a * b * c)
    vs = vs / b - totadjustment
    vr = vr / (b * c) - totadjustment
    v -= totadjustment
    vc = sum([x * x for x in groupsums]) / (a * b) - totadjustment
    vi = vs - vr - vc
    ve = v - (vr + vc + vi)
    print "debug vs, vr, vc=", vs, vr, vc, v

    dfvr = (a - 1)
    dfvc = (c - 1.0)
    dfvi = ((a - 1) * (c - 1))
    dfve = (a * c * (b - 1))
    dfvs = a * c - 1
    dfv = (a * b * c - 1)
    mvr = vr / (dfvr)
    mvc = vc / (dfvc)
    mvi = vi / dfvi
    mve = ve / dfve
    Fr = mvr / mve
    Fc = mvc / mve
    Fi = mvi / mve

    from scipy import stats

    pvalr = 1.0 - stats.f.cdf(Fr, dfvr, dfve)
    pvalc = 1.0 - stats.f.cdf(Fc, dfvc, dfve)
    pvali = 1.0 - stats.f.cdf(Fi, dfvi, dfve)

    if format == "html":
       output = """
    <table border="1">
    <tr><th>Variation  </th><th>Sum of Squares</th><th>  df</th><th>  Mean Sum of Squares</th><th>   F-value</th><th> p-value</th></tr>
    r><td>Rows(treatments) </td><td>%f</td><td>    %d</td><td>     %f</td> <td> %f</td> <td>%f</td></tr>
    <tr><td>Columns(blocks)</td><td>%f</td><td>  %d</td><td>     %f</td> <td> %f</td> <td>%f</td></tr>
    <tr><td>Interaction</td><td>%f</td><td>  %d</td><td>     %f</td> <td> %f</td> <td>%f</td></tr>
    <tr><td>Subtotals </td><td> %f</td><td> %d</td></tr>
    <tr><td>Residuals(random)  </td><td>%f</td><td>  %d</td><td>%f</td></tr>
    <tr><td>Totals</td><td>%f.2 </td><td>%d </td></tr>
    </table>
    """ % (vr, dfvr, mvr, mvr / mve, pvalr,
           vc, dfvc, mvc, mvc / mve, pvalc,
           vi, dfvi, mvi, mvi / mve, pvali,
           vs, dfvs,
           ve, dfve, mve,
           v, dfv)
    else:
        output = [[vr, dfvr, mvr, mvr / mve, pvalr],
         [vc, dfvc, mvc, mvc / mve, pvalc],
         [vi, dfvi, mvi, mvi / mve, pvali],
         [vs, dfvs],
         [ve, dfve, mve],
         [v, dfv]]

    return output

groups = [
           [[64, 72, 74],
            [66, 81, 51],
            [70, 64, 65]
            ],
           [[65, 57, 47],
            [63, 43, 58],
            [58, 52, 67]
            ],
           [[59, 66, 58],
            [68, 71, 39],
            [65, 59, 42]
            ],
           [[58, 57, 53],
            [41, 61, 59],
            [46, 53, 38]
            ]
         ]


# Spiegel, "Probability and Statistics",
#   page 324-326
groups = [
          [[6, 4, 5, 5, 4],
           [5, 7, 4, 6, 8]
           ],
          [[10, 8, 7, 7, 9],
           [7, 9, 12, 8, 8]
           ],
          [[7, 5, 6, 5, 9],
           [9, 7, 5, 4, 6]
           ],
          [[8, 4, 6, 5, 5],
           [5, 7, 9, 7, 10]
           ]
         ]

groups = [
          [[4, 7, 10],
          [6, 13, 12]],

          [[5, 9, 12],
          [6, 15, 13]],

          [[6, 8, 11],
          [4, 12, 10]],

          [[5, 12, 9],
          [4, 12, 13]]

          ]

groups = [
          [[106, 95, 94, 103, 100],
          [110, 98, 100, 108, 105],
          [94, 86, 98, 99, 94]],

          [[110, 100, 107, 104, 102],
          [112, 99, 101, 112, 107],
          [97, 87, 99, 101, 98]]

          ]

g1 = []
for i in range(len(groups)):
    g1.append(list(np.asarray(groups[i]).T))

g1 = groups

output = twoway_interaction(g1, "html")
print output
