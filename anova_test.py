#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Alexander Loew"
__version__ = "0.1"
__date__ = "2012/10/29"
__email__ = "alexander.loew@zmaw.de"

'''
# Copyright (C) 2012 Alexander Loew, alexander.loew@zmaw.de
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





'''
test anova with Data object

validation o.k. for 1-way and 2-way ANOVA
'''

import numpy as np

from pyCMBS import *

from pylab import *

from anova import *






#~ a.add_experiment('e2')


#~ d2=Data(None,None)
#~ d3=Data(None,None)
#~ d4=Data(None,None)


#~ x=rand(10,100,200)
#~ d1.data = np.ma.array(x,mask=x>0.9)
#~ a.add_data('e2',d1)
#~ del x
#~
#~ x=rand(10,100,200)
#~ d2.data = np.ma.array(x,mask=x>0.9)
#~ a.add_data('e2',d2)
#~ del x
#~
#~ x=rand(10,100,200)
#~ d3.data = np.ma.array(x,mask=x>2.)
#~ a.add_data('e1',d3)
#~ del x
#~
#~ x=rand(10,100,200)
#~ d4.data = np.ma.array(x,mask=x>2.)
#~ a.add_data('e1',d4)

if False:
    # TWO WAY ANOVA
    #http://people.richland.edu/james/lecture/m170/ch13-2wy.html

    a = ANOVA()

    a.add_experiment('e1')
    a.add_experiment('e2')
    a.add_experiment('e3')

    #ens #1
    d1=Data(None,None)
    x = np.zeros((5,1,1))
    x[:,0,0] = [106,95,94,103,100]
    d1.data = np.ma.array(x,mask=np.isnan(x))
    a.add_data('e1',d1)

    d2=Data(None,None)
    x = np.zeros((5,1,1))
    x[:,0,0] = [110,98,100,108,105]
    d2.data = np.ma.array(x,mask=np.isnan(x))
    a.add_data('e2',d2)

    d3=Data(None,None)
    x = np.zeros((5,1,1))
    x[:,0,0] = [94,86,98,99,94]
    d3.data = np.ma.array(x,mask=np.isnan(x))
    a.add_data('e3',d3)

    #ens #2
    d4=Data(None,None)
    x = np.zeros((5,1,1))
    x[:,0,0] = [110,100,107,104,102]
    d4.data = np.ma.array(x,mask=np.isnan(x))
    a.add_data('e1',d4)

    d5=Data(None,None)
    x = np.zeros((5,1,1))
    x[:,0,0] = [112,99,101,112,107]
    d5.data = np.ma.array(x,mask=np.isnan(x))
    a.add_data('e2',d5)

    d5=Data(None,None)
    x = np.zeros((5,1,1))
    x[:,0,0] = [97,87,99,101,98]
    d5.data = np.ma.array(x,mask=np.isnan(x))
    a.add_data('e3',d5)

    a.analysis()


if True:
    #one way ANOVA example
    #http://adorio-research.org/wordpress/?p=1102

    B = ANOVA()
    B.add_experiment('A')

    d=Data(None,None)
    x = np.zeros((3,1,1))
    #~ x[:,0,0] = [48, 49, 50, 49]
    x[:,0,0] = [48, 47, 49]

    d.data = np.ma.array(x,mask=np.isnan(x))
    B.add_data('A',d)

    d=Data(None,None)
    x = np.zeros((3,1,1))
    #~ x[:,0,0] = [47, 49, 48, 48]
    x[:,0,0] = [49, 49, 51]
    d.data = np.ma.array(x,mask=np.isnan(x))
    B.add_data('A',d)

    d=Data(None,None)
    x = np.zeros((3,1,1))
    #~ x[:,0,0] = [49, 51, 50, 50]
    x[:,0,0] = [50, 48, 50]
    d.data = np.ma.array(x,mask=np.isnan(x))
    B.add_data('A',d)

    d=Data(None,None)
    x = np.zeros((3,1,1))
    #~ x[:,0,0] = [49, 51, 50, 50]
    x[:,0,0] = [49, 48, 50]
    d.data = np.ma.array(x,mask=np.isnan(x))
    B.add_data('A',d)

    B.analysis(analysis_type='one')




    #~ groups = [[48, 49, 50, 49],
              #~ [47, 49, 48, 48],
              #~ [49, 51, 50, 50]]







