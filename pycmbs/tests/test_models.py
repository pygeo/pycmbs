# -*- coding: utf-8 -*-
"""
This file is part of pyCMBS.
(c) 2012- Alexander Loew
For COPYING and LICENSE details, please refer to the LICENSE file
"""

from pycmbs.benchmarking.models import MeanModel, Model
import numpy as np
from pycmbs.data import Data
import numpy.testing as npt

def test_mean_model():

    #The following code provides a routine that allows to validate the MeanModel() class
    print ('Jetzt gehts los')
    # generate some sample data ---
    x = Data(None, None)
    x.data = np.random.random((10,20,30))
    x.label='nothing'

    y = x.mulc(0.3)
    z = x.mulc(0.5)
    m = x.add(y).add(z).divc(3.)
    r = m.div(x)  # gives 0.6 as reference solution

    # generate Model instances and store Data objects as 'variables' ---
    dic_variables = ['var1', 'var2']
    X = Model(None, dic_variables, name='x', intervals='season')
    X.variables = {'var1': x, 'var2': x}
    Y = Model(None, dic_variables, name='y', intervals='season')
    Y.variables = {'var1': y, 'var2': y}
    Z = Model(None, dic_variables, name='z', intervals='season')
    Z.variables={'var1': z, 'var2': z}

    #... now try multimodel ensemble
    M=MeanModel(dic_variables,intervals='season')
    M.add_member(X)
    M.add_member(Y)
    M.add_member(Z)
    M.ensmean()  # calculate ensemble mean
    # print M.variables['var2'].div(x).data #should give 0.6
    npt.assert_equal(np.all(np.abs(1. - M.variables['var2'].div(x).data/0.6) < 0.00000001), True)

def xxxxtest_median_model():
    x = Data(None, None)
    x.label = 'nothing'
    d = np.random.random((100, 1, 1))
    x.data = np.ma.array(d, mask= d!=d)

    # odd number and no masked values
    a = x.copy()
    a.data[:, 0, 0] = 1.
    b = x.copy()
    b.data[:, 0, 0] = 3.
    c = x.copy()
    c.data[:, 0, 0] = 2.
    d = x.copy()
    d.data[:, 0, 0] = 5.
    e = x.copy()
    e.data[:, 0, 0] = 4.

    m = MedianModel()
    m.add_member(a)
    m.add_member(b)
    m.add_member(c)
    m.add_member(d)
    m.add_member(e)
    m.ensmedian()

    # should give the value of 3. for all timesteps

    del m

    # even number and no masked values
    a = x.copy()
    a.data[:, 0, 0] = 1.
    b = x.copy()
    b.data[:, 0, 0] = 3.
    c = x.copy()
    c.data[:, 0, 0] = 2.
    d = x.copy()
    c.data[:, 0, 0] = 4.

    m = MedianModel()
    m.add_member(a)
    m.add_member(b)
    m.add_member(c)
    m.add_member(d)

    m.ensmedian()

    # should give the value of 2.5 for all timesteps
    del m






