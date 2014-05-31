# -*- coding: UTF-8 -*-
from unittest import TestCase

from nose.tools import assert_raises
from pycmbs.utils import Dict2TXT


class TestStatistic(TestCase):

    def setUp(self):
        d = dict([('a', 5), ('b', 10), ('c', {'x':1,'y': {'AA' : 77, 'BB' : 'test'}})])
        #~ d.update({'c' : {'x' : 3, 'y': 7}})
        self.x = d


    def test_init_false(self):
        x = 'nothing'
        with self.assertRaises(ValueError):
            D = Dict2TXT(x)

    def test_conversion(self):
        print self.x
        D = Dict2TXT(self.x)
        h, s = D.convert()
        print ''
        print 'FINAL RESULT'
        print h
        print s
        stop

