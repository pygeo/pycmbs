# -*- coding: utf-8 -*-

"""
This file is part of pyCMBS. (c) 2012-2014
For COPYING and LICENSE details, please refer to the file
COPYRIGHT.md
"""


"""
This module contains generic utility functions
"""


class Dict2TXT(object):
    """
    class to convert a python dictionary to a TXT file
    """

    def __init__(self, x, fieldsep='\t', tagsep=':'):
        """
        Parameters
        ----------
        x : dict
            dictionary to be converted
        fieldsep : str
            separator between fields
        tagsep : str
            separator for column labels
        """
        if not isinstance(x, dict):
            raise ValueError('The input needs to eb a valid python dictionary!')
        self.x = x
        self.fieldsep = fieldsep
        self.tagsep = tagsep

    def convert(self):
        return self._convert(self.x)

    def _convert(self, d, h='', s='', parent=''):
        """
        convert dictionary to a string in recursive manner

        Parameters
        ----------
        d : some input
        """

        if not isinstance(d, dict):
            raise ValueError('Need to provide dictionary!')

        keys = d.keys()
        keys.sort()

        print ''
        print keys
        print d
        for k in keys:
            print '--'
            print 'k = ', k
            if isinstance(d[k], dict):
                h1, s1 = self._convert(d[k], h='', s='', parent=k)
                h += h1
                s += s1
            else:
                newtag = parent + self.tagsep + k
                if newtag[0] == self.tagsep:
                    newtag = newtag[1:]
                h += newtag
                sep = self.fieldsep
                s += str(d[k])
            s += self.fieldsep
            h += sep

        return h, s



