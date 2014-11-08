# -*- coding: utf-8 -*-

"""
This file is part of pyCMBS. (c) 2012-2014
For COPYING and LICENSE details, please refer to the file
COPYRIGHT.md
"""

import os

"""
This module contains generic utility functions
"""

def get_month_string(m, numeric=False):
    """
    get month as string
    """
    if numeric:
        return str(m).zfill(2)

    d = {}
    d.update({1 : 'JAN'})
    d.update({2 : 'FEB'})
    d.update({3 : 'MAR'})
    d.update({4 : 'APR'})
    d.update({5 : 'MAY'})
    d.update({6 : 'JUN'})
    d.update({7 : 'JUL'})
    d.update({8 : 'AUG'})
    d.update({9 : 'SEP'})
    d.update({10 : 'OCT'})
    d.update({11 : 'NOV'})
    d.update({12 : 'DEC'})

    return d[m]


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
        self.eol = '\n'

    def convert(self, filename=None, mode='w'):
        """
        convert dictionary and store in ASCII file

        Parameters
        ----------
        filename : str
            name of outputfile. If given, then results will be stored
        mode : str
            ['w','a']: access mode for output file
            'w' : write = overwrites already existing file
            'a' : append = appends values to file. Makes only sense if
            the keys in the dictionary are all the same and in same order!
        """
        header, value = self._convert(self.x)
        if filename is not None:
            if mode == 'w':
                if os.path.exists(filename):
                    os.remove(filename)
            f = open(filename, mode=mode)
            if mode == 'w':
                f.write(header + self.eol)
            f.write(value + self.eol)
        return header, value

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

        if len(keys) == 0:
            return '', ''

        for k in keys:
            if parent == '':
                pstr = str(k)
            else:
                pstr = parent + self.tagsep + str(k)
            if isinstance(d[k], dict):
                h1, s1 = self._convert(d[k], h='', s='', parent=pstr)
                h += h1
                s += s1
            else:
                newtag = parent + self.tagsep + str(k)
                if newtag[0] == self.tagsep:
                    newtag = newtag[1:]
                h += newtag
                s += str(d[k])

            s += self.fieldsep
            h += self.fieldsep

        def _remove_trailing_sep(S, sep):
            # remove multiple separators at the end
            if S[-1:] == sep:
                S = S[:-1]
                S = _remove_trailing_sep(S, sep)
            return S

        ho = _remove_trailing_sep(h, self.fieldsep)
        so = _remove_trailing_sep(s, self.fieldsep)

        return ho, so
