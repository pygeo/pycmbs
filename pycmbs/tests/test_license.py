#!/usr/bin/env python
"""
This file is part of pyCMBS.
For COPYRIGHT, LICENSE and AUTHORSHIP please referr to
the pyCMBS licensing details.
"""

import re
import glob
import unittest

class TestCodingStandards(unittest.TestCase):
    """
    test coding standards: check for license,
                                     coda,
                                     encoding
    """
    def test_PythonFiles_HaveLicenseText(self):
        pyfiles = find_python_files()
        files_missing_license = []
        for filename in pyfiles:
            if license_missing(filename) is True:
                files_missing_license.append(filename)
        
        self.assertEquals(len(files_missing_license), 0,
                str(files_missing_license))

def find_python_files():
    """
    find all python files
    """
    python_files = glob.glob('./*.py')
    
    return python_files

def license_missing(filename):
    license_string = \
    "This file is part of pyCMBS. (c) 2012-2014" + "\n" + \
    "For COPYING and LICENSE details, please refer to the files"  + "\n" + \
    "LICENSE.md and COPYRIGHT.md"
    license_missing = True
    fh = open(filename, 'r')
    file_contents = fh.read()
    if license_string in file_contents:
        license_missing = False
    fh.close()
    return license_missing    
