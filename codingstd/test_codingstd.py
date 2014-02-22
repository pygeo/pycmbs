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
    license_text = 'This file is part of pyCMBS.\nFor COPYRIGHT, LICENSE and AUTHORSHIP please referr to\nthe pyCMBS licensing details.'
    license_missing = True
    fh = open(filename, 'r')
    file_contents = fh.read()
    if license_text in file_contents:
        license_missing = False
    fh.close()
    return license_missing
    
