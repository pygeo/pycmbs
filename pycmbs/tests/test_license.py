# -*- coding: utf-8 -*-
"""
This file is part of pyCMBS. (c) 2012-2014"
For COPYING and LICENSE details, please refer to the file"
COPYRIGHT.md
"""

import os
import glob
import unittest

class TestCodingStandards(unittest.TestCase):
    """
    test coding standards: check for license
    """
    def test_PythonFiles_HaveLicenseText(self):
        pyfiles = find_python_files()
        files_missing_license = []
        skip_files = ['netcdftime.py']
        for filename in pyfiles:
            file_basename = os.path.basename(filename)
            if license_missing(filename) is True and file_basename not in skip_files:
                files_missing_license.append(filename)
        
        self.assertEquals(len(files_missing_license), 0,
                str(files_missing_license))

def find_python_files():
    """
    find all python files
    """
    python_files = glob.glob('./scripts/*.py')
    python_files += glob.glob('./scripts/*/*.py')
    python_files += glob.glob('./pycmbs/*.py')
    python_files += glob.glob('./pycmbs/tests.py')
    python_files += glob.glob('./pycmbs/benchmarking/*.py')
    python_files += glob.glob('./pycmbs/benchmarking/tests/*.py')
    return python_files

def license_missing(filename):
    license_string = \
    "This file is part of pyCMBS. (c) 2012-2014" + "\n" + \
    "For COPYING and LICENSE details, please refer to the file"  + "\n" + \
    "COPYRIGHT.md"
    license_missing = True
    fh = open(filename, 'r')
    file_contents = fh.read()
    if license_string in file_contents:
        license_missing = False
    fh.close()
    return license_missing    
