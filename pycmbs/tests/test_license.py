# -*- coding: utf-8 -*-
"""
This file is part of pyCMBS.
(c) 2012- Alexander Loew
For COPYING and LICENSE details, please refer to the LICENSE file
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
        for filename in pyfiles:
            file_basename = os.path.basename(filename)
            if license_missing(filename):  # is True #and file_basename not in skip_files:
                files_missing_license.append(filename)

        self.assertEquals(len(files_missing_license), 0,
                str(files_missing_license))

def find_python_files():
    """
    find all python files in pyCMBS installation
    """

    # get path of current file and set root path relative to it
    path = os.path.dirname(os.path.realpath(__file__)) + os.sep + '..' + os.sep + '..' + os.sep

    # *.py
    pyfiles = [os.path.join(dirpath, f)
        for dirpath, dirnames, files in os.walk(path)
        for f in files if f.endswith('.py')]

    # *.pyx
    pyxfiles = [os.path.join(dirpath, f)
        for dirpath, dirnames, files in os.walk(path)
        for f in files if f.endswith('.pyx')]

    res = pyfiles + pyxfiles

    return res



def license_missing(filename):
    license_string = \
    "This file is part of pyCMBS." + "\n" + \
    "(c) 2012- Alexander Loew" + "\n" + \
    "For COPYING and LICENSE details, please refer to the LICENSE file"

    # check first if directory shall be skipped
    skip_dirs = ['docsrc']
    skip_files = ['emd.py']
    skip_tags = ['cartopy']
    for sd in skip_dirs:  # skip predefined directories
        if sd in os.path.dirname(filename):
            return False


    license_missing = True
    fh = open(filename, 'r')
    file_contents = fh.read()
    if license_string in file_contents:
        license_missing = False
    fh.close()
    if os.path.basename(filename) in skip_files:
        license_missing = False
    for k in skip_tags:
        if k in filename:
            license_missing = False
    return license_missing
