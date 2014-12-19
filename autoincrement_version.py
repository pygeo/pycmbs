# -*- coding: UTF-8 -*-

"""
This file is part of pyCMBS. (c) 2012-2014
For COPYING and LICENSE details, please refer to the file
COPYRIGHT.md
"""

"""
function to autoincrement version
"""
import json
import os
import sys

# read current version from file
cpath = os.path.dirname(os.path.realpath(__file__))

vfile = cpath + os.sep + 'pycmbs' + os.sep + 'version.json'

v = json.load(open(vfile, 'r'))
V = v.split('.')

# retrieve last subversion number
try:
    sub = int(V[-1])
except:
    print('Version number not appropriate for autoincrement!')
    sys.exit(1)

sub += 1
vn = ''
for i in xrange(len(V) - 1):
    vn += V[i] + '.'
vn += str(sub)

if os.path.exists(vfile):
    os.remove(vfile)
json.dump(vn, open(vfile, 'w'))
print('Incremented version to: ' + vn)
