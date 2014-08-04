# -*- coding: utf-8 -*-
"""
This file is part of pyCMBS. (c) 2012-2014
For COPYING and LICENSE details, please refer to the file
COPYRIGHT.md
"""

"""
run all pyCMBS examples in a sequence
"""

import os

examples=['example-01.py','example-02.py','example-03.py','icon_example.py']
for e in examples:
    os.system('python ' + e)

