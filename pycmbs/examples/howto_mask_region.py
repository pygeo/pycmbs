# -*- coding: utf-8 -*-

"""
This file is part of pyCMBS. (c) 2012-2014
For COPYING and LICENSE details, please refer to the file
COPYRIGHT.md
"""

from pycmbs.examples import download
import matplotlib.pyplot as plt
from pycmbs.mapping import map_plot
from pycmbs.data import Data
import numpy as np

plt.close('all')

# read some data as Data object
filename = download.get_sample_file(name='air', return_object=False)
air = Data(filename, 'air', read=True)

# generate some mask that fits geometry of data. Could be also read from some
# file. We mimic here some irregular mask
mask = np.zeros(air.data[0,:,:].shape)
mask[20:30, 50:60] = 1.
mask[40:60, 100:120] = 1.
mask = np.asarray(mask).astype('bool')  # make a bool array




# generate figure for illustration purposes
f = plt.figure()
ax1 = f.add_subplot(2,2,1)
ax2 = f.add_subplot(2,2,2)
ax3 = f.add_subplot(2,2,3)
ax4 = f.add_subplot(2,2,4)

f1 = map_plot(air, ax=ax1, title='This is the original data (unprojected)')
ax2.imshow(mask)
ax2.set_title('This is the mask')

# now we apply the mask to the Data object
air._apply_mask(mask)  # applies a mask to each timestep
f3 = map_plot(air, ax=ax3, title='Masked data')

# if you want you can estimate automatically the bounding box
air.cut_bounding_box()
f4 = map_plot(air, ax=ax4, title='Cutted data')

plt.show()



