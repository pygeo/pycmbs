# -*- coding: utf-8 -*-
import numpy as np
from matplotlib import pyplot as plt

f = plt.figure()
ax = f.add_subplot(111)

x = np.random.random(100)
ax.plot(x)
ax.set_title('This is a test plot')

plt.show()


