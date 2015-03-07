# -*- coding: utf-8 -*-

from pycmbs.plots import GlecklerPlot
import matplotlib.pyplot as plt


# this is just an artificial example illustrating principle usage of the
# Portraet diagram

fig = plt.figure(figsize=(8,6))

G = GlecklerPlot(fig=fig)
# register models
G.add_model('MPI-ESM-MR')
G.add_model('MPI-ESM-LR')
G.add_model('MPI-ESM-P')

# then register variables
G.add_variable('Ta')
G.add_variable('P')

# after that you can add values to be plotted
# pos=1 top triangle, pos2= lower triangle
# different positions are for different observations
G.add_data('Ta','MPI-ESM-MR',0.5,pos=1)
G.add_data('Ta','MPI-ESM-MR',-0.2,pos=2)

G.add_data('Ta','MPI-ESM-LR',0.3,pos=1)
G.add_data('Ta','MPI-ESM-LR',0.2,pos=2)

G.add_data('Ta','MPI-ESM-P',0.05,pos=1)
G.add_data('Ta','MPI-ESM-P',-0.1,pos=2)

G.add_data('P','MPI-ESM-P',0.5,pos=1)
G.add_data('P','MPI-ESM-P',-0.2,pos=2)

G.add_data('P','MPI-ESM-LR',0.3,pos=1)

G.add_data('P','MPI-ESM-MR',-0.2,pos=2)

G.plot() #do plot

plt.show()

