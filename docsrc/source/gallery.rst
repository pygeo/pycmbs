==============
pyCMBS gallery
==============

In the following, examples will be given that introduce different features of pyCMBS.

Basic plotting
--------------

.. plot:: ../../pycmbs/examples/01_basic_plotting.py
  :include-source:

For plotting you have a rich suite of keyword parameters. Please use the
python help() system to get full documentation from docstrings. Some
more illustrations of options are provided::

    map_plot(air,show_timeseries=True, use_basemap=True,title='show_timeseries=True')
    map_plot(air,show_zonal=True, use_basemap=True,title='show_zonal=True')
    map_plot(air,show_histogram=True, use_basemap=True,title='show_histogram=True')

.. plot:: ../../pycmbs/examples/01_basic_plotting_A.py

And a few more details on customizing your map ...::

    map_plot(air, use_basemap=True, title='vmin=-30.,vmax=30.,cmap_data=RdBu_r', vmin=-30., vmax=30., cmap_data='RdBu_r', ax=ax1)
    map_plot(air, contours=True, use_basemap=True, levels=np.arange(-50.,50.,10.), title='contours=True', ax=ax2)
    map_plot(air, show_stat=True, use_basemap=True,title='show_stat=True',ax=ax3)
    map_plot(air, show_stat=True, stat_type='median', use_basemap=True, title='show_stat=True,stat_type="median"', ax=ax4)

.. plot:: ../../pycmbs/examples/01_basic_plotting_B.py


Basic data analysis
-------------------

TBD

Working with different datasets
-------------------------------

TBD




