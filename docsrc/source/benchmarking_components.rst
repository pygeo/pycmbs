Benchmarking components
-----------------------

The benchmarking part of pyCMBS is based on a modular system of diagnostics.
Currently the focus is on the comparison of climate mean states between
observations and models. Below are examples of the currently implemented
components which the user can combine to generate a specific report for a
certain variable.

Map season
~~~~~~~~~~

Creates a figure with seasonal or monthly climatological means of the
investigated variables.

.. plot:: ./figures/fig_map_season.py



Map difference
~~~~~~~~~~~~~~

Create a plot of temporal mean fields of observations and models and their
absolute and relative differences.

.. plot:: ./figures/fig_map_difference.py

Hovmoeller plot
~~~~~~~~~~~~~~~

Generates a standard Hovmoeller plot (= time-latitude plot) of the variable.

.. plot:: ./figures/fig_hovmoeller.py



Pattern correlation
~~~~~~~~~~~~~~~~~~~

The correlation between the spatial patterns of investigated variables between
models and observations can be vizualized in differnt ways. For each month or
season the spatial correlation coefficient is estimated and vizualized as a
timeline.

.. plot:: ./figures/fig_pattern_correlation.py


Portraet diagram
~~~~~~~~~~~~~~~~

The *Portraet diagram* was proposed by `Gleckler et al. (2008) <http://www.agu.org/pubs/crossref/2008/2007JD008972.shtml>`_. It is an
efficient way to vizualize the relative rank of different models  for
comparisons against different observations. While the original Portraet diagram
supports only two different observations, pyCMBS supports up to four different
datasets for each variable.

.. plot:: ./figures/fig_portraet.py


Global Mean Plot
~~~~~~~~~~~~~~~~

Global mean timeseries of a specific variable. The global mean values are
estimated using proper area weighting. The plot object also allows to
automatically plot climatology mean values.

.. plot:: ./figures/fig_globalmean.py




