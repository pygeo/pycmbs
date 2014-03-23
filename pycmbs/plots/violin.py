# -*- coding: utf-8 -*-
"""
This file is part of pyCMBS. (c) 2012-2014
For COPYING and LICENSE details, please refer to the file
COPYRIGHT.md
"""
import matplotlib.pyplot as plt
import numpy as np


class ViolinPlot(object):
    """
    generate ViolinPlot

    References
    ----------
    inspired by
    [1] http://pyinsci.blogspot.de/2009/09/violin-plot-with-matplotlib.html
    [2] http://en.wikipedia.org/wiki/Violin_plot
    [3] http://statsmodels.sourceforge.net/devel/generated/statsmodels.graphics.boxplots.violinplot.html
    [4] http://nbviewer.ipython.org/github/EnricoGiampieri/dataplot/blob/master/statplot.ipynb
    """

    def __init__(self, data, data2=None, labels=None, ax=None,
                 boxplot=True, figsize=(10, 6)):
        """
        Parameters
        ----------
        data : ndarray/dict
            data to be plotted. needs to be of geometry [ngroups, nsamp]
            where ngroups is the the number of different groups to
            be plotted and would correspond also to len(labels)
        data2 : ndarray/dict
            second dataset to be plotted
        labels : list
            labels to be used for xticks
        ax : axis
            when provided, then the plot is generated on the particular
            axis
        boxplot : bool
            plot boxplot
        figsize : tuple
            figure size
        """
        self.data = data
        self.data2 = data2
        if labels is None:
            self.labels = range(len(self.data))
        else:
            self.labels = labels

        if ax is None:
            fig = plt.figure(figsize=figsize)
            rect = [0.1, 0.2, 0.8, 0.7]  # l,b,w,h
            #~ self.ax = fig.add_subplot(1, 1, 1)
            self.ax = fig.add_axes(rect)
        else:
            self.ax = ax
        self.boxplot = boxplot

        # check
        self._check()

    def _check(self):
        """
        routine to check internal consistency
        """
        if self.data is not None:
            if len(self.labels) != len(self.data):
                raise ValueError('Invalid geometry of labels and data!')
        if self.data2 is not None:
            if len(self.data) != len(self.data2):
                raise ValueError('Data arrays need to have same geometry')

    def plot(self, alpha=0.3, classic=False):
        """
        plot ViolinPlot

        Parameters
        ----------
        alpha : float
            alpha value for area fill
        classic : bool
            make classic violin plot
        """
        if self.data is None:
            raise ValueError('Data is None and can therefore not be plotted!')
        if classic:
            self._plot_classic(alpha=alpha)
        else:
            self._plot_two_sided()
        self._set_xticks()

    def _plot_half_violin(self, data, pos, left=False, **kwargs):
        """
        plot half violin
        inspired by [4]

        Parameters
        ----------
        data : ndarray
            array with data (only valid data, no masked array)
        pos : ndarray
            position
        left : bool
            specifies if the plot should be on the left side
        """
        from scipy.stats import gaussian_kde
        amplitude = kwargs.pop('amplitude', 0.33)
        x = np.linspace(min(data), max(data), 101)
        v = gaussian_kde(data).evaluate(x)
        v = v/v.max()*amplitude * (1 if left else -1)
        kwargs.setdefault('facecolor', 'r')
        kwargs.setdefault('alpha', 0.33)
        return self.ax.fill_betweenx(x, pos, pos+v, **kwargs)

    def _plot_two_sided(self, color1='b', color2='b'):
        """
        violin plot with two sides
        inspired by [4]

        Parameters
        ----------
        color1 : str
            color for left plot
        color2 : str
            color for right plot
        """
        positions = self._get_positions()
        data2 = self.data2 if self.data2 is not None else self.data
        labels = self.labels if self.labels is not None else positions
        for pos, key in zip(positions, labels):
            try:
                # in case that data is a dict
                d1, d2 = self.data[key], data2[key]
            except TypeError:
                d1, d2 = self.data[pos], data2[pos]
            if self.data is not data2:
                color2 = 'r'

            # generate plot
            self._plot_half_violin(d1, pos, left=False, facecolor=color1)
            self._plot_half_violin(d2, pos, left=True, facecolor=color2)

            # division line between the two half
            self.ax.plot([pos]*2, [min(min(d1), min(d2)),
                         max(max(d1), max(d2))], '-', color='grey')

    def _set_xticks(self, rotation=30.):
        """
        set ticklabels
        """
        self.ax.set_xticks(self._get_positions())
        self.ax.set_xticklabels(self.labels, rotation=rotation)

    def _get_positions(self):
        """
        get positions of xticks
        """
        return range(len(self.data))

    def _plot_classic(self, alpha):
        """
        create classical violin plots on an axis
        http://pyinsci.blogspot.de/2009/09/violin-plot-with-matplotlib.html

        Parameters
        ----------
        alpha : float
            alpha value for area fill
        """
        from scipy.stats import gaussian_kde
        pos = self._get_positions()
        dist = max(pos)-min(pos)
        w = min(0.15*max(dist, 1.0), 0.5)
        for d, p in zip(self.data, pos):
            if not np.all(d==0.):  # avoid singular matrices
                k = gaussian_kde(d)  # calculates the kernel density
                m = k.dataset.min()  # lower bound of violin
                M = k.dataset.max()  # upper bound of violin
                x = np.arange(m, M, (M-m)/100.)  # support for violin
                v = k.evaluate(x)  # violin profile (density curve)
                v = v/v.max()*w  # scaling the violin to the available space
                self.ax.fill_betweenx(x, p, v+p, facecolor='y', alpha=alpha)
                self.ax.fill_betweenx(x, p, -v+p, facecolor='y', alpha=alpha)
        if self.boxplot:
            self.ax.boxplot(self.data, notch=1, positions=pos,
                            vert=True, sym='')


def _classic_example():
    """
    some example how to do violin plotting
    """
    plt.close('all')
    pos = range(5)
    data = [np.random.normal(size=100) for i in pos]

    V = ViolinPlot(data)
    V.plot(classic=True)

    data = [np.random.normal(size=100) for i in pos]
    V1 = ViolinPlot(data, labels=['A', 'B', 'C', 'D', 'E'])
    V1.plot()

    data = [np.random.normal(size=100) for i in pos]
    data2 = [np.random.normal(size=100) for i in pos]
    V2 = ViolinPlot(data, data2=data2, labels=['A', 'B', 'C', 'D', 'E'])
    V2.plot()

    plt.show()

if __name__ == '__main__':
    _classic_example()
