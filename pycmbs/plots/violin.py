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

    def __init__(self, data, labels=None, ax=None, boxplot=True, figsize=(10,10)):
        """
        Parameters
        ----------
        data : ndarray
            data to be plotted. needs to be of geometry [ngroups, nsamp]
            where ngroups is the the number of different groups to
            be plotted and would correspond also to len(labels)
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
        """ routine to check internal consistency """
        if self.data is not None:
            if len(self.labels) != len(self.data):
                raise ValueError('Invalid geometry of labels and data!')

    def plot(self, alpha=0.3):
        """
        plot Violinplot
        """
        if self.data is None:
            raise ValueError('Data is None and can therefore not be plotted!')
        ### TODO implement also plot like shown in Ref [4]
        self._plot_classic(alpha=alpha)
        self._set_xticks()

    def _set_xticks(self, rotation=90.):
        """ set ticklabels """
        self.ax.set_xticks(self._get_positions())
        self.ax.set_xticklabels(self.labels, rotation=rotation)

    def _get_positions(self):
        """ get positions of xticks """
        return range(len(self.data))

    def _plot_classic(self, alpha):
        """
        create classical violin plots on an axis
        http://pyinsci.blogspot.de/2009/09/violin-plot-with-matplotlib.html
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
            self.ax.boxplot(self.data, notch=1, positions=pos, vert=True, sym='')


def _classic_example():
    """
    some example how to do violin plotting
    """
    plt.close('all')
    pos = range(5)
    data = [np.random.normal(size=100) for i in pos]

    V = ViolinPlot(data)
    V.plot()

    data = [np.random.normal(size=100) for i in pos]
    V1 = ViolinPlot(data, labels=['A', 'B', 'C', 'D', 'E'])
    V1.plot()

    plt.show()


if __name__ == '__main__':
    _classic_example()
