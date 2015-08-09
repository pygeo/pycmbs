# -*- coding: utf-8 -*-
"""
This file is part of pyCMBS.
(c) 2012- Alexander Loew
For COPYING and LICENSE details, please refer to the LICENSE file
"""


class Koeppen(object):
    """
    KOEPPEN CLASS
    class to generate koeppen plot
    """

    def __init__(self, temp=None, precip=None, lsm=None):
        """
        Koeppen class
        This class implements the functionality to generate koeppen plots.

        Parameters
        ----------
        temp : Data
            data objekt of temperature
        precip : Data
            data objekt of precipitation
        lsm : Data
            data objekt of land-sea-mask (0.0 to 1.0)

        EXAMPLES
        ========

        """

        # check consistency
        if temp is None:
            raise ValueError('No temperature given')
        if precip is None:
            raise ValueError('No precipitation given')
        if lsm is None:
            raise ValueError('No land-sea-mask given')

        # set values of class
        self.temp = temp
        self.precip = precip
        self.lsm = lsm

        if not self._check_resolution():
            raise ValueError('ERROR:The three array differe in the resolution')
        if not self._check_units():
            raise ValueError('ERROR:The units of one value is wrong')

        # Create new koeppen Color map
        self.koeppen_cmap()
        self.cmap = cm.get_cmap('koeppen')
        # convert from [kg m-2 s-1] to [kg m-2 day-1] (= [mm day-1])
         # ??? Unklar warum nicht 'precip.mulc(60. * 60. * 24. * 365.)'
        self.precip = precip.mulc(60. * 60. * 24. * 365. / 12., copy=True)
        self.temp = temp.subc(273.15, copy=True)  # ??? Unklar warum nicht 'temp.subc(273.15)'

        Psum = self.precip.timsum(return_object=True)            # Berechnet die Summe der Jahresniederschlag

        nt, ny, nx = self.temp.shape
        nlat = ny
        nlon = nx

        Pmin = self.precip.data.min(axis=0)
        Pmax = self.precip.data.max(axis=0)

        precipHS = self.precip.copy()
        precipHS.data[(0, 1, 2, 3, 4, 5), 0:(nlat / 2 - 1), :] \
            = self.precip.data[(3, 4, 5, 6, 7, 8), 0:(nlat / 2 - 1), :]
        precipHS.data[(6, 7, 8, 9, 10, 11), 0:(nlat / 2 - 1), :] \
            = self.precip.data[(3, 4, 5, 6, 7, 8), 0:(nlat / 2 - 1), :]
        precipHS.data[(0, 1, 2, 3, 4, 5), (nlat / 2):(nlat - 1), :] \
            = self.precip.data[(0, 1, 2, 9, 10, 11), (nlat / 2):(nlat - 1), :]
        precipHS.data[(6, 7, 8, 9, 10, 11), (nlat / 2):(nlat - 1), :] \
            = self.precip.data[(0, 1, 2, 9, 10, 11), (nlat / 2):(nlat - 1), :]

        precipHW = self.precip.copy()
        precipHW.data[(0, 1, 2, 3, 4, 5), 0:(nlat / 2 - 1), :] = self.precip.data[(0, 1, 2, 9, 10, 11), 0:(nlat / 2 - 1), :]
        precipHW.data[(6, 7, 8, 9, 10, 11), 0:(nlat / 2 - 1), :] = self.precip.data[(0, 1, 2, 9, 10, 11), 0:(nlat / 2 - 1), :]
        precipHW.data[(0, 1, 2, 3, 4, 5), (nlat / 2):(nlat - 1), :] = self.precip.data[(3, 4, 5, 6, 7, 8), (nlat / 2):(nlat - 1), :]
        precipHW.data[(6, 7, 8, 9, 10, 11), (nlat / 2):(nlat - 1), :] = self.precip.data[(3, 4, 5, 6, 7, 8), (nlat / 2):(nlat - 1), :]

        PminHS = precipHS.data.min(axis=0)   # Bestimmt den minimalen Monastniederschlag aus PmaxHS
        PmaxHS = precipHS.data.max(axis=0)   # Bestimmt den maximalen Monastniederschlag aus PmaxHS
        PminHW = precipHW.data.min(axis=0)   # Bestimmt den minimalen Monastniederschlag aus PminHW
        PmaxHW = precipHW.data.max(axis=0)   # Bestimmt den maximalen Monastniederschlag aus PminHW

        Tavg = self.temp.data.mean(axis=0)   # Bestimmt die mittlere Jahrestemperatur
        Tmin = self.temp.data.min(axis=0)     # Bestimmt die minimale Monatstemperatur
        Tmax = self.temp.data.max(axis=0)     # Bestimmt die maximale Jahrestemperatur

        self.Clim = self.precip.timmean(return_object=True)
        self.Clim.units = "climate type"

        for lat in range(0, nlat):
            for lon in range(0, nlon):
                psum = Psum.data.data[lat][lon]
                pmin = Pmin[lat][lon]
                pminhs = PminHS[lat][lon]
                pminhw = PminHW[lat][lon]
                pmaxhs = PmaxHS[lat][lon]
                pmaxhw = PmaxHW[lat][lon]
                tavg = Tavg[lat][lon]
                tmin = Tmin[lat][lon]
                tmax = Tmax[lat][lon]
                self.Clim.data.data[lat][lon] = self.set_clim(psum, pmin, pminhs, pminhw, pmaxhs, pmaxhw, tavg, tmin, tmax)

        self.Clim.data.mask[less(self.lsm.data, 0.5)] = True

    def koeppen_cmap(self):
        """
        Create a colormap with 14 discrete colors and register it
        """
        # define individual colors as hex values
        cpool = ['#7f0000', '#ff0000', '#ff4c4c', '#ff9999', '#ffa500',
                 '#ffff4c', '#009900', '#00ff00', '#99ff99', '#990099',
                 '#e500e5', '#ff66ff', '#0000ff', '#9999ff', '#000000']
        cmap3 = col.ListedColormap(cpool[0:14], 'koeppen')
#       plt.cm.register_cmap(cmap=cmap3,name='koeppen',lut=15)
        plt.cm.register_cmap(cmap=cmap3, name='koeppen')
        return cmap3

    def set_clim(self, psum, pmin, pminhs, pminhw, pmaxhs, pmaxhw, tavg, tmin, tmax):
        clim = -999

        if tmin > 18:
            if pmin > 60:                 # A(B)
                clim = 1                    # Af
            else:
                if pmin > (0.04 * (2500 - psum)):      # A(B)-msw
                    clim = 2                  # Am
                else:
                    if (pminhs < 40) and (pminhs < (pmaxhw / 3)):   # A(B)-sw
                        if (psum / 10) < (2 * tavg):          # A(B)-s
                            if (psum / 10) < (tavg):            # B
                                clim = 6                    # BW
                            else:
                                clim = 5                    # BS
                        else:
                            clim = 3                      # As
                    else:
                        if (psum / 10) < (2 * (tavg + 14)):       # A(B)-w
                            if (psum / 10) < (tavg + 14):       # B
                                clim = 6                    # BW
                            else:
                                clim = 5                    # BS
                        else:
                            clim = 4                      # Aw
        else:
            if (pminhs < 40) and (pminhs < (pmaxhw / 3)):   # CDE(B)
                if (psum / 10) < (2 * tavg):          # CDE(B)-s
                    if (psum / 10) < (tavg):            # B
                        clim = 6                    # BW
                    else:
                        clim = 5                    # BS
                else:
                    if tmax < 10:                # CDE-s
                        if tmax < 0:                # E
                            clim = 14                 # EF
                        else:
                            clim = 13                 # ET
                    else:
                        if (tmin > -3):             # CD-s
                            clim = 8                  # Cs
                        else:
                            clim = 11                 # Ds
            else:
                if pminhw < (pmaxhs / 10):            # CDE(B)-fw
                    if (psum / 10) < (2 * (tavg + 14)):     # CDE(B)-w
                        if (psum / 10) < (tavg + 14):         # B
                            clim = 6                  # BW
                        else:
                            clim = 5                  # BS
                    else:
                        if tmax < 10:               # CDE-w

                            if (tmax < 0):                # E
                                clim = 14               # EF
                            else:
                                clim = 13              # ET
                        else:
                            if (tmin > -3):               # CD-w
                                clim = 9                # Cw
                            else:
                                clim = 12               # Dw
                else:
                    if (psum / 10) < (2 * (tavg + 7)):      # CDE(B)-f
                        if (psum / 10) < (tavg + 7):          # B
                            clim = 6                  # BW
                        else:
                            clim = 5                  # BS
                    else:
                        if (tmax < 10):             # CDE-f
                            if (tmax < 0):                # E
                                clim = 14               # EF
                            else:
                                clim = 13              # ET
                        else:
                            if (tmin > -3):               # CD-f
                                clim = 7                # Cf
                            else:
                                clim = 10              # Df
        return clim

    def _check_resolution(self):
        """
        This routine just checks if all three array have a equal number of ny and nx values
        """
        nt_t, ny_t, nx_t = self.temp.shape
        nt_p, ny_p, nx_p = self.precip.shape
        ny_l, nx_l = self.lsm.shape

        if (ny_t != ny_p) or (ny_t != ny_l):
            sys.exit('ERROR: The resolution ot the three arrays differ in \
       Y-dimension: \n' + str(ny_t) + "(temp)  " + str(ny_p)
                     + "(precip) " + str(ny_l) + "(lsm) ")
            return False

        if (nx_t != nx_p) or (nx_t != nx_l):
            sys.exit('ERROR: The resolution ot the three arrays differ in \
       X-dimension: \n' + str(nx_t) + "(temp)  " + str(nx_p)
                     + "(precip) " + str(nx_l) + "(lsm) ")
            return False

        return True

    def _check_units(self):
        """
        This routine just checks if all three array have a equal number of ny and nx values
        """
        if self.precip.unit != "kg/m^2s":
            raise ValueError('ERROR: The unit of the precip is not [kg/m^2s] its set to [' + self.precip.unit + "]")

        if self.temp.unit != "K":
            raise ValueError('ERROR: The unit of the temperature is not [K] its set to [' + self.temp.unit + "]")

        if self.lsm.unit != "fractional":
            raise ValueError('ERROR: The unit of the temperature is not [fractional] its set to [' + self.temp.unit + "]")

        return True

    def copy(self):
        """
        Returns the Clim Data as an Data variable
        """
        return self.Clim

    def climfrac(self):
        """
        This routine calculats the fraction of each type in per centum.
        ToDo:
        Unclear id the print is OK or if the values should given back as an array.
        """
        climfrac = [0] * 14

        ny, nx = self.Clim.data.data.shape
        for ny in range(0, ny - 1):
            Aweight = cos(self.Clim.lat[ny][0] / 180 * 3.14159265359)
            for nx in range(0, nx - 1):
                clim = int(self.Clim.data.data[ny][nx])
                climfrac[clim - 1] = climfrac[clim - 1] + Aweight

        s = sum(climfrac)
        climfrac[:] = [x / s for x in climfrac]

        print "Af: " + str(climfrac[0])
        print "Am: " + str(climfrac[1])
        print "As: " + str(climfrac[2])
        print "Aw: " + str(climfrac[3])
        print "BS: " + str(climfrac[4])
        print "BW: " + str(climfrac[5])
        print "Cf: " + str(climfrac[6])
        print "Cs: " + str(climfrac[7])
        print "Cw: " + str(climfrac[8])
        print "Df: " + str(climfrac[9])
        print "Ds: " + str(climfrac[10])
        print "Dw: " + str(climfrac[11])
        print "ET: " + str(climfrac[12])
        print "EF: " + str(climfrac[13])

    def legend(self):
        """
        This routine prints a legend of the geiger-koeppen types.
        The description is taken from:
        MARKUS KOTTEK, JUERGEN GRIESER, CHRISTOPH BECK , BRUNO RUDOLF and FRANZ RUBEL
        World Map of the Koeppen-Geiger climate classification updated
        Meteorologische Zeitschrift, Vol. 15, No. 3, 259-263 (June 2006)

        """

        print "|================= Class legend =================|"
        print "| Af: Equatorial rainforest, fully humid         |"
        print "| Am: Equatorial climates                        |"
        print "| As: Equatorial monsoon                         |"
        print "| Aw: Equatorial savannah with dry winter        |"
        print "|------------------------------------------------|"
        print "| BS: Steppe climate                             |"
        print "| BW: Desert climate                             |"
        print "|------------------------------------------------|"
        print "| Cf: Warm temperate climate, fully humid        |"
        print "| Cs: Warm temperate climate with dry summer     |"
        print "| Cw: Warm temperate climate with dry winter     |"
        print "|------------------------------------------------|"
        print "| Df: Snow climate, fully humid                  |"
        print "| Ds: Snow climate with dry summer               |"
        print "| Dw: Snow climate with dry winter               |"
        print "|------------------------------------------------|"
        print "| ET: Tundra climate                             |"
        print "| EF: Frost climate                              |"
        print "|================================================|"

    def plot(self, **kwargs):
        """
        This routine plots the data of his own geiger-koeppen data by
        using the plot-routine map_plot.
        It use the own created color-map and sets the color-bar to a
        horizontal orientation.
        It set the range of values between 0.5 and 14.5. Which are the
        possible values of geiger-koeppen.
        ToDo:
        At the moment the label of the geiger-koeppen types are missing
        at the color-bar
        """
        map_plot(self.Clim, cmap_data=self.cmap, colorbar_orientation='horizontal', vmin=0.5, vmax=14.5,
                 cticks=[1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 11., 12., 13., 14.],
                 cticklabels=["Af", "Am", "As", "Aw", "BS", "BW", "Cf",
                              "Cs", "Cw", "Df", "Ds", "Dw", "ET", "EF"],
                 nclasses=15, **kwargs)
