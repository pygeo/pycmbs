# implicit import
# from pyCMBS import *
import matplotlib.pyplot as plt

plt.close('all')

#//// read data ////
gridfile ='../example_data/icon/r2b4_amip.nc'
datafile = '../example_data/icon/rms0006_atm_phy_DOM01_ML_0001.nc'

#/// read data ///
IC = Icon(datafile,gridfile,'rsns'); IC.read() #todo: read as standard for Data
IC.label='shortwave net flux at surface'; IC.unit='$W/^2$'

#/// Do plotting ///
print('Doing first plot ...')
map_plot(IC,use_basemap=True)
print('Doing second plot ...')
map_plot(IC,use_basemap=False)

plt.show()
