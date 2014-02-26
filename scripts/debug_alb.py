"""
This file is part of pyCMBS. (c) 2012-2014
For COPYING and LICENSE details, please refer to the files
LICENSE.md and COPYRIGHT.md
"""

from pycmbs import *
from pylab import *


def get_data_pool_directory():
    if 'SEP' in os.environ.keys():
        data_pool_directory = os.environ['SEP'] #get directory of pool/SEP
    else:
        data_pool_directory = '/pool/SEP/'

    return data_pool_directory

def get_T63_landseamask(shift_lon,mask_antarctica=True):
    """
    get JSBACH T63 land sea mask
    the LS mask is read from the JSBACH init file

    @todo: put this to the JSBACH model class
    """
    ls_file = get_data_pool_directory() + 'variables/land/land_sea_mask/jsbach_T63_GR15_4tiles_1992.nc'
    ls_mask = Data(ls_file,'slm',read=True,label='T63 land-sea mask',lat_name='lat',lon_name='lon',shift_lon=shift_lon)
    msk=ls_mask.data>0.; ls_mask.data[~msk] = 0.; ls_mask.data[msk] = 1.
    ls_mask.data = ls_mask.data.astype('bool') #convert to bool
    if mask_antarctica:
        ls_mask.data[ls_mask.lat < -60.] = False

    return ls_mask

close('all')

mo = Data('/home/m300028/shared/data/SEP/variables/land/surface_albedo/modis/with_snow/monmean.nc','surface_albedo_WSA',read=True)
ce = Data('/home/m300028/shared/dev/svn/pyCMBS/framework/tmp_processing/CERES_EBAF-Surface__Ed2.6r__sfc_sw_up_all_mon__1x1__200003-201002_albedo_monmean.nc','sfc_sw_up_all_mon',read=True)
mo_org = Data('/home/m300028/shared/data/SEP/variables/land/surface_albedo/modis/with_snow/T63_MCD43C3-QC_merged.nc','surface_albedo_WSA',read=True)
me = Data('/home/m300028/shared/data/SEP/variables/land/surface_albedo/meteosat/METEOSAT_MSA_BHR_GAUSS_work.nc','bhr_mean_all',read=True)
me.data = np.ma.array(me.data.data,mask=me.data.data < 0.)


#jsb = Data('/home/m300028/shared/dev/svn/pyCMBS/framework/tmp_processing/rsus_Amon_MPI-ESM-LR_amip_ensmean_1983-01-01_1999-12-31_T63_monmean.nc','sfc_sw_up_all_mon',read=True)

mo.adjust_time(day=15)

me.time -=  1.#if we assume that the albedo values are wrong by half a month!!! --> albedo study !!!
me.adjust_time(day=15)


msk=get_T63_landseamask(True)



#l/s mask
mo._apply_mask(msk)
mo_org._apply_mask(msk)
ce._apply_mask(msk)
me._apply_mask(msk)




#~ map_plot(ce)
#~ map_plot(mo)
#~ map_plot(mo_org,title='MODIS original')

#global means
figure()
plot(num2date(ce.time),ce.fldmean(),label='ceres-monthly')
plot(num2date(mo.time),mo.fldmean(),label='modis-monthly')
plot(num2date(mo_org.time),mo_org.fldmean(),label='modis-original')
grid()
legend()


sahel_all    = Region(-20.,45.,10.,13.,'sahel_all',type='latlon') #specify a region with coordinates
amazon = Region(-90.,-30,-25.,10.,'amazon',type='latlon') #specify a region with coordinates


sub = amazon

mo .get_aoi_lat_lon(sub)
mo_org .get_aoi_lat_lon(sub)
ce .get_aoi_lat_lon(sub)
me .get_aoi_lat_lon(sub)

map_plot(ce,use_basemap=True)

#map_plot(ce.sub(mo))



figure()
plot(num2date(ce.time),ce.fldmean(),label='ceres-monthly')
plot(num2date(mo.time),mo.fldmean(),label='modis-monthly')
plot(num2date(mo_org.time),mo_org.fldmean(),label='modis-original')
plot(num2date(me.time),me.fldmean(),label='Meteosat')
title(sub.label)
grid()
legend()

show()


#1) original MODIS Daten fldmean() rechnen
#2) unmaskierte MODIS Daten ???
#3) sind die Meteosat Albedodaten geshiftet ???? --> originaldaten anschauen ???

#plot for amazon !!!


#- upward flux in the paper is weighted with MPI-ESM downward radiation flux
#- phase shift in CERES/MODIS data?
#- amazon: CERES is darker by 0.02 on average
#- CERES shows rather different shape of phenology over the amazon
#- all sky vs. clear sky

