#!/usr/bin/env python
# -*- coding: utf-8 -*-

from pyCMBS import *

from utils import *

class Model(Data):
    '''
    This class is the main class, specifying a climate model or a particular run
    Sub-classes for particular models or experiments are herited from this class
    '''
    def __init__(self,data_dir,dic_variables,name='',**kwargs):
        '''
        constructor for Model class

        INPUT
        -----
        filename: name of the file to read data from (single file currently)
        could be improved later by more abstract class definitions

        dic_variables: dictionary specifiying variable names for a model
        e.g. 'rainfall','var4'
        '''

        #--- set a list with different datasets for different models
        self.dic_vars = dic_variables

        #--- set some metadata
        self.name = name
        self.data_dir = data_dir

        if 'start_time' in kwargs.keys():
            self.start_time = kwargs['start_time']
        else:
            self.start_time = None
        if 'stop_time' in kwargs.keys():
            self.stop_time = kwargs['stop_time']
        else:
            self.stop_time = None


    def get_data(self):
        '''
        central routine to extract data for all variables
        using functions specified in derived class
        '''

        self.variables={}
        for k in self.dic_vars.keys():
            routine = self.dic_vars[k] #get name of routine to perform data extraction
            cmd = 'dat = self.' + routine
            #~ print cmd
            #~ if hasattr(self,routine):
            if hasattr(self,routine[0:routine.index('(')]): #check if routine name is there
                exec(cmd)
                self.variables.update({ k : dat })
            else:
                print 'WARNING: unknown function to read data (skip!) ', routine
                self.variables.update({ k : None })
                #~ sys.exit()




class CMIP5Data(Model):
    def __init__(self,data_dir,model,experiment,dic_variables,name='',shift_lon=False,**kwargs):
        Model.__init__(self,None,dic_variables,name=model,shift_lon=shift_lon,**kwargs)

        self.model = model
        self.experiment = experiment
        self.data_dir = data_dir
        self.shift_lon = shift_lon




    def get_surface_shortwave_radiation_down(self):
        filename1 = self.data_dir + 'rsds/' +  self.model + '/' + 'rsds_Amon_' + self.model + '_' + self.experiment + '_ensmean.nc'



        if self.start_time == None:
            raise ValueError, 'Start time needs to be specified'
        if self.stop_time == None:
            raise ValueError, 'Stop time needs to be specified'

        s_start_time = str(self.start_time)[0:10]
        s_stop_time = str(self.stop_time)[0:10]

        tmp  = pyCDO(filename1,s_start_time,s_stop_time).seldate()
        tmp1 = pyCDO(tmp,s_start_time,s_stop_time).seasmean()
        filename = pyCDO(tmp1,s_start_time,s_stop_time).yseasmean()


        if not os.path.exists(filename):
            return None

        sis = Data(filename,'rsds',read=True,label=self.model,unit='W/m**2',lat_name='lat',lon_name='lon',shift_lon=False)
        print 'Data read!'

        return sis

    def get_surface_shortwave_radiation_up(self):
        filename1 = self.data_dir + 'rsus/' +  self.model + '/' + 'rsus_Amon_' + self.model + '_' + self.experiment + '_ensmean.nc'

        if self.start_time == None:
            raise ValueError, 'Start time needs to be specified'
        if self.stop_time == None:
            raise ValueError, 'Stop time needs to be specified'

        s_start_time = str(self.start_time)[0:10]
        s_stop_time  = str(self.stop_time)[0:10]

        tmp  = pyCDO(filename1,s_start_time,s_stop_time).seldate()
        tmp1 = pyCDO(tmp,s_start_time,s_stop_time).seasmean()
        filename = pyCDO(tmp1,s_start_time,s_stop_time).yseasmean()

        if not os.path.exists(filename):
            return None

        sis = Data(filename,'rsus',read=True,label=self.model,unit='W/m**2',lat_name='lat',lon_name='lon',shift_lon=False)
        print 'Data read!'

        return sis

    def get_albedo_data(self):
        '''
        calculate albedo as ratio of upward and downwelling fluxes
        first the monthly mean fluxes are used to calculate the albedo,
        '''


        #--- read land-sea mask
        ls_mask = get_T63_landseamask(self.shift_lon)

        if self.start_time == None:
            raise ValueError, 'Start time needs to be specified'
        if self.stop_time == None:
            raise ValueError, 'Stop time needs to be specified'

        s_start_time = str(self.start_time)[0:10]
        s_stop_time = str(self.stop_time)[0:10]


        file_down = self.data_dir + 'rsds/' +  self.model + '/' + 'rsds_Amon_' + self.model + '_' + self.experiment + '_ensmean.nc'
        file_up   = self.data_dir + 'rsus/' +  self.model + '/' + 'rsus_Amon_' + self.model + '_' + self.experiment + '_ensmean.nc'

        if not os.path.exists(file_down):
            print 'File not existing: ', file_down
            return None
        if not os.path.exists(file_up):
            print 'File not existing: ', file_up
            return None

        #/// calculate ratio on monthly basis
        # CAUTION: it might happen that latitudes are flipped. Therefore always apply remapcon!

        #select dates
        Fu = pyCDO(file_up,s_start_time,s_stop_time).seldate()
        Fd = pyCDO(file_down,s_start_time,s_stop_time).seldate()

        #remap to T63
        tmpu = pyCDO(Fu,s_start_time,s_stop_time).remap()
        tmpd = pyCDO(Fd,s_start_time,s_stop_time).remap()

        #calculate monthly albedo
        albmon = pyCDO(tmpu,s_start_time,s_stop_time).div(tmpd,output=self.model + '_' + self.experiment + '_albedo_tmp.nc')

        #calculate seasonal mean albedo
        tmp1 = pyCDO(albmon,s_start_time,s_stop_time).seasmean()
        albfile = pyCDO(tmp1,s_start_time,s_stop_time).yseasmean()

        alb = Data(albfile,'rsus',read=True,label=self.model + ' albedo',unit='-',lat_name='lat',lon_name='lon',shift_lon=True)

        alb._apply_mask(ls_mask.data)

        return alb





#####################################################




class JSBACH_BOT(Model):
    '''

    '''

    def __init__(self,filename,dic_variables,experiment,name='',shift_lon=False,**kwargs):

        Model.__init__(self,filename,dic_variables,name=name,**kwargs)
        self.experiment = experiment
        self.shift_lon = shift_lon
        self.get_data()

    def get_albedo_data(self):
        '''
        get albedo data for JSBACH

        returns Data object
        '''

        v = 'var176'

        filename = self.data_dir + 'data/model1/' + self.experiment + '_echam6_BOT_mm_1979-2006_albedo_yseasmean.nc' #todo: proper files
        ls_mask = get_T63_landseamask()

        albedo = Data(filename,v,read=True,
        label='MPI-ESM albedo ' + self.experiment, unit = '-',lat_name='lat',lon_name='lon',
        shift_lon=shift_lon,
        mask=ls_mask.data.data)

        return albedo



    def get_tree_fraction(self):
        '''
        todo implement this for data from a real run !!!
        '''

        ls_mask = get_T63_landseamask()

        filename = '/home/m300028/shared/dev/svn/trstools-0.0.1/lib/python/pyCMBS/framework/external/vegetation_benchmarking/VEGETATION_COVER_BENCHMARKING/example/historical_r1i1p1-LR_1850-2005_forest_shrub.nc'
        v = 'var12'
        tree = Data(filename,v,read=True,
        label='MPI-ESM tree fraction ' + self.experiment, unit = '-',lat_name='lat',lon_name='lon',
        shift_lon=shift_lon,
        mask=ls_mask.data.data,start_time = pl.num2date(pl.datestr2num('2001-01-01')),stop_time=pl.num2date(pl.datestr2num('2001-12-31')))

        return tree

    def get_grass_fraction(self):
        '''
        todo implement this for data from a real run !!!
        '''

        ls_mask = get_T63_landseamask()

        filename = '/home/m300028/shared/dev/svn/trstools-0.0.1/lib/python/pyCMBS/framework/external/vegetation_benchmarking/VEGETATION_COVER_BENCHMARKING/example/historical_r1i1p1-LR_1850-2005_grass_crop_pasture_2001.nc'
        v = 'var12'
        grass = Data(filename,v,read=True,
        label='MPI-ESM tree fraction ' + self.experiment, unit = '-',lat_name='lat',lon_name='lon',
        shift_lon=shift_lon,
        mask=ls_mask.data.data,start_time = pl.num2date(pl.datestr2num('2001-01-01')),stop_time=pl.num2date(pl.datestr2num('2001-12-31')) , squeeze=True  )


        return grass









    def get_surface_shortwave_radiation_down(self,interval = 'season'):
        '''
        get surface shortwave incoming radiation data for JSBACH

        returns Data object
        '''

        v = 'var176'

        y1 = '1979-01-01'; y2 = '2006-12-31'
        rawfilename = self.data_dir + 'data/model/' + self.experiment + '_echam6_BOT_mm_1979-2006_srads.nc'

        if not os.path.exists(rawfilename):
            return None


        #--- read data
        cdo = pyCDO(rawfilename,y1,y2)
        if interval == 'season':
            seasfile = cdo.seasmean(); del cdo
            print 'seasfile: ', seasfile
            cdo = pyCDO(seasfile,y1,y2)
            filename = cdo.yseasmean()
        else:
            raise ValueError, 'Invalid interval option ', interval

        #--- read land-sea mask
        ls_mask = get_T63_landseamask()

        #--- read SIS data
        sis = Data(filename,v,read=True,
        label='MPI-ESM SIS ' + self.experiment, unit = '-',lat_name='lat',lon_name='lon',
        shift_lon=shift_lon,
        mask=ls_mask.data.data)

        return sis


    def get_rainfall_data(self,interval='season'):
        '''
        get rainfall data for JSBACH
        returns Data object
        '''

        #todo: implement preprocessing here

        v = 'var4'
        if interval == 'season':
            filename = self.data_dir + 'data/model/' + self.experiment + '_echam6_BOT_mm_1982-2006_sel_yseasmean.nc'
        else:
            raise ValueError, 'Invalid value for interval: ' + interval

        ls_mask = get_T63_landseamask(self.shift_lon)

        if not os.path.exists(filename):
            stop
            return None

        try: #todo this is silly
            rain = Data(filename,v,read=True,scale_factor = 86400.,
            label='MPI-ESM ' + self.experiment, unit = 'mm/day',lat_name='lat',lon_name='lon',
            shift_lon=self.shift_lon,
            mask=ls_mask.data.data)
        except:
            v='var142'
            rain = Data(filename,v,read=True,scale_factor = 86400.,
            label='MPI-ESM ' + self.experiment, unit = 'mm/day',lat_name='lat',lon_name='lon',
            shift_lon=self.shift_lon,
            mask=ls_mask.data.data)



        return rain



#-----------------------------------------------------------------------

class JSBACH_RAW(Model):
    '''
    RAW JSBACH model output
    '''

    def __init__(self,filename,dic_variables,experiment,name='',shift_lon=False,**kwargs):

        Model.__init__(self,filename,dic_variables,name=name,**kwargs)
        self.experiment = experiment
        self.shift_lon = shift_lon
        self.get_data()

    def get_albedo_data(self):
        '''
        calculate albedo as ratio of upward and downwelling fluxes
        first the monthly mean fluxes are used to calculate the albedo,
        '''

        if self.start_time == None:
            raise ValueError, 'Start time needs to be specified'
        if self.stop_time == None:
            raise ValueError, 'Stop time needs to be specified'




        sw_down = self.get_surface_shortwave_radiation_down()
        sw_up   = self.get_surface_shortwave_radiation_up()
        alb     = sw_up.div(sw_down)
        alb.label = self.experiment + ' albedo'
        alb.unit = '-'




        #~ s_start_time = str(self.start_time)[0:10]
        #~ s_stop_time = str(self.stop_time)[0:10]
#~
        #~ file_down = self.data_dir + 'rsds/' +  self.model + '/' + 'rsds_Amon_' + self.model + '_' + self.experiment + '_ensmean.nc'
        #~ file_up   = self.data_dir + 'rsus/' +  self.model + '/' + 'rsus_Amon_' + self.model + '_' + self.experiment + '_ensmean.nc'
#~
        #~ if not os.path.exists(file_down):
            #~ print 'File not existing: ', file_down
            #~ return None
        #~ if not os.path.exists(file_up):
            #~ print 'File not existing: ', file_up
            #~ return None

        #/// calculate ratio on monthly basis
        # CAUTION: it might happen that latitudes are flipped. Therefore always apply remapcon!

        #select dates
        #~ Fu = pyCDO(file_up,s_start_time,s_stop_time).seldate()
        #~ Fd = pyCDO(file_down,s_start_time,s_stop_time).seldate()

        #remap to T63
        #~ tmpu = pyCDO(Fu,s_start_time,s_stop_time).remap()
        #~ tmpd = pyCDO(Fd,s_start_time,s_stop_time).remap()

        #calculate monthly albedo
        #~ albmon = pyCDO(tmpu,s_start_time,s_stop_time).div(tmpd,output=self.model + '_' + self.experiment + '_albedo_tmp.nc')

        #calculate seasonal mean albedo
        #~ tmp1 = pyCDO(albmon,s_start_time,s_stop_time).seasmean()
        #~ albfile = pyCDO(tmp1,s_start_time,s_stop_time).yseasmean()

        #~ alb = Data(albfile,'rsus',read=True,label=self.model + ' albedo',unit='-',lat_name='lat',lon_name='lon',shift_lon=True)
        return alb




    def get_surface_shortwave_radiation_down(self,interval = 'season'):
        '''
        get surface shortwave incoming radiation data for JSBACH

        returns Data object
        '''

        v = 'swdown_acc'

        y1 = '1992-01-01'; y2 = '2001-12-31'
        rawfilename = self.data_dir + 'yseasmean_' + self.experiment + '_jsbach_' + y1[0:4] + '_' + y2[0:4] + '.nc'

        if not os.path.exists(rawfilename):
            return None


        #--- read data
        #~ cdo = pyCDO(rawfilename,y1,y2)
        #~ if interval == 'season':
            #~ seasfile = cdo.seasmean(); del cdo
            #~ print 'seasfile: ', seasfile
            #~ cdo = pyCDO(seasfile,y1,y2)
            #~ filename = cdo.yseasmean()
        #~ else:
            #~ raise ValueError, 'Invalid interval option ', interval
        filename = rawfilename

        #--- read land-sea mask
        ls_mask = get_T63_landseamask(self.shift_lon)

        #--- read SIS data
        sw_down = Data(filename,v,read=True,
        label=self.experiment + ' ' + v, unit = 'W/m**2',lat_name='lat',lon_name='lon',
        shift_lon=self.shift_lon,
        mask=ls_mask.data.data)

        return sw_down


#-----------------------------------------------------------------------

    def get_surface_shortwave_radiation_up(self,interval = 'season'):
        '''
        get surface shortwave upward radiation data for JSBACH

        returns Data object

        todo CDO preprocessing of seasonal means
        todo temporal aggregation of data --> or leave it to the user!
        '''

        v = 'swdown_reflect_acc'

        y1 = '1992-01-01'; y2 = '2001-12-31'
        rawfilename = self.data_dir + 'yseasmean_' + self.experiment + '_jsbach_' + y1[0:4] + '_' + y2[0:4] + '.nc'

        if not os.path.exists(rawfilename):
            return None


        #--- read data
        #~ cdo = pyCDO(rawfilename,y1,y2)
        #~ if interval == 'season':
            #~ seasfile = cdo.seasmean(); del cdo
            #~ print 'seasfile: ', seasfile
            #~ cdo = pyCDO(seasfile,y1,y2)
            #~ filename = cdo.yseasmean()
        #~ else:
            #~ raise ValueError, 'Invalid interval option ', interval
        filename = rawfilename

        #--- read land-sea mask
        ls_mask = get_T63_landseamask(self.shift_lon)

        #--- read SW up data
        sw_up = Data(filename,v,read=True,
        label=self.experiment + ' ' + v, unit = 'W/m**2',lat_name='lat',lon_name='lon',
        shift_lon=self.shift_lon,
        mask=ls_mask.data.data)

        return sw_up

#-----------------------------------------------------------------------

    def get_rainfall_data(self,interval = 'season'):
        '''
        get surface rainfall data for JSBACH

        returns Data object

        todo CDO preprocessing of seasonal means
        todo temporal aggregation of data --> or leave it to the user!
        '''

        v = 'precip_acc'

        y1 = '1992-01-01'; y2 = '2001-12-31' #todo
        rawfilename = self.data_dir + 'yseasmean_' + self.experiment + '_jsbach_' + y1[0:4] + '_' + y2[0:4] + '.nc'

        if not os.path.exists(rawfilename):
            return None


        #--- read data
        #~ cdo = pyCDO(rawfilename,y1,y2)
        #~ if interval == 'season':
            #~ seasfile = cdo.seasmean(); del cdo
            #~ print 'seasfile: ', seasfile
            #~ cdo = pyCDO(seasfile,y1,y2)
            #~ filename = cdo.yseasmean()
        #~ else:
            #~ raise ValueError, 'Invalid interval option ', interval
        filename = rawfilename

        #--- read land-sea mask
        ls_mask = get_T63_landseamask(self.shift_lon)

        #--- read SW up data
        rain = Data(filename,v,read=True,
        label=self.experiment + ' ' + v, unit = 'mm/day',lat_name='lat',lon_name='lon',
        shift_lon=self.shift_lon,
        mask=ls_mask.data.data,scale_factor = 86400.)

        return rain


#-----------------------------------------------------------------------






