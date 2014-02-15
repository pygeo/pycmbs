"""
preprocess CMIP5 data offline
The data is assumed to have been extracted already from the archive
while the archive structure from the MiKlip server has been preseverd
except that the VERSION directory has been removed

This script does
a) merge data from different timeslices into a single timesclice
b) extract data for specified timeperiods
c) calculate ensemble means and stdv
"""

from pycmbs.benchmarking import preprocessor
import datetime as dt

output_dir = '/data/share/mpiles/TRS/PROJECT_RESULTS/EvaClimod/CMIP5_RAWDATA_NEW/radiation/dummy_out/amip_rsus/'
data_dir = '/data/share/mpiles/TRS/PROJECT_RESULTS/EvaClimod/CMIP5_RAWDATA_NEW/radiation/'

# specifiy experiment here
the_experiment = 'amip'
the_variable = 'rsds'

# init parser that returns a list of institutes and models
CP = CMIP5ModelParser(data_dir)
model_list = CP.get_all_models()

# perform for each institute and model the calculation of ensemble means etc.
for institute in model_list.keys():
    for model in model_list[institute]:
        output_file = output_dir + the_variable + '_Amon_' + model + '_' + the_experiment + '_ensmean.nc'
        E = CMIP5Preprocessor(data_dir, output_file, the_variable, model, the_experiment, institute=institute)
        E.get_ensemble_files()
        E.ensemble_mean(delete=False, start_time=dt.datetime(1979,1,1), stop_time=dt.datetime(2012,12,31))
