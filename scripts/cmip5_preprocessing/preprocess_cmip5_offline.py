"""
This file is part of pyCMBS. (c) 2012-2014
For COPYING and LICENSE details, please refer to the file
COPYRIGHT.md
"""

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

def main():
    from pycmbs.benchmarking import preprocessor
    import datetime as dt
    import sys

    if len(sys.argv) == 3:
        the_experiment = sys.argv[1]
        the_variable = sys.argv[2]
    else:
        # specifiy experiment here
        the_experiment = 'historical'
        the_variable = 'rsds'

    print('Experiment: %s' % the_experiment)
    print('Variable: %s' % the_variable)

    output_dir = '/data/share/mpiles/TRS/PROJECT_RESULTS/EvaClimod/CMIP5_RAWDATA_NEW/radiation/dummy_out/' + the_experiment + '_' + the_variable + '/'
    data_dir = '/data/share/mpiles/TRS/PROJECT_RESULTS/EvaClimod/CMIP5_RAWDATA_NEW/radiation/'

    # init parser that returns a list of institutes and models
    CP = preprocessor.CMIP5ModelParser(data_dir)
    model_list = CP.get_all_models()

    # perform for each institute and model the calculation of ensemble means etc.
    for institute in model_list.keys():
        for model in model_list[institute]:
            output_file = output_dir + the_variable + '_Amon_' + model + '_' + the_experiment + '_ensmean.nc'
            E = preprocessor.CMIP5Preprocessor(data_dir, output_file, the_variable, model, the_experiment, institute=institute)
            E.get_ensemble_files()
            E.ensemble_mean(delete=False, start_time=dt.datetime(1979,1,1), stop_time=dt.datetime(2012,12,31))

if __name__ == '__main__':
    main()
