"""
Generate a dummy configuration file that can be used
for pycmbs processing of CMIP5 arechive data

The script needs the path of the rootdirectory where the data is located.
It then creates a list of model items like required by the *.cfg file.

e.g.:
BCC:bcc-csm1-1,CMIP5RAW,amip,/data/share/mpiles/.../PATH/TO/RAWDATA
"""

from pycmbs.benchmarking import preprocessor
import os

def main():

    ofile = 'test.txt'
    if os.path.exists(ofile):
        os.remove(ofile)
    o = open(ofile, 'w')

    data_dir = '/data/share/mpiles/TRS/PROJECT_RESULTS/EvaClimod/CMIP5_RAWDATA_NEW/radiation'
    if data_dir[-1] != os.sep:
        data_dir += os.sep
    CMP = preprocessor.CMIP5ModelParser(data_dir)
    model_list = CMP.get_all_models()

    experiments = ['amip', 'historical']

    for experiment in experiments:
        o.write('\n')
        o.write('##############################\n')
        o.write('# ' + experiment + '\n')
        o.write('##############################\n')
        for institute in model_list.keys():
            for model in model_list[institute]:
                s = institute + ':' + model + ',' + 'CMIP5RAW' + ',' + experiment + ',' + data_dir + '\n'
                o.write(s)
    o.close()





if __name__ == "__main__":
    main()
