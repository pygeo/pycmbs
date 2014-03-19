"""
This file is part of pyCMBS. (c) 2012-2014
For COPYING and LICENSE details, please refer to the file
COPYRIGHT.md
"""

"""
Combine latex tables that were written using gleckler.write_table routine
"""

import os
from pycmbs.benchmarking.report import Report

# specify directory name
rdir = '/home/m300028/shared/publications/01_ongoing/2014_EvaClimod_SSI/results/'

files = []
files.append({'Albedo' : 'ranking_table_albedo.tex'})
files.append({'SIS' : 'ranking_table_sis.tex'})
files.append({'Rup' : 'ranking_table_surface_upward_flux.tex'})





#|        Albedo       |     RN      |
#| CLARA | MODIS | ... | CLARA | XXX |


class Table(object):
    def __init__(self, groupname):
        self.name = groupname
        #~ self.categories = Category(groupname)
        self.cols = []
        self.data = {}

    def read(self, filename):
        if not os.path.exists(filename):
            raise ValueError('File not existing!')

        cnt = 0
        for l in open(filename).readlines():
            if 'begin{' in l:
                continue
            if 'end{' in l:
                continue
            if '\hline' in l:
                continue
            s = l.replace('\n', '').replace('\\','').split('&')
            if cnt == 0:  # header
                self.cols = []
                for i in xrange(1, len(s)):
                    self.cols.append(s[i].strip())
            else:
                v = {}  # construct data vector for single model
                for i in xrange(1, len(s)):
                    v.update({self.cols[i-1] : s[i].strip()})
                self.data.update({s[0].strip() : v})
            cnt += 1

    def get_data(self, key, col, default='-'):
        if key in self.data.keys():
            r = self.data[key][col]
            r = r.replace('{bf', '{\\bf')
            return r
        else:
            return default


class Writer(object):
    def __init__(self, fname):
        self.fname = fname
        self._open()

    def _open(self):
        if os.path.exists(self.fname):
            os.remove(self.fname)
        self.handler = open(self.fname, 'w')

    def close(self):
        self.handler.close()

    def write(self, s):
        s += '\n'
        self.handler.write(s)


####################
# open report
####################
R = Report('merging_results', 'Merging', 'Alex')
R.open(landscape=False)

####################
# read all tables
####################
# read ALL latex tables
Tdata = []
for g in files:
    k = g.keys()[0]
    T = Table(k)
    T.read(rdir + g[k])
    Tdata.append(T)

# combine tables
# unique key list for all models
mkeys = []
for t in Tdata:
    for k in t.data.keys():
        if k not in mkeys:
            mkeys.append(k)

mkeys.sort()

ofilename = 'merged_radiation.tex'
W = Writer(ofilename)

# estimate number of output columns and wirte HEADERS
ncols = 1
sep = ' & '
header1 = '    \multirow{2}{*}{model} ' + sep
header = '        ' + sep

for T in Tdata:
    ncols += len(T.cols)

c=1
for T in Tdata:
    header1 += '\multicolumn{' + str(len(T.cols)) + '}{|c|}{' + T.name + '}'  # todo multicol here
    if c < len(Tdata):
        header1 += sep
    c += 1
header1 += ' \\\\'

c=1
short_name={}
short_name.update({'MODIS': 'MOD'})
short_name.update({'CLARASAL': 'CLA'})
short_name.update({'GLOBALBEDO-BHR': 'G-BHR'})
short_name.update({'GLOBALBEDO-DHR': 'G-DHR'})
short_name.update({'CERES2.7': 'CER'})
short_name.update({'ISCCP': 'ISCCP'})
short_name.update({'SRBv3.0': 'SRB'})
for T in Tdata:
    for col in T.cols:
        header += short_name[col]  # use only first characters
        if c < ncols-1:
            header += sep
        c += 1
header += ' \\\\'


coltok=''
for i in xrange(ncols):
    coltok += '|c|'

W.write('\\begin{tabular}{l' + coltok + '}')
W.write('    \hline')
W.write(header1)
W.write(header)
W.write('    \hline')

for model in mkeys:
    try:
        s = '    ' + model.split(':')[1] + sep  # one line per model
    except:
        s = '    ' + model + sep  # one line per model
    s = s.replace('-historical', '')
    c = 1
    for T in Tdata:
        for col in T.cols:
            s += T.get_data(model, col)   # get table cell value for model/obs combination
            if c < ncols-1:
                s += sep
            c += 1
    s += ' \\\\'
    W.write(s)
W.write('    \hline')
W.write('\end{tabular}')
W.close()

R.open_table()
R.input(ofilename)
R.close_table(caption='Merged radiation ranking')

R.close()




















