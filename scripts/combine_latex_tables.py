"""
This file is part of pyCMBS. (c) 2012-2014
For COPYING and LICENSE details, please refer to the file
COPYRIGHT.md
"""

"""
Combine latex tables that were written using gleckler.write_table routine
"""

import os

# specify directory name
rdir = '/home/m300028/shared/temp/ensmean/reports/'

files = ['ranking_table_albedo.tex', 'ranking_table_sis.tex', 'ranking_table_surface_upward_flux.tex']
group_labels = ['ALBEDO', 'SIS', 'Rup']

if len(files) != len(group_labels):
    raise ValueError('Size of files and labels needs to be similar.')


#|        Albedo       |     RN      |
#| CLARA | MODIS | ... | CLARA | XXX |


class Table(object):
    def __init__(self, groupname):
        self.group = Group(groupname)

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
                for x in s:
                    self.group._add_column(x.lstrip())
            else:
                v = {}  # construct data vector
                print s
                for i in xrange(len(h)):
                    v.update({h[i] : s[i].strip()})


                r.update({s[0].strip() : v})
            cnt += 1

        return r




class Group(object):
    def __init__(self, k):
        self.columns = {}
        self.name = k

    def _add_column(self, k):
        """ add an item e.g. an observation """

    store also sequence !!!

        if k not in self.columns.keys():
            self.columns.update({k: {}})

    def _add_value(self, k, id, v):
        if k not in self.columns.keys():
            raise ValueError('This key is not existing: %s' % k)
        if id in self.columns[k].keys():
            raise ValueError('An entry is already existing: %s' % id)
        else:
            self.columns[k].update({id: v})



T = Table('Albedo')
T.read(rdir+files[0])




#~ def parse_file(filename):
    #~ if not os.path.exists(filename):
        #~ raise ValueError('File not existing!')
    #~ cnt = 0
    #~ for l in open(filename).readlines():
        #~ if 'begin{' in l:
            #~ continue
        #~ if 'end{' in l:
            #~ continue
        #~ if '\hline' in l:
            #~ continue
        #~ s = l.replace('\n', '').replace('\\','').split('&')
        #~ if cnt == 0:  # header
            #~ r = {}
            #~ h = []
            #~ for x in s:
                #~ h.append(x.lstrip())
        #~ else:
            #~ v = {}  # construct data vector
            #~ for i in xrange(len(h)):
                #~ v.update({h[i] : s[i].strip()})
#~
#~
            #~ r.update({s[0].strip() : v})
        #~ cnt += 1
#~
    #~ return r


# read latex tables
T = []
for file in files:
    T.append(parse_file(rdir + file))

# combine tables

# unique key list for all models
keys = []
for t in T:
    for k in t.keys():
        if k not in keys:
            keys.append(k)

# create table for each key

for k in keys:
    for t in T:
        x = t[k]  # dict with obs
        print k, x

# write results per group
for i in xrange(len(group_labels)):
    print 'Group: ', group_labels[i]

    #~ obs_keys =
    #~ for o_k in

    obs_keys = T[i]

    print 'obs_keys: ', obs_keys

    #for model in keys():
    #    print model









