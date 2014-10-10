# -*- coding: utf-8 -*-
"""
This file is part of pyCMBS. (c) 2012-2014
For COPYING and LICENSE details, please refer to the file
COPYRIGHT.md
"""

import os


class ExternalAnalysis():
    def __init__(self, executable, template, tags, output_directory='./', options=''):
        """
        external analysis to run arbitrary scripts using the command line

        Parameters
        ----------
        executable : str
            name of system executable to be called to run analysis.
            This string need to include all options and have a tag <INPUTFILE>
            to specify the filename to be called e.g. 'matlab -r <INPUTFILE>'
        template : str
            name of template file to be modified
        tags : dict
            dictionary that specifies which variables should be replaced with
            which value. e.g.
            d = {'tag1':'value1','tag2':'value2' ...}
        output_directory : str
            output directory to write output of external analysis to
        options : str
            options that will be appended to the filename in the end
        """

        if not '<INPUTFILE>' in executable:
            raise ValueError('The commmand needs to include <INPUTFILE> tag')
        if not isinstance(tags, dict):
            raise ValueError('Tags need to be provided in dictionary!')

        self.exe = executable
        self.tags = tags
        self.output_directory = output_directory
        if not os.path.exists(self.output_directory):
            os.makedirs(self.output_directory)
        if self.output_directory[-1] != os.sep:
            self.output_directory += '/'
        self.template = template
        if not os.path.exists(self.template):
            print self.template
            raise ValueError('Template file not existing!')
        self.options = options

        self._check()

    def _check(self):
        for k in self.tags.keys():
            if not isinstance(self.tags[k], str):
                print self.tags[k]
                raise ValueError('ERROR: all tags are required to be strings!')

    def _create_script(self):
        """
        replace tags in a template file

        Returns
        -------
        returns filename of file to process
        """
        # check if template file existing
        if not os.path.exists(self.template):
            raise ValueError('ERROR: Template file not existing! %s' % self.template)

        # copy template file
        filename = self.output_directory + os.path.basename(self.template)
        if os.path.exists(filename):
            os.remove(filename)

        # replace template content
        f = open(self.template, 'r')
        s = []
        for l in f.readlines():
            d = l
            for k in self.tags.keys():  # replace all tags
                d = d.replace('<' + k + '>', self.tags[k])
            s.append(d)
        f.close()

        # write output
        of = open(filename, 'w')
        of.writelines(s)
        of.close()

        return filename

    def run(self, execute=True, remove_extension=False):
        """
        run external program

        the execution of the program will be done in the directory where
        the modified script is located. this is done, because programs
        like e.g. matlab have problems when the scriptname is provided with
        a pathname

        Parameters
        ----------
        execute : bool
            execute command = run in shell
            if false, then the command is simply printed
            (good for debugging)
        remove_extension : bool
            remove extension from filename.
            This is for instance needed for matlab calls from the
            command line, as matlab does NOT accept script names with
            a fileextension
        """

        #/// generate run script (returns filename) ///
        filename = self._create_script()
        if remove_extension:
            filename = os.path.splitext(filename)[0]

        #/// split filename into path and
        thedir = os.path.dirname(filename)
        thefile = os.path.basename(filename)
        curdir = os.getcwd()

        #/// run script
        cmd = self.exe.replace('<INPUTFILE>', thefile) + self.options

        rstatus = False
        if execute:
            os.chdir(thedir)  # change to directory where script is located
            r = os.system(cmd)  # run it
            os.chdir(curdir)  # go back
        else:
            print cmd
