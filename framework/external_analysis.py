'''
Class for running External scripts in pyCMBS
'''
import os

class ExternalAnalysis():
    def __init__(self,executable,template,tags,output_directory='./'):
        '''
        constructor for ExternalAnalysis class

        @param executable: name of system executable to be called to run analysis. This string need to include all options and have a tag <INPUTFILE>
                           to specify the filename to be called e.g. 'matlab -r <INPUTFILE>'
        @type executable: str

        @param template: name of template file to be modified
        @type template: str

        @param tags: tag that specifies which variables should be replaced with which value.
                     The dictionary must have the list consists of dictonary entries like the following example
                     d = [{'tag':'myfile','value':'myvaluestring'},{'tag':'myfile1','value':'myvaluestring1'}]
        @type tags: list

        @param output_directory: output directory to write output of external analysis to
        @type output_directory: str
        '''

        if not '<INPUTFILE>' in executable:
            raise ValueError, 'The commmand needs to include <INPUTFILE> tag'


        self.exe = executable
        self.tags = tags
        self.output_directory = output_directory
        if not os.path.exists(self.output_directory):
            os.makedirs(self.output_directory)
        self.template = template


    def _create_script(self):
        '''
        replace tags in a template file

        returns filename
        '''

        #--- check if template file existing
        if not os.path.exists(self.template):
            raise ValueError, '*** Template file not existing! ' + self.template

        #--- copy template file
        filename = self.output_directory + os.path.basename(self.template) #target file
        if os.path.exists(filename):
            os.remove(filename)

        #--- replace template content
        f = open(self.template,'r')
        s = []
        for l in f.readlines():
            d = l
            for tag in self.tags: #replace all tags
                d = d.replace(tag['tag'],tag['value'])
            s.append(d)
        f.close()

        #--- write output
        of = open(filename,'w')
        of.writelines(s); of.close()

        return filename


    def run(self):
        '''
        run external program
        '''

        #/// generate run script (returns filename) ///
        filename = self._create_script()

        #/// run script
        cmd = self.exe.replace('<INPUTFILE>',filename)
        print cmd
        #~ r = os.system(self.exe) #@todo: use subprocess




#~ template='test.txt'
#~ tags = [{'tag':'ALEX','value':'peakXXXmat'}]
#~
#~ E=ExternalAnalysis('matlab -r <INPUTFILE>',template,tags,output_directory='./tmp/')
#~ E.run()

