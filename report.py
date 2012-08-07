#!/usr/bin/python
# -*- coding: utf-8 -*-

__author__ = "Alexander Loew"
__version__ = "0.0"
__date__ = "0000/00/00"

'''
module to generate a LaTex report
'''


import os,sys

from matplotlib import pylab as pl

class Report():
    '''
    A class to generate latex based report

    @todo: example how to use report class
    '''

    def __init__(self,filename,title,author,format='png',outdir='./'):
        '''
        constructor for Latex report class

        @param filename: name of output file
        @type filename: str

        @param title: name of author
        @type title: str

        @param format: output format for figures (e.g. png, pdf)
        @type format: str

        @param outdir: output directory to write report and images to
        @type outdir: str
        '''
        ext = ''
        if filename[:-4] != '.tex':
            ext='.tex'
        self.filename=outdir + filename+ext
        self.format = format
        self.title=title
        self.author=author
        self.outdir=outdir
        self.open()
        self.figure_counter = 0

#-----------------------------------------------------------------------

    def open(self):
        '''
        open report
        '''
        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)
        if os.path.exists(self.filename):
            os.remove(self.filename)
        self.file = open(self.filename,'w')
        self._write_header()

#-----------------------------------------------------------------------

    def close(self):
        '''
        close report
        '''
        self._write_footer()
        self.file.close()

#-----------------------------------------------------------------------

    def _write_header(self):
        '''
        write document header
        '''
        self.write('\documentclass{article}')
        self.write('\usepackage{fancyhdr}')
        self.write('\usepackage{graphicx}')
        #self.write('\usepackage{todonotes}')

        self.write('\pagestyle{fancy}')
        self.write('\fancyhf{}')

        self.write('\lhead{\nouppercase{\leftmark}}')  #%writes section header
        self.write('\lfoot{\today}')
        self.write('\rfoot{\thepage}')

        self.write('\\begin{document}')
        self.write('\\title{' + self.title.replace('_',' ') + '}')
        self.write('\\author{' + self.author +'}')

        self.write('\maketitle')
        self.write('\\newpage')
        self._write_separator()

#-----------------------------------------------------------------------

    def _write_footer(self):
        '''
        write document footer
        '''
        self._write_separator()
        self.write('\end{document}')

#-----------------------------------------------------------------------

    def _write_separator(self):
        '''
        write line with comments (useful to structure document sections)
        '''
        self.write('')
        self.write('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
        self.write('')

#-----------------------------------------------------------------------

    def figure(self,f,caption=''):
        '''
        add a figure string to the report

        @param f: figure that will be incuded into the report
        @type f: matplotlib figure object
        '''

        self.figure_counter +=1
        figname = 'fig_' + str(self.figure_counter).zfill(5) + '.' + self.format

        self._write_separator()
        self.write('\\begin{figure}[htp]')
        self.write('   \centering')
        #self.write('   \includegraphics[width=12cm]{' + figname + '} \\\ ')
        self.write('   \includegraphics[width=12cm]{' + figname + '} ')
        self.write('   \caption{' + caption + '}')
        self.write('   \label{fig:' + str(self.figure_counter) + '}')
        self.write('\\end{figure}')
        self._write_separator()

        f.savefig(self.outdir + figname, bbox_inches='tight')

    def section(self,s):
        '''
        write section header

        @param s: title of section
        @type s: str
        '''
        self.write('\clearpage')
        self.write('\section{' + s.replace('_',' ') + '}')

#-----------------------------------------------------------------------

    def capture_figures(self):
        '''
        captures all figures that are plotted and
        store them in the report
        '''

        print 'Capturing figures and writing to report ...'
        for i in pl.get_fignums():
            f=pl.figure(i)
            self.figure(f)
            self.newpage()

#-----------------------------------------------------------------------

    def newpage(self):
        '''
        create a new page
        '''
        self.write('\clearpage')
        self.write('\\newpage')

#-----------------------------------------------------------------------

    def write(self,s):
        '''
        write a string to the file

        @param s: string to be written to the file
        @type: str
        '''
        self.file.write(s.replace('\f','\\f').replace('\n','\\n').replace('\t','\\t').replace('\r','\\r') + '\n')

#-----------------------------------------------------------------------

    def compile(self):
        '''
        compile latex document
        '''
        os.system('pdflatex ' + self.filename)
