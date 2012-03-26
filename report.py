#!/usr/bin/python
# -*- coding: utf-8 -*-

__author__ = "Alexander Loew"
__version__ = "0.0"
__date__ = "0000/00/00"

import os,sys

from matplotlib import pylab as pl

class Report():
    '''
    A class to generate latex based report
    '''
    
    def __init__(self,filename,title,author,format='png',outdir='./'):
        ext = ''
        if filename[:-4] != '.tex':
            ext='.tex'
        self.filename=outdir + filename+ext
        self.format = format
        self.title=title
        self.author=author
        self.outdir=outdir
        self.open()
        
    def open(self):
        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)
        if os.path.exists(self.filename):
            os.remove(self.filename)
        self.file = open(self.filename,'w')
        self._write_header()
        
    def close(self):
        self._write_footer()
        self.file.close()
        
    def _write_header(self):
        self.write('\documentclass{article}')
        self.write('\usepackage{graphicx}')
        #self.write('\usepackage{todonotes}')

        self.write('\\begin{document}')
        self.write('\\title{' + self.title.replace('_',' ') + '}')
        self.write('\\author{' + self.author +'}')

        self.write('\maketitle')
        self.write('\\newpage')
        self._write_separator()
        
        
    def _write_footer(self):
        self._write_separator()
        self.write('\end{document}')
        
    def _write_separator(self):
        self.write('')
        self.write('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
        self.write('')
        
    def figure(self,f,caption=''):
        '''
        add a figure string to the report
        
        f: matplotlib figure object
        '''
        
        figname = 'fig_' + str(f.number).zfill(4) + '.' + self.format
        
        self._write_separator()
        self.write('\\begin{figure}[htp]')
        self.write('   \centering')
        self.write('   \includegraphics[width=12cm]{' + figname + '} \\\ ')
        self.write('   \caption{' + caption + '}')
        self.write('   \label{fig:' + str(f.number) + '}')
        self.write('\\end{figure}')
        self._write_separator()
        
        f.savefig(self.outdir + figname, bbox_inches='tight')
        
    def section(self,s):
        self.write('\section{' + s + '}')
        
        
    def capture_figures(self):
        '''
        captures all figures that are plotted and
        stores them in the report
        '''
        
        print 'Capturing figures and writing to report ...'
        for i in pl.get_fignums():
            f=pl.figure(i)
            self.figure(f)
            self.newpage()
    
    def newpage(self):
        self.write('\clearpage')
        self.write('\\newpage')
        
    def write(self,s):
        '''
        write a string to the file
        '''
        self.file.write(s + '\n')
        
    def compile(self):
        os.system('pdflatex ' + self.filename)
