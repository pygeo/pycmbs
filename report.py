#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Alexander Loew"
__version__ = "0.1"
__date__ = "2012/10/29"
__email__ = "alexander.loew@zmaw.de"

'''
# Copyright (C) 2012 Alexander Loew, alexander.loew@zmaw.de
# See COPYING file for copying and redistribution conditions.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; version 2 of the License.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FI/TNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
'''

'''
module to generate a LaTex report
'''


import os,sys

from matplotlib import pylab as pl

class Report():
    """
    A class to generate latex based report

    @todo: example how to use report class
    """

    def __init__(self,filename,title,author,format='png',outdir='./',dpi=300,logofile='Phytonlogo5.pdf'):
        """
        constructor for Latex report class

        @param filename: name of output file
        @type filename: str

        @param title: name of author
        @type title: str

        @param format: output format for figures (e.g. png, pdf)
        @type format: str

        @param outdir: output directory to write report and images to
        @type outdir: str

        @param dpi: specify dots per inch for graphic output
        @type dpi: int

        @param logofile: name of file for logo on first page; if None or file not existing, then nothing will be plotted
        @type logofile: str
        """

        ext = ''
        if filename[:-4] != '.tex':
            ext='.tex'
        self.filename=outdir + filename+ext
        self.format = format
        self.title=title
        self.author=author
        self.outdir=outdir
        self.logofile = logofile #needs to be before open()
        self.open()
        self.figure_counter = 0
        self.dpi = dpi

#-----------------------------------------------------------------------

    def open(self):
        """
        open report
        """

        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)

        if os.path.exists(self.filename):
            os.remove(self.filename)
        self.file = open(self.filename,'w')
        self._write_header()






#-----------------------------------------------------------------------

    def close(self):
        """
        close report
        """
        self._write_footer()
        self.file.close()

#-----------------------------------------------------------------------

    def _write_header(self):
        """
        write document header
        """
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

        #if self.logofile != None:
        #    self._write_single_figure(self.logofile,None) #logo for report


        self.write('\\newpage')

        self.write('\\tableofcontents')
        self.write('\\newpage')
        self._write_separator()

#-----------------------------------------------------------------------

    def _write_single_figure(self,figpath,caption):
        #/// LOGO ///
        self._write_separator()
        self.write('\\begin{figure}[htp]')
        self.write('   \centering')
        self.write('   \includegraphics[width=4cm]{'+ figpath + '} ')
        if caption != None:
            self.write('   \caption{' + caption.replace('_','-') + '}')
        self.write('\\end{figure}')
        self._write_separator()

#-----------------------------------------------------------------------

    def _write_footer(self):
        """
        write document footer
        """
        self._write_separator()
        self.write('\end{document}')

#-----------------------------------------------------------------------

    def _write_separator(self):
        """
        write line with comments (useful to structure document sections)
        """
        self.write('')
        self.write('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
        self.write('')

#-----------------------------------------------------------------------

    def figure(self,f,caption='',width=None,bbox_inches='tight'):
        """
        add a figure string to the report

        @param f: figure that will be incuded into the report
        @type f: matplotlib figure object

        @param caption: caption for the figure to be put in the report
        @type caption: str

        @param width: width as string like in latext e.g. width='12cm'
        @type width: str

        @param bbox_inches: option for savefig
        @type bbox_inches: str
        """

        if f == None:
            return

        self.figure_counter +=1
        figname = 'fig_' + str(self.figure_counter).zfill(5) + '.' + self.format

        if width == None:
            width = '12cm'


        self._write_separator()
        self.write('\\begin{figure}[htp]')
        self.write('   \centering')
        self.write('   \includegraphics[width=' + width + ']{' + figname + '} ')
        if len(caption)>0:
            self.write('   \caption{' + caption.replace('_','-') + '}')
            self.write('   \label{fig:' + str(self.figure_counter) + '}')
        self.write('\\end{figure}')
        self._write_separator()

        f.savefig(self.outdir + figname, bbox_inches=bbox_inches,dpi=self.dpi)

#-----------------------------------------------------------------------

    def section(self,s):
        """
        write section header

        @param s: title of section
        @type s: str
        """
        self.write('\clearpage')
        self.write('\section{' + s.replace('_',' ') + '}')

#-----------------------------------------------------------------------

    def subsection(self,s):
        """
        write subsection header

        @param s: title of subsection
        @type s: str
        """
        self.write('\clearpage')
        self.write('\subsection{' + s.replace('_',' ') + '}')

#-----------------------------------------------------------------------

    def subsubsection(self,s):
        """
        write subsection header

        @param s: title of subsection
        @type s: str
        """
        self.write('\clearpage')
        self.write('\subsubsection{' + s.replace('_',' ') + '}')

#-----------------------------------------------------------------------

    def capture_figures(self):
        """
        captures all figures that are plotted and
        store them in the report
        """

        print 'Capturing figures and writing to report ...'
        for i in pl.get_fignums():
            f=pl.figure(i)
            self.figure(f)
            self.newpage()

#-----------------------------------------------------------------------

    def newpage(self):
        """
        create a new page
        """
        self.write('\clearpage')
        self.write('\\newpage')

#-----------------------------------------------------------------------

    def clearpage(self):
        """
        create a new page
        """
        self.write('\clearpage')

#-----------------------------------------------------------------------

    def write(self,s):
        """
        write a string to the file

        @param s: string to be written to the file
        @type: str
        """
        self.file.write(s.replace('\f','\\f').replace('\n','\\n').replace('\t','\\t').replace('\r','\\r') + '\n')

#-----------------------------------------------------------------------

    def compile(self):
        """
        compile latex document
        """
        os.system('pdflatex ' + self.filename)
