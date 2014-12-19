# -*- coding: utf-8 -*-
"""
This file is part of pyCMBS. (c) 2012-2014
For COPYING and LICENSE details, please refer to the file
COPYRIGHT.md
"""

import os
import sys

from matplotlib import pylab as pl


class Report(object):
    """
    A class to generate latex based report
    """

    def __init__(self, filename, title, author, format='png',
                 outdir='./', dpi=300, logofile='Phytonlogo5.pdf',
                 usehyperlinks=True, autocompile=True):
        """
        constructor for Latex report class

        filename : str
            name of output file
        title : str
            name of author
        format : str
            output format for figures (e.g. png, pdf)
        outdir : str
            output directory to write report and images to
        dpi : int
            specify dots per inch for graphic output
        logofile : str
            name of file for logo on first page
            if None or file not existing, then nothing will
            be plotted
        usehyperlinks : bool
            use hyperlinks for table of contents
        author : str
            Author of the document
        autocompile : bool
            ensure automatic PDF creation when
            report is closed
        """

        ext = ''
        if filename[:-4] != '.tex':
            ext = '.tex'
        self.filename = outdir + filename + ext
        self.format = format
        self.title = title
        self.author = author
        self.outdir = outdir
        # needs to be before open()
        self.logofile = logofile
        self.usehyperlinks = usehyperlinks
        self.open()
        self.figure_counter = 0
        self.dpi = dpi
        self.autocompile = autocompile

    def open(self, landscape=False):
        """ open report """
        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)
        if os.path.exists(self.filename):
            os.remove(self.filename)
        self.landscape = landscape
        self.file = open(self.filename, 'w')
        self._write_header()

    def close(self):
        """ close report """
        self._write_footer()
        self.file.close()
        if self.autocompile:
            print 'Compiling REPORT ...'
            self.compile()

    def _write_header(self):
        """ write document header """
        if self.landscape:
            landscape = 'landscape'
        else:
            landscape = ''
        self.write('\documentclass[' + landscape + ']{article}')
        self.write('\usepackage{fancyhdr}')
        self.write('\usepackage{graphicx}')
        self.write('\usepackage{multirow}')
        self.write('\usepackage{multicol}')
        #facilitates handling of floating environments
        self.write('\usepackage{placeins}')
        self.write('\usepackage{tabularx}')
        #self.write('\usepackage{todonotes}')

        self.write('\pagestyle{fancy}')
        self.write('\fancyhf{}')

        # writes section header
        self.write('\lhead{\nouppercase{\leftmark}}')
        self.write('\lfoot{\today}')
        self.write('\rfoot{\thepage}')

        if self.usehyperlinks:
            self.write('\usepackage{hyperref}')
            self.write('\hypersetup{colorlinks,citecolor=black,\
                       filecolor=black,linkcolor=black,urlcolor=black}')

        self.write('\\begin{document}')
        self.write('\\title{' + self.title.replace('_', ' ') + '}')
        self.write('\\author{' + self.author + '}')

        self.write('\maketitle')

        if self.logofile is not None:
            # logo for report
            if os.path.exists(self.logofile):
                self._write_single_figure(self.logofile, None)

        self.write('\\newpage')

        self.write('\\tableofcontents')
        self.write('\\newpage')
        self._write_separator()

    def _write_single_figure(self, figpath, caption):
        #/// LOGO ///
        self._write_separator()
        self.write('\\begin{figure}[!htp]')
        self.write('   \centering')
        self.write('   \includegraphics[width=4cm]{' + figpath + '} ')
        if caption is not None:
            self.write('   \caption{' + caption.replace('_', '-').replace('#', '-') + '}')
        self.write('\\end{figure}')
        self._write_separator()

    def _write_footer(self):
        """ write document footer """
        self._write_separator()
        self.write('\end{document}')

    def _write_separator(self):
        """
        write line with comments (useful to structure document sections)
        """
        self.write('')
        self.write('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
        self.write('')

    def figure(self, f, caption='', width='\\textwidth', height='\\textheight,keepaspectratio', bbox_inches='tight'):
        """
        add a figure string to the report

        Parameters
        ----------
        f : figure
            figure that will be incuded into the report
        caption : str
            caption for the figure to be put in the report
        width : str
            width as string like in latext e.g. width='12cm'
        bbox_inches : str
            option for savefig
        """

        if f is None:
            return

        self.figure_counter += 1
        figname = 'fig_' + str(self.figure_counter).zfill(5) + '.' + self.format
        self._include_figure(figname, caption=caption, width=width, height=height)

        print('Saving figure %s' % self.outdir + figname)
        f.savefig(self.outdir + figname, bbox_inches=bbox_inches, dpi=self.dpi)

    def _include_figure(self, figname, caption='', width='\\textwidth', height='\\textheight,keepaspectratio'):
        """
        include figure in latex file
        """
        self._write_separator()
        self.write('\\begin{figure}[htp]')
        self.write('   \centering')
        self.write('   \includegraphics[width=' + width + ', height=' + height + ']{'
                   + figname + '} ')
        if len(caption) > 0:
            self.write('   \caption{' + caption.replace('_', '-').replace('#', '-') + '}')
            self.write('   \label{fig:' + str(self.figure_counter) + '}')
        self.write('\\end{figure}')
        self._write_separator()

    def section(self, s):
        """
        write section header

        s : str
            title of section
        """
        self.clearpage()
        self.write('\section{' + s.replace('_', ' ') + '}')

    def subsection(self, s):
        """
        write subsection header

        s : str
            title of subsection
        """
        #self.write('\clearpage')
        self.barrier()
        self.write('\subsection{' + s.replace('_', ' ') + '}')

    def subsubsection(self, s):
        """
        write subsection header

        Parameters
        ----------
        s : str
            title of subsection
        """
        self.barrier()
        self.write('\subsubsection{' + s.replace('_', ' ') + '}')

    def capture_figures(self):
        """
        captures all figures that are plotted and
        store them in the report
        """

        print 'Capturing figures and writing to report ...'
        for i in pl.get_fignums():
            f = pl.figure(i)
            self.figure(f)
            self.newpage()

    def newpage(self):
        """ create a new page """
        self.clearpage()
        self.write('\\newpage')

    def clearpage(self):
        """
        create a new page
        as an alternative, one could also use the placeins package

        http://tug.ctan.org/tex-archive/info/l2picfaq/german/l2picfaq.pdf

        \usepackage{placeins}

        ...

        \FloatBarrier

        This ensures that all Figures/tables before the breakpoint
        are put to paper
        WITHOUT generating a new page. It is thus the opposite to
        \clearpage

        """
        self.write('\clearpage')

    def barrier(self):
        self.write('\FloatBarrier')

    def write(self, s):
        """
        write a string to the file

        Parameters
        ----------
        s : str
            string to be written to the file
        """
        self.file.write(s.replace('\f', '\\f').replace('\n', '\\n')
                        .replace('\t', '\\t')
                        .replace('\r', '\\r') + '\n')

    def open_table(self):
        """ opens a table """
        self.write('\\begin{table}[htp]')
        self.write('    \centering')

    def close_table(self, caption='Put a figure caption here'):
        """ closes a table """
        self.write('    \caption{' + caption.replace('_', '-').replace('#', '-') + '}')
        self.write('\end{table}')

    def input(self, filename):
        """ write an input statement """
        self.write('\\input{' + filename + '}')
        if not os.path.exists(filename):
            print('WARNING: output file used in report not yet existing!')

    def compile(self):
        """
        compile latex document
        """
        curdir = os.getcwd()

        pdfdir = os.path.dirname(self.filename)
        texfile = os.path.basename(self.filename)
        os.chdir(pdfdir)

        # compile report twice
        os.system('pdflatex ' + texfile)
        os.system('pdflatex ' + texfile)
        os.chdir(curdir)
