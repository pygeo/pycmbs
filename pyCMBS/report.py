# -*- coding: utf-8 -*-
"""
This file is part of pyCMBS.
For COPYRIGHT, LICENSE and AUTHORSHIP please referr to
the pyCMBS licensing details.
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

        @param logofile: name of file for logo on first page;
        if None or file not existing, then nothing will be plotted
        @type logofile: str

        @param usehyperlinks: use hyperlinks for table of contents
        @type usehyperlinks: bool

        @param author: Author of the document
        @type author: str

        @param autocompile: ensure automatic PDF creation when
        report is closed
        @type autocompile: bool
        """

        ext = ''
        if filename[:-4] != '.tex':
            ext = '.tex'
        self.filename = outdir + filename+ext
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

#-----------------------------------------------------------------------

    def open(self):
        """
        open report
        """

        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)

        if os.path.exists(self.filename):
            os.remove(self.filename)
        self.file = open(self.filename, 'w')
        self._write_header()

#-----------------------------------------------------------------------

    def close(self):
        """
        close report
        """
        self._write_footer()
        self.file.close()

        if self.autocompile:
            print 'Compiling REPORT ...'
            self.compile()

#-----------------------------------------------------------------------

    def _write_header(self):
        """
        write document header
        """
        self.write('\documentclass{article}')
        self.write('\usepackage{fancyhdr}')
        self.write('\usepackage{graphicx}')
        #facilitates handling of floating environments
        self.write('\usepackage{placeins}')
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
            self._write_single_figure(self.logofile, None)

        self.write('\\newpage')

        self.write('\\tableofcontents')
        self.write('\\newpage')
        self._write_separator()

#-----------------------------------------------------------------------

    def _write_single_figure(self, figpath, caption):
        #/// LOGO ///
        self._write_separator()
        self.write('\\begin{figure}[!htp]')
        self.write('   \centering')
        self.write('   \includegraphics[width=4cm]{' + figpath + '} ')
        if caption is not None:
            self.write('   \caption{' + caption.replace('_', '-') + '}')
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

    def figure(self, f, caption='', width=None, bbox_inches='tight'):
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

        if f is None:
            return

        self.figure_counter += 1
        figname = 'fig_' + str(self.figure_counter).zfill(5) + '.' + self.format

        if width is None:
            width = '12cm'

        self._write_separator()
        self.write('\\begin{figure}[htp]')
        self.write('   \centering')
        self.write('   \includegraphics[width=' + width + ']{'
                   + figname + '} ')
        if len(caption) > 0:
            self.write('   \caption{' + caption.replace('_', '-') + '}')
            self.write('   \label{fig:' + str(self.figure_counter) + '}')
        self.write('\\end{figure}')
        self._write_separator()

        print('Saving figure %s' % self.outdir + figname)
        f.savefig(self.outdir + figname, bbox_inches=bbox_inches, dpi=self.dpi)

#-----------------------------------------------------------------------

    def section(self, s):
        """
        write section header

        @param s: title of section
        @type s: str
        """
        self.clearpage()
        self.write('\section{' + s.replace('_', ' ') + '}')

#-----------------------------------------------------------------------

    def subsection(self, s):
        """
        write subsection header

        @param s: title of subsection
        @type s: str
        """
        #self.write('\clearpage')
        self.barrier()
        self.write('\subsection{' + s.replace('_', ' ') + '}')

#-----------------------------------------------------------------------

    def subsubsection(self, s):
        """
        write subsection header

        @param s: title of subsection
        @type s: str
        """
        self.barrier()
        self.write('\subsubsection{' + s.replace('_', ' ') + '}')

#-----------------------------------------------------------------------

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

#-----------------------------------------------------------------------

    def newpage(self):
        """
        create a new page
        """
        self.clearpage()
        self.write('\\newpage')

#-----------------------------------------------------------------------

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

#-----------------------------------------------------------------------

    def write(self, s):
        """
        write a string to the file

        @param s: string to be written to the file
        @type: str
        """
        self.file.write(s.replace('\f', '\\f').replace('\n', '\\n')
                        .replace('\t', '\\t')
                        .replace('\r', '\\r') + '\n')

#-----------------------------------------------------------------------

    def open_table(self):
        """
        opens a table
        """
        self.write('\\begin[htp]{table}')
        self.write('    \centering')

#-----------------------------------------------------------------------

    def close_table(self, caption='Put a figure caption here'):
        """
        opens a table
        """
        #~ self.write('    \caption{' + caption + '}')
        self.write('\end{table}')

#-----------------------------------------------------------------------

    def input(self, filename):
        """ write an input statement """
        self.write('\\input{' + filename + '}')
        if not os.path.exists(filename):
            print('WARNING: output file used in report not yet existing!')

#-----------------------------------------------------------------------

    def compile(self):
        """
        compile latex document
        """

        curdir = os.getcwd()

        pdfdir = os.path.dirname(self.filename)
        texfile = os.path.basename(self.filename)

        os.chdir(pdfdir)

        #--- compile report twice
        os.system('pdflatex ' + texfile)
        os.system('pdflatex ' + texfile)

        os.chdir(curdir)
