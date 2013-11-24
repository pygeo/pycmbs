from distutils.core import setup
import os

import pyCMBS

install_requires=["numpy>0.1","cdo>1.2",]


setup(name='pyCMBS',
      version=pyCMBS.__version__,
      py_modules=['pyCMBS',],
      packages=['pyCMBS','pyCMBS/framework','pyCMBS/tests','pyCMBS/logo','pyCMBS/misc','pyCMBS/examples','pyCMBS/tools'],
      author = "Alexander Loew",
      author_email='alexander.loew@mpimet.mpg.de',
      maintainer='Alexander Loew',
      maintainer_email='alexander.loew@mpimet.mpg.de',
      url='https://code.zmaw.de/projects/pycmbs',
      description='pyCMBS - python Climate Model Benchmarking Suite',
      long_description='The pyCMBS project is a suite of tools to process, analyze, visualize and benchmark scientific model output against each other or against observational data. It is in particular useful for analyzing in an efficient way output from climate model simulations.',
      install_requires=install_requires,
      keywords = ["data","science","climate","meteorology","model evaluation","benchmarking","metrics"],
      scripts = ["pyCMBS/framework/pycmbs.py"],
      license = "XXXX"      )

#print('Compressing sample data ...')
#ofile='./dist/example_data.tgz'
#cmd='tar -cvzf ' + ofile + ' example_data'
#if os.path.exists(ofile):
#    os.remove(ofile)
#os.system(cmd)





#PIP
#~ python setup.py register
#~ python setup.py sdist
#~ python setup.py upload





