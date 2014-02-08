#######################################
# This file is part of pyCMBS
#######################################

# XXX: machine specific paths
# pep8, ignoring some errors like e.g. indention errors or linelength error
PEP = pep8 --ignore=E501,E128,E111,E127
TDIR = /tmp
VERSION = 0.1.6

clean :
	find . -name "*.pyc" -exec rm -rf {} \;
	find . -name "data_warnings.log" -exec rm -rf {} \;
	rm -rf build
	rm -rf MANIFEST
	rm -rf cover

ship : dist
	rm -rf $(TDIR)/*
	cp ./dist/pyCMBS-$(VERSION).tar.gz $(TDIR)
	tar -C $(TDIR) -xvf $(TDIR)/pyCMBS-$(VERSION).tar.gz

coverage:
	nosetests --with-coverage --cover-package=benchmarking --cover-package=pycmbs pycmbs/benchmarking/tests pycmbs/tests --cover-html

dist : clean
	python setup.py sdist

pep8 :
	$(PEP) *.py
	$(PEP) ./pycmbs/*.py
	$(PEP) ./pycmbs/benchmarking/*.py



