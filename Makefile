#######################################
# This file is part of pyCMBS
#######################################

# XXX: machine specific paths
# pep8, ignoring some errors like e.g. indention errors or linelength error
PEP = pep8
TDIR = ./tmp
VERSION = 0.1.6

clean :
	find . -name "*.pyc" -exec rm -rf {} \;
	find . -name "data_warnings.log" -exec rm -rf {} \;
	rm -rf build
	rm -rf MANIFEST
	rm -rf cover
	rm -rf tmp
	rm -rf docs
	rm -rf dist

ship : dist
	# XXX: should not remove non-pycmbs files
	rm -rf $(TDIR)/*
	rm -rfv $(TDIR)/pycmbs.*.tar.gz
	cp ./dist/pycmbs-$(VERSION).tar.gz $(TDIR)
	tar -C $(TDIR) -xvf $(TDIR)/pycmbs-$(VERSION).tar.gz

coverage:
	nosetests --with-coverage --cover-package=benchmarking --cover-package=pycmbs pycmbs/benchmarking/tests pycmbs/tests --cover-html

dist : clean
	python setup.py sdist

pep8 :
	$(PEP) *.py
	$(PEP) ./pycmbs/*.py
	$(PEP) ./pycmbs/benchmarking/*.py



