#######################################
# This file is part of pyCMBS
#######################################

# XXX: machine specific paths
# pep8, ignoring some errors like e.g. indention errors or linelength error
PEP = pep8 --ignore=E501
TDIR = ./tmp
VERSION = 1.0.0-dev
TESTDIRS = pycmbs/benchmarking/tests pycmbs/tests


clean :
	find . -name "*.pyc" -exec rm -rf {} \;
	find . -name "y.pkl" -exec rm -rf {} \;
	find . -name "*.so" -exec rm -rf {} \;
	find . -name "data_warnings.log" -exec rm -rf {} \;
	rm -rf C:*debuglog.txt
	rm -rf build
	rm -rf MANIFEST
	rm -rf cover
	rm -rf tmp
	rm -rf doc
	rm -rf dist
	rm -rf pycmbs.egg-info

ship : dist
	# XXX: should not remove non-pycmbs files
	mkdir -p $(TDIR)
	rm -rf $(TDIR)/*
	rm -rfv $(TDIR)/pycmbs.*.tar.gz
	cp ./dist/pycmbs-$(VERSION).tar.gz $(TDIR)
	tar -C $(TDIR) -xvf $(TDIR)/pycmbs-$(VERSION).tar.gz

coverage: dependencies
	nosetests --with-coverage --cover-package=benchmarking --cover-package=pycmbs $(TESTDIRS) --cover-html

tests: dependencies
	nosetests $(TESTDIRS)

dist : clean
	python setup.py sdist

build_docs:
	python setup.py build_sphinx

upload_docs:
	python setup.py upload_sphinx

update_version:
	python autoincrement_version.py

upload_pip: update_version
	# ensure that pip version has always counterpart on github
	git push origin master
	# note that this requres .pypirc file beeing in home directory
	python setup.py sdist upload

dependencies : clean
	sh compile_extensions.sh

pep8 :
	$(PEP) *.py
	$(PEP) ./pycmbs/*.py
	#$(PEP) ./pycmbs/benchmarking/*.py
	$(PEP) ./pycmbs/colormaps/*.py
	$(PEP) ./pycmbs/geostatistic/*.py
	$(PEP) ./pycmbs/diagnostic/*.py
	$(PEP) ./pycmbs/plots/*.py
	$(PEP) ./pycmbs/statistic/*.py



