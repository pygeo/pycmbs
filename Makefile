#######################################
# This file is part of pyCMBS
#######################################

clean :
	find . -name "*.pyc" -exec rm -rf {} \;
	rm -rf build
	rm -rf MANIFEST

dist : clean
	python setup.py sdist





