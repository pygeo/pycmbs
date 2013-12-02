#######################################
# This file is part of pyCMBS
#######################################


clean :
	find . -name "*.pyc" -exec rm -rf {} \;
	rm -rf build

dist : clean
	python setup.py sdist
	rm MANIFEST




