Roadmap for pycmbs release on github
====================================

github move
-----------

* [MI] familiarize with github tools for code deployment
* [MI] Jenkins integration? check whihc options for automatic code testing,
  coverage testing and testing for different environments are possible
  * ** MI: ** Jenkins require our own running server. For publically available github projects
    there is [coveralls.io][http://coveralls.io] for unittests and coverage, 
    and [travis-ci.org](http://travis-ci.org) service builds the software and uploads it to pip.
  
* decide for a version number: v1.0 ???
  * **MI**: Semantic versioning would be probably best. 
    It uses the following version format: MAJOR.MINOR.PATCH. Check it out: [semver.org](http://semver.org) 
    I would assign it version 1.0.0, there is nothing wrong giving it a solid number, 
    and it tells something about the current API version.

### Prio 2

* deployment using pip (information also using testing server can be found
  here: https://pypi.python.org/pypi
* generation of tagged versions?
* adapt ZMAW WIKI

Code clean up and completion
----------------------------

### Prio 1

* Remove machine/user specific information (.cfg, .ini, etc), use/develop sample config generator
* [MI] unittests for large parts of the code (at least as dummy tests)
* objective of tests: avoid that there are still too much not working code
  parts because of missing import statements (not everything covered yet by tests)
* put example data as TARBALL somewhere for download
 * data required for testing
 * data required for examples (note that this can also be automatically
   downloaded)

### Prio 2

* [AL,MI] consistency of docstrings


Documentation
-------------

### Prio 1

* [AL] review current sphinx documentation and put placeholders where still needed
  ; MI to review
* [MI] test again the examples scripts in /pycmbs/examples; would be good if Mikhail
  would do that as he was not included in its development
* [AL] expert installation descritpion for developers; MI to review
* test installation in clean environment (virtualenv??)
* [AL] review general README, MI to review
* [MI] put a copyright notice in every file. Perhaps git could do this
  automatically?? [MI to check]

### Prio 2

* improved installation using setup.py for general users
