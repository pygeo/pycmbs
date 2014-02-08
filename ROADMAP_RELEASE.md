Roadmap for pycmbs release on github
====================================

github move
-----------

* :-) create github organization account [AL]
 * :-) Done: https://github.com/pygeo/pycmbs.git
* familiarize with github tools for code deployment [MI]
* Jenkins integration? check whihc options for automatic code testing,
  coverage testing and testing for different environments are possible [MI]
* move pycmbs from assembla to github [AL,MI]


Code clean up and completion
----------------------------

### Prio 1

* unittests for large parts of the code (at least as dummy tests) [MI]
* objective of tests: avoid that there are still too much not working code
  parts because of missing import statements (not everything covered yet by tests)

### Prio 2

* PEP8 compliance of code [AL, MI]
* consistency of docstrings [AL, MI]


Documentation
-------------

### Prio 1

* review current sphinx documentation and put placeholders where still needed
* [AL]; MI to review
* expert installation descritpion for developers [AL]; MI to review
* review general README [AL], MI to review
* review LICENSE [AL]

### Prio 2

* improved installation using setup.py for general users
