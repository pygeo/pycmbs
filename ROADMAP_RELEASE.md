Roadmap for pycmbs release on github
====================================

github move
-----------

* [AL] revise documentation on REDMINE Wiki


### Prio 2

* deployment using pip (information also using testing server can be found
  here: https://pypi.python.org/pypi
* pip installation covering also external dependencies recursevely
* adapt ZMAW WIKI

Code clean up and completion
----------------------------

### Prio 1

* [MI] ensure that all files have conssitent file header with link to copyright
* file
* Remove machine/user specific information (.cfg, .ini, etc), use/develop sample config generator
* [MI] unittests for benchmarking
* put example data as TARBALL somewhere for download
 * data required for testing
 * data required for examples (note that this can also be automatically
   downloaded)


Documentation
-------------

### Prio 1

* [AL] fill file INSTALL.md with content
* [AL] review current sphinx documentation and put placeholders where still needed
  ; MI to review
* [MI] test again the examples scripts in /pycmbs/examples; would be good if Mikhail
  would do that as he was not included in its development
* [MI] review general README, which has been already adapted by AL
* [MI] put a copyright notice in every file. Perhaps git could do this
  automatically?? [MI to check]

### Prio 2

* improved installation using setup.py for general users
