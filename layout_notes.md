# Layouts
## Nested layout
 * `benchmarking` is always a part of `pycmbs`
 * Imports look like `from pycmbs.benchmarking import model`
 * Unittests reside in the same directory for easy packaging

        .
        ├── configuration
        │   ├── analysis.json
        │   ├── cfc.ini
        │   ├── dummy_project.cfg
        │   └── models.json
        ├── docs
        ├── pycmbs
        │   ├── benchmarking (currently "framework")
        │   │   ├── config.py
        │   │   ├── __init__.py
        │   │   ├── model.py
        │   │   ├── test_model.py
        │   │   └── test_config.py
        │   ├── data.py
        │   ├── __init__.py
        │   ├── test_data.py
        │   ├── test_utils.py
        │   └── utils.py
        ├── apidocs
        ├── benchmark_models.py (currently "pycmbs.py")
        └── Makefile

## Side by side layout
 * `benchmarking` module resides in the same root as pycmbs.
    * That way its easy to split the packages apart at any point.
 * Imports look like: `from benchmarking import model`
 * Unittests reside in the same directory for easy packaging

        .
        ├── benchmarking (currently "framework")
        │   ├── config.py
        │   ├── __init__.py
        │   ├── test_model.py
        │   ├── model.py
        │   └── test_config.py
        ├── configuration
        │   ├── analysis.json
        │   ├── cfc.ini
        │   ├── dummy_project.cfg
        │   └── models.json
        ├── docs
        ├── pycmbs
        │   ├── data.py
        │   ├── __init__.py
        │   ├── test_data.py
        │   ├── test_utils.py
        │   └── utils.py
        ├── apidocs
        ├── benchmark_models.py (currently "pycmbs.py")
        └── Makefile

# Notes
 * `MI`: according to PEP8 module names should have short lowercase names. Name `pyCMBS` could be the name of the package, but importable module shall be in lowercase: `pycmbs`
 * `MI`: There can be a dedicated folder just for integration/acceptance tests. It also makes running unittests much faster.
 * `MI`: Python applications, like `pycmbs.py` could be just placed in the root space. Projects like [`virtualenv`](https://github.com/pypa/virtualenv/blob/develop/virtualenv.py) (with standalone application) do it that way.
