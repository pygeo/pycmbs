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
 * `MI`: There can be a dedicated folder just for integration/acceptance tests. It also makes running unittests much faster as they are separated from the acceptance tests which normally are computationally more intensive.
 * `MI`: Python applications, like `pycmbs.py` could be just placed in the root space. Projects like [`virtualenv`](https://github.com/pypa/virtualenv/blob/develop/virtualenv.py) (with standalone application) do it that way.

# Comments (AL)

I think I am in favor of the side-by-side layout as it is cleaner. Thanks for the idea. Some remarks below

 * What is purpose of having `docs` and `apidocs` ?
 * I suggest that we replace the `pycmbs` subfolder by `core`
 * The structure of the *whole* package would look then e.g. like

        MANIFEST.in
        Makefile
        setup.py
        pyCMBS──|
                | __init__.py   (new!)
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
                ├── core
                │   ├── data.py
                │   ├── __init__.py
                │   ├── test_data.py
                │   ├── test_utils.py
                │   └── utils.py
                ├── apidocs
                ├── benchmark_models.py (currently "pycmbs.py")
                └── Makefile (purpose of this??)

The installation would then by like `python setup.py install` which would copy the `pyCMBS` directory in `dist-packages`. Import statements in any script would then look like.
`from pyCMBS.core import Data` or `from pyCMBS.benchmarking import Model`. I am still not very familar with best way of python module structuring, thus there might be a cleaner way to do it.

 * I don't like too much the `benchmark_models.py`. We should try to find a better name for this. What  about `pycmbs-benchmarking.py` or `pycmbs-score.py`
 * What is the purpose of the `Makefile` on the same level like `core` and `benchmarking`?
 * Having the `benchmark_models.py` on the uppermost level would be o.k. with me. However I don't see a strong need here, as the installation procedure will basically ensure that it will be found in the system path anyhow. Thus for the user it doesn't matter. It's more that it might matter from a developer perspective.








