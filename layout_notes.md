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
        pyCMBS──│
                │ __init__.py   (new!)
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


----------------
Comments from MI:

> I  dont like too much the `benchmark_models.py`. We should try to find a better name for this. What  about `pycmbs-benchmarking.py` or `pycmbs-score.py`
`pycmbs-benchmarking.py` sounds good (`benchmark_models.py` was just a stub)

> What is the purpose of the `Makefile` on the same level like `core` and `benchmarking`?
I think there was a bit of confusion. There must be one more directory level that defines the lib. Yes, there should be only one Makefile, in the project root. I will paste an example below.

>  Having the `benchmark_models.py` on the uppermost level would be o.k. with me. However I dont see a strong need here, as the installation procedure will basically ensure that it will be found in the system path anyhow. Thus for the user it doesnt matter. Its more that it might matter from a developer perspective.

Yes, that should make the life of the developers much easier. Basically the aim is to avoid setting PYTHONPATH for development branches, as it causes collisions with production code. Placed in the root level the application will have exactly the same import statements as if it would be outside of the project, and it will use modules from the current (development) branch. 

The cool thing is that `nosetests` will automatically track where is the root of the project and collect the full path to the modules. It allows parts of the library to communicate just as they would normally do when invoked by user. Also tests will be able to use appropriate import statements, practically eliminating the need to include the project under development in the system path. Here is what I mean:
 


    pyCMBS-v1.00
    ├── configuration (not part of the library, but part of the package)
    │   ├── models.json
    │   ├── analysis.json
    │   └── parameter.ini
    ├── docs
    ├── docsrs
    ├── pycmbs (module and library names shall be always in lowercase)
    │   ├── benchmarking
    │   │   ├── analysis.py
    │   │   ├── __init__.py
    │   │   ├── test_analysis.py
    │   ├── core
    │   │   ├── data.py
    │   │   ├── grid.py
    │   │   ├── __init__.py
    │   │   └── test_data.py
    │   └── __init__.py (this should be the first/deepest occurence of __init__.py)
    ├── INSTALL
    ├── Makefile
    ├── MANIFEST
    ├── pycmbs-benchmarking.py
    ├── README
    └── setup.py

If `grid.py` or `test_data.py`  need to import `data.py`, with nosetests a statement like `from pycmbs.core import data` will just work. it wont work otherwise unless the PYTHONPATH is set somewhere. 

 ** Layout **: if users are going to frequently use `pycmbs.core` in their scripts or for interactive work its probably better to make the path shorter and skip `core` completely, unless you see that it may cause name collisions in the future. Now it looks pretty safe. Just to improve users workflow I would probably use such layout:


    pyCMBS-v1.00
    ├── configuration
    │   ├── analysis.json
    │   ├── models.json
    │   └── parameter.ini
    ├── docs
    ├── docsrs
    ├── pycmbs
    │   ├── benchmarking
    │   │   ├── analysis.py
    │   │   ├── __init__.py
    │   │   ├── models.py
    │   │   ├── test_analysis.py
    │   │   └── test_models.py
    │   ├── data.py
    │   ├── grid.py
    │   ├── __init__.py
    │   ├── test_data.py
    │   └── test_grid.py
    ├── Makefile
    ├── pycmbs-benchmarking.py
    └── setup.py

** Benchmarking **: If `benchmarking` wont become totally independent of the `pycmbs` library, it may be easier just to keep it inside the `pycmbs` module as displayed above. 

Actually my suggestion is to change as little as absolutely necessary and to make pycmbs easy to maintain and to use. It seems to me that many parts are already isolated quite well, I would mostly like to move the application and configuration parts out of the library space and make the package more compliant to the standards (thus making it more attractive to the others).
