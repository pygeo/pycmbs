# CHANGES

 * `examples` -- need explicit import statements
 * `framework` has been moved out of pycmbs library and renamed to `benchmarking`
 * `tools` contents have been moved to `../scripts` directory as they are effectively scripts
 *  tests have to be written for ~17 pycmbs submodules

# ISSUES
 * a lot of comments in the code
 * some multitask tests can be split into single method test with single assertion statement
 * manual copying updates from pyCMBS to pycmbs does not preserve history -- not optimal solution

# QUESTIONS
 * shall `benchmarking` be a part of pycmbs or a standalone project?
 * shall `colormaps` be a regular pycmbs submodule?

# SUGGESTIONS
 * increase test coverage (only four submodules are test covered, also incompletely)


# FEATURE FREEZE?

# TESTS, TESTS, TESTS

# Github, PIP?


# NOTES 
## Two major problems: 
 * not preserving original history when moving updates from pyCMBS into pycmbs -- we need to stick to common naming/structure
 * test coverage is really small, makes pycmbs development kind of "blind"
 * lots of new feature development

## Minor problems:
 * lots of implicit statements, make hard to develop -- I don't know to which module belong classes and functions
 * test methods (e.g. test_data) often include more than one assertion in one methods -- hard to maintain
 * machine/user specific configurations needs to be taken out of package (sample config generator would help)

# ROADMAP?
 1. Pick up one naming convention/module structure
 2. Make temporary feature freeze
 3. cover code with unit tests
 4. release version 1. possibly without benchmarking  
