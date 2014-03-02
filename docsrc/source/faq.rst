FAQ
===

How to mask a specific area in a whole datafile?
------------------------------------------------

Let us assume that you have a mask which specifies some region in your data.
You now want to apply a mask to that region and probably even only extract the
region you are interested in. This implies two steps:

a. apply the mask to the data
b. cut the data using its bounding box.

Doing so is straight forward in pyCMBS as the following examples illustrates.


.. plot:: ../../pycmbs/examples/howto_mask_region.py
    :include-source:


