New PAW setups
==============

:Who:
    Marcin and Jens Jørgen

The goal is to produce new PAW setup for all elements.  Compared to our
old collection of PAW setups we will focus on:

* higher accuracy - more semi-core states might be added
* reduced eggbox error
* faster convergence of energy with number of grid-points - if possible
* tesing the setups more carfully agains a bigger set of all-electron results

The code is in :svn:`~gpaw/gpaw/atom/generator2.py` and it is based on
a new and more robust atomic solver: :svn:`~gpaw/gpaw/atom/aeatom.py`.
