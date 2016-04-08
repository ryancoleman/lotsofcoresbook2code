Bethe-Salpeter equation
=======================

:Who:
    Jun

Documentation can be found at :ref:`bse`.

The BSE module can work very well for small bulk systems. See gpaw/test/bse_silicon.py as an example.
Recently I have implemented using scalapack library to diagonalize the BSE matrix, so in priciple, it can also work for relatively bigger systems. Bigger tests coming soon.
