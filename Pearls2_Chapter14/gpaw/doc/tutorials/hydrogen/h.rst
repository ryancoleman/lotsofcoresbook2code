===================================
Bare Coulomb potential for hydrogen
===================================


Using plane waves
=================

GPAW has a special PAW setup for hydrogen called ``ae`` for
all-electron.  It's not really a PAW setup because it doesn't make use
of any PAW magic at all --- it's just the bare Coulomb potential
(`-1/r`).

The convergence of the energy as a function of plane-wave cutoff energy
will be very slow due to the divergence of the potential at the
hydrogen nucleus and also because of the cusp in the wave function:

.. literalinclude:: h.py

You can look at the energy convergence with this command::

    $ ase-gui H.ae.txt

Let's do the same calculation with a PAW setup.  Replace the ``h.calc
=`` line with this::

    h.calc = GPAW(txt='H.paw.txt')

Now the energy is converged much quicker:

.. image:: h.png


Using a 1-d radial grid
=======================

Since the H atom is spherically symmetric, one can solve the problem
on a 1-d grid.  GPAW has a program to do this called ``aeatom.py``.
It can be used like this::

    $ python -m gpaw.atom.aeatom -p H
