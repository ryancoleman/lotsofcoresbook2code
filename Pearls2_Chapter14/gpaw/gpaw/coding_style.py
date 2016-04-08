# Copyright (C) 2008  CAMd
# Please see the accompanying LICENSE file for further information.

"""This module is an example of good coding style.

This docstring should begin with a one-line description followed by a
blank line, and then this paragraph describing in more word what kind
of functionality this module implements.

After this docstring we have import statements in this order:

1. From the Python standard library.
2. Other libraries (numpy, ase, ...).
3. GPAW stuff.
"""

from math import pi

import numpy as np
from ase.units import kJ, Hartree

from gpaw import debug
from gpaw.fd_operators import Gradient
import gpaw.mpi as mpi


class SimpleExample:
    """A simple example class.

    A headline, a blank line and then this longer description of the
    class.

    Here one could put an example of how to use the class::

      ex = SimpleExample('Test', (2, 3), int, verbose=False)
      ex.run(7, verbose=True)

    """

    def __init__(self, name, shape, dtype=float, verbose=True):
        """Create an example object.

        Again, headline, blank line, ... .  If there are many
        parameters, there should be a parameter section (see below).
        If there only a few possible arguments, then the parameter
        section can be left out and the arguments can be described in
        the section folowing the headline and blank line (see the
        `run` method).  If a method is real simple and
        self-explanatory, the docstring can be the headline only (see
        the `reset` method).

        Parameters:

        name : string
            Name of the example.
        shape:  tuple
            Shape of the ndarray.
        dtype: ndarray datatype
            The datatype of the ndarray.  Here, the description can go
            on to a second line if needed.  Make sure that the
            indentation is like shown here, and remember to end with a
            period.
        verbose: boolean
            Print information about this and that.

        Other sections:

        There can be other sections - see bolow and here:

          http://scipy.org/...

        """

        self.name = name
        if verbose:
            print(name)
        self.a = np.zeros(shape, dtype)
        self.verbose = verbose

    def method_with_long_name(self, b, out=None):
        """Do something very complicated.

        Long story with all details here ...

        Parameters:

        b : ndarray
            Add this array.
        out : ndarray
            Optional output array.

        Returns:

        The sum of ...
        """

        if out is none:
            return self.a + b
        else:
            return np.add(self.a, b, out)

    def run(self, n):
        """Do something.

        Do it n times, where n must be a positive integer.  The final
        result bla-bla is returned.
        """

        for i in range(n):
            self.a += i
            if self.verbose:
                print(self.a)

        return pi * self.a / n + 1

    def reset(self):
        """Simple method - no explanation needed."""
        self.a[:] = 0


def function(a, b):
    """Headline.

    Long story ..."""

    return a + b
