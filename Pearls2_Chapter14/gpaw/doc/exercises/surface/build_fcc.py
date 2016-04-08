"""Module for building atomic structures"""

from math import sqrt
from ase import Atoms, Atom
from ase.visualize import view


def fcc100(symbol, a, layers, L):
    """Build an fcc(100) surface

    Parameters
    ----------
    symbol: string
        Chemical symbol ('H', 'Li', ...).
    a: float
        Lattice constant.
    layers: int
        Number of layers.
    L: float
        Height of unit cell.

    """

    # Distance between atoms:
    d = a / sqrt(2)

    # Distance between layers:
    z = a / 2.

    assert L > layers * z, 'Unit cell too small!'

    # Start with an empty Atoms object:
    atoms = Atoms(cell=(d, d, L),
                  pbc=(True, True, False))

    # Fill in the atoms:
    for n in range(layers):
        position = [d / 2 * (n % 2),
                    d / 2 * (n % 2),
                    n * z]
        atoms.append(Atom(symbol, position))

    atoms.center(axis=2)

    return atoms

if __name__ == '__main__':

    fcc = fcc100('Al', 4.0, 4, 15.0)
    view(fcc * (4, 4, 2))
