"""Module for building atomic structures"""

from ase import Atom, Atoms


def bcc100(symbol, a, layers, L):
    """Build a bcc(100) surface

    symbol: chemical symbol ('H', 'Li', ...)
    a     : lattice constant
    layers: number of layers
    L     : height of unit cell"""

    a = float(a)

    # Distance between layers:
    z = a / 2

    assert L > layers * z, 'Unit cell too small!'

    # Start with an empty Atoms object with an orthorhombic unit cell:
    atoms = Atoms(pbc=(True, True, False), cell=(a, a, L))

    # Fill in the atoms:
    for n in range(layers):
        position = [a / 2 * (n % 2), a / 2 * (n % 2), n * z]
        atoms.append(Atom(symbol, position))

    atoms.center(axis=2)
    return atoms

if __name__ == '__main__':
    bcc = bcc100('Al', 4.0, 4, 15.0)
    from ase.visualize import view
    view(bcc * (4, 4, 2))
