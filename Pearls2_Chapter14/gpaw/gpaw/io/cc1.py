from math import pi, cos, sin, sqrt, acos

import numpy as np

from ase.atoms import Atoms
from ase.parallel import paropen


def read_cc1(fileobj, index=-1):
    """Read Chem3D(CC1) format"""
    if isinstance(fileobj, str):
        fileobj = open(fileobj)

    lines = fileobj.readlines()
    L1 = lines[0].split()
    del lines[0]
    natoms = int(L1[0])

    images = []
    while len(lines) > 0:
        positions = []
        symbols = []
        for line in lines[:natoms]:
            symbol, number, x, y, z = line.split()[:5]
            symbols.append(symbol)
            positions.append([float(x), float(y), float(z)])
        images.append(Atoms(symbols=symbols, positions=positions))
        del lines[:natoms + 1]
    return images[index]

