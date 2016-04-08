from ase.units import Bohr
from gpaw.analyse.vdwradii import vdWradii

# data from A. Bondi, J. Phys. Chem. 68 (1964) 441
data_Bondi = { # units Anstrom
    'He' : 1.40,
    'Ne' : 1.54,
    'Ar' : 1.88,
    'Kr' : 2.02,
    'Xe' : 2.16
}
# data from Felix Hanke (FHI-AIMS ?)
data_Hanke = { # units Anstrom
    'H' : 1.640449351,
    'C' : 1.8997461838999998,
    'N' : 1.7674518813999998,
    'O' : 1.6880752998999997,
    'Cu': 1.9897063095999996,
}    
for symbol in ['He', 'Ne', 'Ar', 'Kr', 'H', 'C', 'N', 'O', 'Cu']:
    R = vdWradii([symbol], 'PBE')[0]
    if symbol in data_Bondi:
        Rref = data_Bondi[symbol]
    else:
        Rref = data_Hanke[symbol]
    error = abs(R - Rref)
    print("symbol, R, Rref, error:", symbol, R, Rref, error)
    assert(error < 0.05)
