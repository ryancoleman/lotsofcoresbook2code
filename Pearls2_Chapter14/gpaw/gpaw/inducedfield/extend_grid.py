import numpy as np
from gpaw.grid_descriptor import GridDescriptor
from ase.units import Bohr

def extend_grid(gd, N_cd):
    N_cd = np.array(N_cd)
    
    N_c = gd.N_c + N_cd.sum(axis=1)
    cell_cv = gd.h_cv * N_c
    
    move_c = gd.get_grid_spacings() * N_cd[:,0]
    
    egd = GridDescriptor(N_c, cell_cv, gd.pbc_c, gd.comm)
    egd.extend_N_cd = N_cd
    
    return egd, cell_cv*Bohr, move_c*Bohr


def extend_array(d_g, gd, d_e, egd):
    
    big_d_g = gd.collect(d_g)
    big_d_e = egd.collect(d_e)
    
    N_cd = egd.extend_N_cd
    
    if gd.comm.rank == 0:
        N1_c = N_cd[:,0]
        N2_c = N1_c + gd.N_c - 1 # implicit zero
        big_d_e[N1_c[0]:N2_c[0],N1_c[1]:N2_c[1],N1_c[2]:N2_c[2]] = big_d_g
    egd.distribute(big_d_e, d_e)


def deextend_array(d_g, gd, d_e, egd):
    
    big_d_g = gd.collect(d_g)
    big_d_e = egd.collect(d_e)
    
    N_cd = egd.extend_N_cd
    
    if gd.comm.rank == 0:
        N1_c = N_cd[:,0]
        N2_c = N1_c + gd.N_c - 1 # implicit zero
        big_d_g[:] = big_d_e[N1_c[0]:N2_c[0],N1_c[1]:N2_c[1],N1_c[2]:N2_c[2]] 
    gd.distribute(big_d_g, d_g)
    
    
def move_atoms(atoms, move_c):
    pos_a = atoms.get_positions()
    for pos in pos_a:
        pos += move_c
    atoms.set_positions(pos_a)
