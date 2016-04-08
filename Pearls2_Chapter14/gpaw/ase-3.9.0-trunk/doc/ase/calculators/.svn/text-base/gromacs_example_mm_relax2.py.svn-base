""" 
    An example for using gromacs calculator in ase.
    Atom positions are relaxed.
    If qm atoms are found in index file, they are kept fixed in the relaxation.
    QM atoms are defined in index file (numbering from 1)
    in sets containing QM or Qm or qm in their name.

    A sample call:

   ./gromacs_example_mm_relax.py 1GM2_sol.gro
"""

from ase.calculators.gromacs import Gromacs

import sys
from ase.io import read

infile_name = sys.argv[1]

CALC_MM_RELAX = Gromacs(
    init_structure_file=infile_name,
    structure_file = 'gromacs_mm-relax.g96',
    force_field='oplsaa', 
    water_model='tip3p',    
    base_filename = 'gromacs_mm-relax',
    doing_qmmm = False, freeze_qm = True,
    index_filename = 'index.ndx',
    extra_mdrun_parameters = ' -nt 1 ',
    define = '-DFLEXIBLE',
    integrator = 'cg',
    nsteps = '10000',
    nstfout = '10',
    nstlog = '10',
    nstenergy = '10',
    nstlist = '10',
    ns_type = 'grid',
    pbc = 'xyz',
    rlist = '1.15',
    coulombtype = 'PME-Switch',
    rcoulomb = '0.8',
    vdwtype = 'shift',
    rvdw = '0.8',
    rvdw_switch = '0.75',
    DispCorr = 'Ener')
CALC_MM_RELAX.generate_topology_and_g96file()
CALC_MM_RELAX.generate_gromacs_run_file()
CALC_MM_RELAX.run()
