""" An example for using gromacs calculator in ase.
    Atom positions are relaxed.
    A sample call:

   python ./gromacs_example_mm_relax.py his.pdb
"""

from ase.calculators.gromacs import Gromacs

import sys
from ase.io import read

infile_name = sys.argv[1]

CALC_MM_RELAX = Gromacs(clean=True)
CALC_MM_RELAX.set_own_params_runs(
    'extra_pdb2gmx_parameters','-ignh')
CALC_MM_RELAX.set_own_params_runs(
    'init_structure',infile_name)
CALC_MM_RELAX.generate_topology_and_g96file()
CALC_MM_RELAX.write_input()
CALC_MM_RELAX.set_own_params_runs(
    'extra_editconf_parameters','-bt cubic -c -d 0.8')
CALC_MM_RELAX.run_editconf()
CALC_MM_RELAX.set_own_params_runs(
    'extra_genbox_parameters','-cs spc216.gro')
CALC_MM_RELAX.run_genbox()
CALC_MM_RELAX.generate_gromacs_run_file()
CALC_MM_RELAX.run()
