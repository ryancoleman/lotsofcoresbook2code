""" demo run for ase_qmmm_manyqm calculator """

#  ./test_ase_qmmm_manyqm.py gromacs_mm-relax.g96

from ase.calculators.gromacs import Gromacs
from ase.calculators.aims import Aims
from ase.calculators.ase_qmmm_manyqm import AseQmmmManyqm
from ase.optimize import BFGS

import sys
from ase.io.gromos import read_gromos

RUN_COMMAND = '/home/mka/bin/aims.071711_6.serial.x'
SPECIES_DIR = '/home/mka/Programs/fhi-aims.071711_6/species_defaults/light/'

LOG_FILE = open("ase-qm-mm-output.log","w")
sys.stdout = LOG_FILE

infile_name = sys.argv[1]

CALC_QM1 = Aims(charge = 0,
                xc = 'pbe',
                sc_accuracy_etot = 1e-5,
                sc_accuracy_eev = 1e-2,
                sc_accuracy_rho = 1e-5,
                sc_accuracy_forces = 1e-3,
                species_dir = SPECIES_DIR,
                run_command = RUN_COMMAND)
CALC_QM1.set(output = 'hirshfeld')

CALC_QM2 = Aims(charge = 0,
                xc = 'pbe',
                sc_accuracy_etot = 1e-5,
                sc_accuracy_eev = 1e-2,
                sc_accuracy_rho = 1e-5,
                sc_accuracy_forces = 1e-3,
                species_dir = SPECIES_DIR,
                run_command = RUN_COMMAND)
CALC_QM2.set(output = 'hirshfeld')

CALC_QM3 = Aims(charge = 0,
                xc = 'pbe',
                sc_accuracy_etot = 1e-5,
                sc_accuracy_eev = 1e-2,
                sc_accuracy_rho = 1e-5,
                sc_accuracy_forces = 1e-3,
                species_dir = SPECIES_DIR,
                run_command = RUN_COMMAND)
CALC_QM3.set(output = 'hirshfeld')

CALC_MM = Gromacs(
    init_structure_file = infile_name,
    structure_file = 'gromacs_qm.g96', \
    force_field = 'oplsaa', 
    water_model = 'tip3p',
    base_filename = 'gromacs_qm',
    doing_qmmm = True, freeze_qm = False,
    index_filename = 'index.ndx',
    define = '-DFLEXIBLE',
    integrator = 'md',
    nsteps = '0',
    nstfout = '1',
    nstlog = '1',
    nstenergy = '1',
    nstlist = '1',
    ns_type = 'grid',
    pbc = 'xyz',
    rlist = '1.15',
    coulombtype = 'PME-Switch',
    rcoulomb = '0.8',
    vdwtype = 'shift',
    rvdw = '0.8',
    rvdw_switch = '0.75',
    DispCorr = 'Ener')
CALC_MM.generate_topology_and_g96file()
CALC_MM.generate_gromacs_run_file()

CALC_QMMM = AseQmmmManyqm(nqm_regions = 3, 
                          qm_calculators = [CALC_QM1, CALC_QM2, CALC_QM3], 
                          mm_calculator = CALC_MM,
                          link_info = 'byQM')
#                         link_info = 'byFILE') 

SYSTEM = read_gromos('gromacs_qm.g96')
SYSTEM.set_calculator(CALC_QMMM)
DYN = BFGS(SYSTEM)
DYN.run(fmax = 0.05)

print('exiting fine')
LOG_FILE.close()
