.. module:: ase.calculators.gromacs

=======
Gromacs
=======

Introduction
============

Gromacs is a free classical molecular dynamics package. It is mainly 
used in modeling of biological systems. It is part of the 
ubuntu-linux distribution.
http://www.gromacs.org/

.. warning:: 1) Ase-Gromacs calculator works properly only with gromacs version 4.5.6 or a newer one. (fixed bug related to .g96 format)

.. warning:: 2) It only makes sense to use ase-gromacs for qm/mm or for testing. For pure MM production runs the native gromacs is much much faster (at the moment ase-gromacs has formatted io using .g96 format which is slow).

Gromacs Calculator
==================
This ASE-interface is a preliminary one and it is VERY SLOW so 
do not use it for production runs. It is here because of 
we'll have a QM/MM calculator which is using gromacs as the 
MM part.

For example: (setting for the MM part of a QM/MM run, 
parameter '-nt 1' for serial run)::

  CALC_MM = Gromacs(doing_qmmm = True)
  CALC_MM.set_own_params_runs('extra_mdrun_parameters', ' -nt 1 ')

Here default values for MM input are::
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
    DispCorr = 'Ener'

The input values can be changed by::
    CALC_MM.set_own_params(
        'nsteps','99999' 'Number of steps')

The arguments for gromacs programs can be changed by::
    CALC_MM.set_own_params_runs(
        'extra_pdb2gmx_parameters','-ignh')
    CALC_MM.set_own_params_runs(
        'init_structure','1GM2.pdb')
    CALC_MM.set_own_params_runs(
        'extra_mdrun_parameters', ' -nt 1 ')
    CALC_MM.set_own_params_runs(
        'extra_grompp_parameters', ' ')
    CALC_MM.set_own_params_runs(
        'extra_editconf_parameters', ' ')
    CALC_MM.set_own_params_runs(
        'extra_genbox_parameters', ' ')


Parameters
==========
The description of the parameters can be found in the Gromacs manual:
http://www.gromacs.org/Documentation/Manual

and extra (ie. non-gromacs) parameter: 

do_qmmm: logical      (default False)
    If true we run only single step of gromacs 
    (to get MM forces and energies in QM/MM)
freeze_qm: logical    (default False)
    If true, the qm atoms will be kept fixed
    (The list of qm atoms is taken from file 'index_filename', below)
clean:       logical  (default True)
    If true old gromacs files are cleaned
force_field: str      (default oplsaa)
    Name of the force field for gromacs
water_model: str      (default tip3p)
    Name of the water model for gromacs


Environmental variables:
========================
  - GMXCMD the name of the main gromacs executable (usually 'mdrun').
    If GMXCMD is not set gromacs test is not run, but in the calculator 
    works using 'mdrun'.
  - GMXCMD_PREF prefix for all gromacs commands (default '')
  - GMXCMD_POST postfix (ie suffix) for all gromacs commands (default '') 
    
Example: MM-only geometry optimization of a histidine molecule
==============================================================
THIS IS NOT A PROPER WAY TO SETUP YOUR MD SIMULATION.
THIS IS JUST A DEMO THAT DOES NOT CRASH. 
(the box size should be iterated by md-NTP, total charge should be 0).

Initial pdb coordinates (file his.pdb):

.. literalinclude:: his.pdb

The sample file for relaxation:

.. literalinclude:: gromacs_example_mm_relax.py
