.. module:: ase.calculators.ase_qmmm_manyqm

===============
ase_qmmm_manyqm
===============

Introduction
============

This an general interface to run QM/MM calculations with a ase-QM and ase-MM calculator. Currently QM can be FHI-aims and MM can be Gromacs.

QM and MM region can be covalently bonded. In principle you can cut wherever you like, but other cutting schemes are probably better than others. The standard recommendation is to cut only the C-C bonds.

There can be many QM regions embedded in the MM system.


Ase_qmmm_manyqm Calculator
==========================
QM/MM interface with QM=FHI-aims, MM=gromacs

QM could be any QM calculator, but you need to read in QM-atom charges
from the QM program (in method 'get_qm_charges'). 

One can have many QM regions, each with a different calculator.
There can be only one MM calculator, which is calculating the whole
system. 


Non-bonded interactions:
""""""""""""""""""""""""
Generally energies and forces are treated by:
  - Within the same QM-QM: by QM calculator
  - MM-MM: by MM calculator
  - QM-MM: by MM using MM vdw parameters and QM charges.
  - Different QM-different QM: by MM using QM and MM charges and MM-vdw parameters

The Hirschfeld charges (or other atomic charges) 
on QM atoms are calculated by QM in a H terminated cluster in vacuum. 
The charge of QM atom next to MM atom (edge-QM-atom) 
and its H neighbors are set as in the classical force field. 

The extra(missing) charge results from: 
  - linkH atoms
  - The edge-QM atoms, and their singly bonded neighbors (typically -H or =O). These have their original MM charges.
  - From the fact that the charge of the QM fraction is not usually an integer when using the original MM charges.
  - The extra/missing charge is added equally to all QM atoms (not being linkH and not being edge-QM-atom or its singly bonded neighbor) so that the total charge of the MM-fragment involving QM atoms will be the same as in the original MM-description.

Vdw interactions are calculated by MM-gromacs for MM and MM-QM interactions.
The QM-QM vdw interaction s could be done by the FHI-aims if desired
(by modifying the input for QM-FHI-aims input accordingly.

Bonded interactions:
""""""""""""""""""""
E = E_qm(QM-H) + E_mm(ALL ATOMS), where
  - E_qm(QM-H): qm energy of H terminated QM cluster(s) 
  - E_mm(ALL ATOMS): MM energy of all atoms, except for terms in which all MM-interacting atoms are in the same QM region.

Forces do not act on link atoms but they are positioned by scaling.
Forces on link atoms are given to their QM and MM neighbors by chain rule.
(see Top Curr Chem (2007) 268: 173–290, especially pages 192-194).
At the beginning of a run, the optimal edge-qm-atom-linkH bond length(s)
is (are) calculated by QM in 'get_eq_qm_atom_link_h_distances'
or they are read from a file (this is determined by the argument 'link_info').

The covalent bonds between QM and MM are found automatically based on 
ase-covalent radii in neighbor search. Link atoms are positioned 
automatically (according to . 

Questions & Comments markus.kaukonen@iki.fi

Interesting applications could be cases when we need two or more 
QM regions. For instance two redox centers in a protein, 
cathode and anode of a fuel cell ... you name it!

Arguments
==========

==================   =========  ==============  =============================
keyword              type       default value   description
==================   =========  ==============  =============================
``nqm_regions``      ``int``                    The number of QM regions

``qm_calculators``   ``list``                   List of members of a Class 
                                                defining a ase-QM calculator 
                                                for each QM region

``mm_calculator``    ``Calc``                   A member of a Class 
                                                defining a ase-MM calculator 

``link_info``        ``str``    ``'byQM'``      Can be either
                                                'byQM': the 
                                                edge_qm_atom-link_h_atom 
                                                distances are calculated 
                                                by QM
                                                or 
                                                'byFile':the 
                                                edge_qm_atom-link_h_atom 
                                                distances are read from 
                                                a file
==================   =========  ==============  =============================


Example
=======

1. Prepare classical input with Gromacs

  - Get THE INDEPENDENT STRUCTURE OF THE ANTITRYPTIC REACTIVE SITE LOOP OF BOWMAN-BIRK INHIBITOR AND SUNFLOWER TRYPSIN INHIBITOR-1 (pdb code 1GM2) from pdb data bank http://www.rcsb.org/pdb/home/home.do, name it to 1GM2.pdb

  - In file 1GM2.pdb take only MODEL1
    replace 'CYS ' by 'CYS2'
    in order to get deprotonated CYS-CYS  bridge.

  - Generate gromacs coordinates and topology 

    >>> pdb2gmx -ff oplsaa -f 1GM2.pdb -o 1GM2.gro -p 1GM2.top -water tip3p -ignh


  - Generate the simulation box (cubic is not the most efficient one...)

    >>> editconf -bt cubic -f 1GM2.gro -o 1GM2_box.gro -c -d 0.2

  - Solvate the protein in water box

    >>> genbox -cp 1GM2_box.gro -cs spc216.gro -o 1GM2_sol.gro -p 1GM2.top

  - Generate index file (needed later also for defining QM atoms)

    >>> make_ndx -f 1GM2_sol.gro -o index.ndx
    >>> select 'q' <ENTER>

  - add to first lines to index.ndx (we define here 3 QM regions, indexing from 1)::

      [qm_ss]
      18 19 20 21 139 140 141 142

      [qm_thr]
      28 29 30 31 32 33 34 35

      [qm_ser]
      64 65 66 67 68 


2. Relax MM system by fixing QM atoms in their (crystallographic?) positions (see documentation for Gromacs MM calculator).

.. literalinclude:: gromacs_example_mm_relax2.py

3. QM/MM relaxation (all atoms are free)

.. literalinclude:: test_ase_qmmm_manyqm.py

