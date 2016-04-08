.. _setups:

=================
Atomic PAW Setups
=================

A setup is to the PAW method what a pseudo-potential is to the
pseudo-potential method.  All available setups are contained in this
tar-file: gpaw-setups-0.9.11271.tar.gz_.  There are setups for the LDA,
PBE, revPBE, RPBE and GLLBSC functionals.  Install them as described
in the :ref:`installationguide_setup_files`.  The setups are stored as
compressed :ref:`pawxml` files.

Setup releases
==============

===========  =======  ========  =============================
Date         Version  Revision  Tarfile                
===========  =======  ========  =============================
Mar 27 2014  0.9      11271     gpaw-setups-0.9.11271.tar.gz_
Oct 26 2012  0.9      9672      gpaw-setups-0.9.9672.tar.gz_
Apr 13 2011  0.8      7929      gpaw-setups-0.8.7929.tar.gz_
Apr 19 2010  0.6      6300      gpaw-setups-0.6.6300.tar.gz_
Jul 22 2009  0.5      3574      gpaw-setups-0.5.3574.tar.gz_
===========  =======  ========  =============================

.. _gpaw-setups-0.9.11271.tar.gz:
    https://wiki.fysik.dtu.dk/gpaw-files/gpaw-setups-0.9.11271.tar.gz

.. _gpaw-setups-0.9.9672.tar.gz:
    https://wiki.fysik.dtu.dk/gpaw-files/gpaw-setups-0.9.9672.tar.gz

.. _gpaw-setups-0.8.7929.tar.gz:
    https://wiki.fysik.dtu.dk/gpaw-files/gpaw-setups-0.8.7929.tar.gz

.. _gpaw-setups-0.6.6300.tar.gz:
    https://wiki.fysik.dtu.dk/gpaw-files/gpaw-setups-0.6.6300.tar.gz

.. _gpaw-setups-0.5.3574.tar.gz:
    https://wiki.fysik.dtu.dk/gpaw-files/gpaw-setups-0.5.3574.tar.gz


Periodic table
==============

=== === === === === === === === === === === === === === === === === ===
H_                                                                  He_
Li_ Be_                                         B_  C_  N_  O_  F_  Ne_ 
Na_ Mg_                                         Al_ Si_ P_  S_  Cl_ Ar_  
K_  Ca_ Sc_ Ti_ V_  Cr_ Mn_ Fe_ Co_ Ni_ Cu_ Zn_ Ga_ Ge_ As_ Se_ Br_ Kr_
Rb_ Sr_ Y_  Zr_ Nb_ Mo_ Tc  Ru_ Rh_ Pd_ Ag_ Cd_ In_ Sn_ Sb_ Te_ I_  Xe_ 
Cs_ Ba_ La_ Hf_ Ta_ W_  Re_ Os_ Ir_ Pt_ Au_ Hg_ Tl_ Pb_ Bi_ Po  At  Rn_ 
=== === === === === === === === === === === === === === === === === ===

See also `NIST Atomic Reference Data`_, `Computational Chemistry
Comparison and Benchmark DataBase`_, `Dacapo pseudo potentials`_, and
`Vasp pseudo potentials`_.

.. _NIST Atomic Reference Data: http://physics.nist.gov/PhysRefData/DFTdata/Tables/ptable.html
.. _Computational Chemistry Comparison and Benchmark DataBase: http://srdata.nist.gov/cccbdb/
.. _Dacapo pseudo potentials: https://wiki.fysik.dtu.dk/dacapo/Pseudopotential_Library
.. _Vasp pseudo potentials: http://cms.mpi.univie.ac.at/vasp/vasp/Pseudopotentials_supplied_with_VASP_package.html


.. toctree::
   :maxdepth: 2

   g2_1
   dcdft
   molecule_tests
   bulk_tests
   generation_of_setups
   pawxml

.. from gpaw.atom.configurations import parameters
   for s in parameters:
       print '.. %3s: %2s.html' % ('_' + s, s)

.. from gpaw.atom.configurations import parameters
   from ase.data import atomic_numbers
   from ase.data import chemical_symbols
   anums=[]
   for s in parameters:
       anums.append(atomic_numbers[s])
   anums.sort()
   for a in anums:
       s = chemical_symbols[a]
       print '.. %3s: %2s.html' % ('_' + s, s)

..  _H:  H.html
.. _He: He.html
.. _Li: Li.html
.. _Be: Be.html
..  _B:  B.html
..  _C:  C.html
..  _N:  N.html
..  _O:  O.html
..  _F:  F.html
.. _Ne: Ne.html
.. _Na: Na.html
.. _Mg: Mg.html
.. _Al: Al.html
.. _Si: Si.html
..  _P:  P.html
..  _S:  S.html
.. _Cl: Cl.html
.. _Ar: Ar.html
..  _K:  K.html
.. _Ca: Ca.html
.. _Sc: Sc.html
.. _Ti: Ti.html
..  _V:  V.html
.. _Cr: Cr.html
.. _Mn: Mn.html
.. _Fe: Fe.html
.. _Co: Co.html
.. _Ni: Ni.html
.. _Cu: Cu.html
.. _Zn: Zn.html
.. _Ga: Ga.html
.. _Ge: Ge.html
.. _As: As.html
.. _Se: Se.html
.. _Br: Br.html
.. _Kr: Kr.html
.. _Rb: Rb.html
.. _Sr: Sr.html
..  _Y: Y.html
.. _Zr: Zr.html
.. _Nb: Nb.html
.. _Mo: Mo.html
.. _Ru: Ru.html
.. _Rh: Rh.html
.. _Pd: Pd.html
.. _Ag: Ag.html
.. _Cd: Cd.html
.. _In: In.html
.. _Sn: Sn.html
.. _Sb: Sb.html
.. _Te: Te.html
..  _I:  I.html
.. _Xe: Xe.html
.. _Cs: Cs.html
.. _Ba: Ba.html
.. _La: La.html
.. _Hf: Hf.html
.. _Ta: Ta.html
..  _W:  W.html
.. _Re: Re.html
.. _Os: Os.html
.. _Ir: Ir.html
.. _Pt: Pt.html
.. _Au: Au.html
.. _Hg: Hg.html
.. _Tl: Tl.html
.. _Pb: Pb.html
.. _Bi: Bi.html
.. _Rn: Rn.html

.. toctree::
   :glob:
   :maxdepth: 1

   [A-Z]*
