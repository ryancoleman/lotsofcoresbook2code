===============================================================
GPAW: DFT and beyond within the projector-augmented wave method
===============================================================

GPAW is a density-functional theory (DFT) Python_ code based on the
projector-augmented wave (:ref:`PAW <literature>`) method and the
atomic simulation environment (ASE_).  It uses real-space uniform
grids and multigrid methods, atom-centered basis-functions or plane-waves.
Read more about the :ref:`algorithms` used.

.. _Python: http://www.python.org
.. _ASE: https://wiki.fysik.dtu.dk/ase


Quick links to all features:
    
.. list-table::
    
    * - :ref:`Plane-waves <manual_mode>`
      - :ref:`Finite-difference <manual_stencils>`
      - :ref:`LCAO <lcao>`
    * - :ref:`XC-functionals <xc>`
      - :ref:`DFT+U <hubbardu>`
      - :ref:`GLLB-SC <bandgab>`
    * - :ref:`DOS <pdos>`
      - :ref:`STM <stm tutorial>`
      - :ref:`Wannier functions <wannier>`
    * - :ref:`delta-SCF <dscf>`
      - :ref:`XAS <xas>`
      - :ref:`Jellium <jellium>`
    * - :ref:`TDDFT <timepropagation>`
      - :ref:`LRTDDFT (molecules) <lrtddft>`
      - :ref:`LRTDDFT (extended systems) <df_theory>`
    * - :ref:`Transport <transport_exercise>`
      - :ref:`NEGF-transport <transport>`
      - :ref:`Keldysh GF-transport ... <keldyshgf>`
    * - :ref:`RPA-correlation <rpa>`
      - :ref:`GW <gw_theory>`
      - :ref:`BSE <bse>`
    * - :ref:`Parallelization <parallel_runs>`
      -
      -

Watch us on GPAWTV
==================

.. include:: ytp.txt

.. _news:

News
====

* :ref:`GPAW version 0.10 <releasenotes>` released (Apr 8 2014)

* GPAW is part of the `PRACE Unified European Application Benchmark Suite`_
  (October 17 2013)

* May 21-23, 2013: :ref:`GPAW workshop <workshop>` at the Technical
  University of Denmark (Feb 8 2013)

* Prof. Häkkinen has received `18 million CPU hour grant`_ for GPAW based 
  research project (Nov 20 2012)

* A new :ref:`setups` bundle released (Oct 26 2012)

* :ref:`GPAW version 0.9 <releasenotes>` released (March 7 2012)

* Help!  The :ref:`todolist` has been updated.  Maybe there is
  something you can do (May 25 2011)
 
* :ref:`GPAW version 0.8 <releasenotes>` released (May 25 2011)

* GPAW is part of benchmark suite for `CSC's supercomputer procurement`_ 
  (Apr 19 2011)

* New features: Calculation of the linear :ref:`dielectric response
  <df_theory>` of an extended system (RPA and ALDA kernels) and
  calculation of :ref:`rpa` (Mar 18 2011)

* Massively parallel GPAW calculations presented at `PyCon 2011`_.
  See William Scullin's talk here: `Python for High Performance
  Computing`_ (Mar 12 2011)

* :ref:`GPAW version 0.7.2 <releasenotes>` released (Aug 13 2010)

* :ref:`GPAW version 0.7 <releasenotes>` released (Apr 23 2010)

* GPAW is :math:`\Psi_k` `scientific highlight of the month`_ (Apr 3 2010)

* A third GPAW code sprint was successfully hosted at CAMD (Oct 20 2009)

* :ref:`GPAW version 0.6 <releasenotes>` released (Oct 9 2009)

* `QuantumWise <http://www.quantumwise.com>`_ adds GPAW-support to
  `Virtual NanoLab`_ (Sep 8 2009)

* Join the new IRC channel ``#gpaw`` on FreeNode (Jul 15 2009)

* :ref:`GPAW version 0.5 <releasenotes>` released (Apr 1 2009)

* A new :ref:`setups` bundle released (Mar 27 2009)

* A second GPAW code sprint was successfully hosted at CAMD (Mar 20 2009)

* :ref:`GPAW version 0.4 <releasenotes>` released (Nov 13 2008)

* The :ref:`exercises` are finally ready for use in the `CAMd summer
  school 2008`_ (Aug 15 2008)

* This site is now powered by Sphinx_ (Jul 31 2008)

* GPAW is now based on numpy_ instead of of Numeric (Jan 22 2008)

* :ref:`GPAW version 0.3 <releasenotes>` released (Dec 19 2007)

* CSC_ is organizing a `GPAW course`_: "Electronic structure
  calculations with GPAW" (Dec 11 2007)

* The `code sprint 2007`_ was successfully finished (Nov 16 2007)

* The source code is now in the hands of :ref:`svn` and Trac_ (Okt 22 2007)

* A GPAW Sprint will be held on November 16 in Lyngby (Okt 18 2007)

* Work on atomic basis-sets begun (Sep 25 2007)

.. _numpy: http://numpy.scipy.org/
.. _CSC: http://www.csc.fi
.. _GPAW course: http://www.csc.fi/english/csc/courses/archive/gpaw-2008-01
.. _Trac: https://trac.fysik.dtu.dk/projects/gpaw
.. _Sphinx: http://sphinx.pocoo.org
.. _CAMd summer school 2008: http://www.camd.dtu.dk/English/Events/CAMD_Summer_School_2008/Programme.aspx
.. _code sprint 2007: http://www.dtu.dk/Nyheder/Nyt_fra_Institutterne.aspx?guid={38B92D63-FB09-4DFA-A074-504146A2D678}
.. _Virtual NanoLab: http://www.quantumwise.com/products/12-products/28-atk-se-200906#GPAW
.. _scientific highlight of the month: http://www.psi-k.org/newsletters/News_98/Highlight_98.pdf
.. _pycon 2011: http://us.pycon.org/2011/schedule/presentations/226/
.. _Python for High Performance Computing: http://pycon.blip.tv/file/4881240/
.. _CSC's supercomputer procurement: http://www.csc.fi/english/pages/hpc2011
.. _18 million CPU hour grant: http://www.prace-ri.eu/PRACE-5thRegular-Call
.. _PRACE Unified European Application Benchmark Suite: http://www.prace-ri.eu/ueabs
