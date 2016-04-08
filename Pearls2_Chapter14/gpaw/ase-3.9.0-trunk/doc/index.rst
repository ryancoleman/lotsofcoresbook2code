=============================
Atomic Simulation Environment
=============================

The Atomistic Simulation Environment (ASE) is a set of tools and Python_
modules for setting up, manipulating, running, visualizing and analyzing
atomistic simulations.  The code is freely available under the `GNU LGPL
license`_.

.. _Python: http://www.python.org
.. _GNU LGPL license: https://wiki.fysik.dtu.dk/ase/licenseinfo.html

Alphabetical list of all modules:
    
.. list-table::

  * - :mod:`~ase.atom`
    - :mod:`~ase.atoms`
    - :mod:`~ase.calculators`
    - :mod:`~ase.constraints`
  * - :mod:`~ase.db`
    - :mod:`~ase.dft`
    - :mod:`~ase.data`
    - :mod:`~ase.ga`
  * - :mod:`~ase.gui`
    - :mod:`~ase.infrared`
    - :mod:`~ase.io`
    - :mod:`~ase.lattice`
  * - :mod:`~ase.md`
    - :mod:`~ase.neb`
    - :mod:`~ase.optimize`
    - :mod:`~ase.parallel`
  * - :mod:`~ase.phonons`
    - :mod:`~ase.lattice.spacegroup`
    - :mod:`~ase.structure`
    - :mod:`~ase.lattice.surface`
  * - :mod:`~ase.transport`
    - :mod:`~ase.thermochemistry`
    - :mod:`~ase.units`
    - :mod:`~ase.utils`
  * - :mod:`~ase.vibrations`
    - :mod:`~ase.visualize`
    - :mod:`~ase.visualize.vtk`
    -

:mod:`Calculators <ase.calculators>`:

|abinit| |Asap| |CASTEP| |dftb| |elk| |exciting| |EMT| |fhi-aims| |fleur|
|gpaw| |gromacs| |hotbit| |jacapo| |jdftx| |lammps| |nwchem| |siesta|
|turbomole| |vasp| Gaussian_ Mopac_


.. |abinit| image:: _static/abinit.png
   :target: ase/calculators/abinit.html
   :align: middle
.. |Asap| image:: _static/asap.png
   :target: http://wiki.fysik.dtu.dk/asap
   :align: middle
.. |CASTEP| image:: _static/castep.png
   :target: ase/calculators/castep.html
   :align: middle
.. |elk| image:: _static/elk.png
   :target: http://elk.sourceforge.net/
   :align: middle
.. |EMT| image:: _static/emt.png
   :target: ase/calculators/emt.html
   :align: middle
.. |exciting| image:: _static/exciting.png
   :target: ase/calculators/exciting.html
   :align: middle
.. |dftb| image:: _static/dftb.png
   :target: ase/calculators/dftb.html
   :align: middle
.. |fhi-aims| image:: _static/fhi-aims.png
   :target: ase/calculators/FHI-aims.html
   :align: middle
.. |fleur| image:: _static/fleur.png
   :target: ase/calculators/fleur.html
   :align: middle
.. |gpaw| image:: _static/gpaw.png
   :target: http://wiki.fysik.dtu.dk/gpaw
   :align: middle
.. |gromacs| image:: _static/gromacs.png
   :target: http://www.gromacs.org/
   :align: middle
.. |hotbit| image:: _static/hotbit.png
   :target: https://trac.cc.jyu.fi/projects/hotbit
   :align: middle
.. |jacapo| image:: _static/jacapo.png
   :target: ase/calculators/jacapo.html
   :align: middle
.. |jdftx| image:: _static/jdftx.png
   :target: http://sourceforge.net/p/jdftx/wiki/ASE%20Interface
   :align: middle
.. |lammps| image:: _static/lammps.png
   :target: ase/calculators/lammps.html
   :align: middle
.. |nwchem| image:: _static/nwchem.png
   :target: http://www.nwchem-sw.org
   :align: middle
.. |siesta| image:: _static/siesta.png
   :target: ase/calculators/siesta.html
   :align: middle
.. |turbomole| image:: _static/tm_logo_l.png
   :target: ase/calculators/turbomole.html
   :align: middle
.. |vasp| image:: _static/vasp.png
   :target: ase/calculators/vasp.html
   :align: middle

.. _Gaussian: http://www.gaussian.com/
.. _Mopac: http://openmopac.net/

.. _news:

News
====

* :ref:`ASE version 3.8.0 <releasenotes>` released (22 October 2013).

* :ref:`ASE version 3.7.0 <releasenotes>` released (13 May 2013).

* :ref:`ASE version 3.6.0 <releasenotes>` released (24 February 2012).

* Bugfix release: :ref:`ASE version 3.5.1 <releasenotes>` (24 May 2011).

* :ref:`ASE version 3.5.0 <releasenotes>` released (13 April 2011).

* :ref:`ASE version 3.4.1 <download_and_install>` released (11 August 2010).

* :ref:`ASE version 3.4 <download_and_install>` released (23 April 2010).

* :ref:`ASE version 3.3 <download_and_install>` released (11 January 2010).

* :ref:`ASE version 3.2 <download_and_install>` released (4 September 2009).

* ASE has reached revision 1000 (16 July 2009).

* :ref:`ASE version 3.1.0 <download_and_install>` released (27 March 2009).

* Improved :mod:`ase.vibrations` module: More accurate and
  possibility to calculate :mod:`infrared intensities <ase.infrared>` (13
  March 2009).

* :ref:`ASE version 3.0.0 <download_and_install>` released (13 November 2008).

* Asap_ version 3.0.2 released (15 October 2008).

* An experimental abinit interface released (9 June 2008).

* Thursday April 24 will be ASE documentation-day.  Ten people from
  CAMd/Cinf will do a "doc-sprint" from 9 to 16.  (17 Apr 2008)

* The new ASE-3.0 Sphinx_ page is now up and running!  (2 Apr 2008)

* A beta version of the new ASE-3.0 will be used for the
  electronic structure course at CAMd_.  (10 Jan 2008)


.. _Sphinx: http://sphinx.pocoo.org
.. _Asap: http://wiki.fysik.dtu.dk/asap
.. _CAMd: http://www.camd.dtu.dk
