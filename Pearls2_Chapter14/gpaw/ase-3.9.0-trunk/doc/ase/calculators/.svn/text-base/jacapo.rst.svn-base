.. module::  ase.calculators.jacapo
   :synopsis: ASE python interface for Dacapo

==========================================================
Jacapo - ASE python interface for Dacapo
==========================================================

Introduction
============

Jacapo_ is an ASE interface for Dacapo_ that is fully compatible with ASE. It
replaces the old Dacapo interface using Numeric python and ASE2.
The code was originally developed by John Kitchin and detailed documentation
as well as many examples are available online:

http://gilgamesh.cheme.cmu.edu/doc/software/jacapo/index.html

Jacapo is included as an optional calculator in ASE and small differences to the
above documentation may occur, and the documentation is no longer maintained.

.. _Jacapo: http://gilgamesh.cheme.cmu.edu/doc/software/jacapo/index.html
.. _Dacapo: http://wiki.fysik.dtu.dk/dacapo

Jacapo calculator
================= 

The Jacapo interface is automatically installed with ase and can be imported using::

  from ase.calculators.jacapo import Jacapo

(You will need to have a working installation of Dacapo, however.)

.. class:: Jacapo()
    
Here is a list of available keywords to initialize the calculator:

============== ============ =====================================
keyword        type         description
============== ============ =====================================
``nc``         ``str``      Output NetCDF file, or input file if nc already exists.
``outnc``      ``str``      Output file. By default equal to nc.
``atoms``      ``object``   ase atoms object
``pw``         ``float``    Planewave cutoff in eV
``dw``         ``float``    Density cutoff in eV
``xc``         ``str``      Exchange-correlation functional. One of ['PZ','VWN','PW91','PBE','RPBE','revPBE']
``nbands``     ``int``      Number of bands
``ft``         ``float``    Fermi temperature
``kpts``       ``list``     K-point grid, e.g. kpts = (2,2,1)
``spinpol``    ``boolean``  Turn on/off spin-polarization
``fixmagmom``  ``str``      Magnetic moment of the unit cell
``symmetry``   ``boolean``  Turn on/off symmetry reduction
``stress``     ``boolean``  Turn on/off stress calculation
``dipole``     ``boolean``  Turn on/off dipole correction
``ados``       ``dict``     Atom-projected density of states
``stay_alive`` ``boolean``  Turn on/off stay alive
``debug``      ``int``      Set debug level (0=off, 10=extreme)
``deletenc``   ``boolean``  If the nc file exists, delete it (to ensure a fresh run). Default is False.
============== ============ =====================================

Example
=======

Here is an example of how to calculate the total energy of a H atom.

.. warning:: This is an example only - the parameters are not physically meaningful!

.. literalinclude:: ../../../ase/test/jacapo/jacapo.py
   :start-after: os
   :end-before: os.system
        
Note that all calculator parameters should be set in the calculator definition
itself. Do not attempt to use the calc.set_* commands as they are intended to
be internal to the calculator. Note also that Dacapo can only operate with
periodic boundary conditions, so be sure that pbc is set to True.

Restarting from an old calculation
==================================

If the file you specify to Jacapo with the ``nc`` keyword exists, Jacapo will
assume you are attempting to restart an existing calculation. If you do not
want this behavior, turn the flag ``deletenc`` to True in your calculator
definition.

For example, it is possible to continue a geometry optimization with something
like this::

  calc = Jacapo('old.nc', stay_alive=True)
  atoms = calc.get_atoms()
  dyn = QuasiNewton(atoms, logfile='qn.log')
  dyn.run(fmax=0.05)

Note, that the stay_alive flag is not stored in the .nc file and must be set
when the calculator instance is created.  

Atom-projected density of states
================================

To find the atom-projected density of states with Jacapo, first specify the
ados dictionary in your calculator definition, as in::

  calc = Jacapo( ... ,
                ados={'energywindow': (-10., 5.),
                      'energywidth': 0.2,
                      'npoints': 250,
                      'cutoff': 1.0})

After this is established, you can use the get_ados command to get the
desired ADOS data. For example::

  energies, dos = calc.get_ados(atoms=[0],
                                orbitals=['d'],
                                cutoff='short',
                                spin=[0])
