.. module:: ase.db

====================
A database for atoms
====================

ASE has its own database that can be used for storing and retrieving atoms and
associated data in a compact and convenient way.
    
.. note::

    This is work in progress.  Use at your own risk!
    
There are currently three back-ends:

JSON_:
    Simple human-readable text file with a ``.json`` extension.
SQLite3_:
    Self-contained, server-less, zero-configuration database.  Lives in a file
    with a ``.db`` extension.
PostgreSQL_:
    Server based database.
MariaDB_ or MySQL:
    Server based database.

The JSON and SQLite3 back-ends work "out of the box", whereas PostgreSQL
requires a server.

There is a `command-line tool`_ called :program:`ase-db` that can be
used to query and manipulate databases and also a `Python interface`_.

.. _JSON: http://www.json.org/
.. _SQLite3: http://www.sqlite.org/
.. _PostgreSQL: http://www.postgresql.org/
.. _MariaDB: http://mariadb.org/

.. contents::
    

What's in the database?
=======================

Every row in the database contains:

* all the information stored in the :class:`~ase.atoms.Atoms` object
  (positions, atomic numbers, ...)
* calculator name and parameters (if a calculator is present)
* already calculated properties such as energy and forces
  (if a calculator is present)
* key-value pairs (for finding the calculation again)
* an integer ID (unique for each database) starting with 1 and always
  increasing for each new row
* a unique ID which is a 128 bit random number which should be globally
  unique (at least in the lifetime of our universe)
* constraints (if present)
* user-name
* creation and modification time


.. _ase-db command:
    
Command-line tool
=================

The :program:`ase-db` command can be used to query databases and for
manipulating key-value pairs.  Try::
    
    $ ase-db --help
    
Example: Show all rows of SQLite database abc.db:
    
.. literalinclude:: ase-db.out
    
Show all details for a single row:
    
.. literalinclude:: ase-db-long.out

    
Querying
--------

Here are some example query strings:
    
.. list-table::
    :widths: 25 75
    
    * - Cu
      - contains copper
    * - H<3
      - less than 3 hydrogen atoms
    * - Cu,H<3
      - contains copper and has less than 3 hydrogen atoms
    * - v3
      - has 'v3' key
    * - abc=bla-bla
      - has key 'abc' with value 'bla-bla'
    * - v3,abc=bla-bla
      - both of the above
    * - calculator=nwchem
      - calculations done with NWChem
    * - 2.2<bandgap<3.0
      - 'bandgap' key has value between 2.2 and 3.0
    * - natoms>=10
      - 10 or more atoms
    * - formula=H2O
      - Exactly two hydrogens and one oxygen
    * - id=2345
      - specific id
    * - age<1h
      - not older than 1 hour
    * - age>1y
      - older than 1 year

These names are special:

.. list-table::
    :widths: 25 75
    
    * - id
      - integer identifier
    * - natoms
      - number of atoms
    * - formula
      - formula
    * - energy
      - potential energy
    * - charge
      - total charge
    * - magmom
      - total magnetic moment
    * - calculator
      - name of calculator
    * - user
      - who did it
    * - age
      - age of calculation (use s, m, h, d, w, M and y for second, minute,
        hour, day, week, month and year respectively)

        
Integration with other parts of ASE
===================================

ASE's :func:`ase.io.read` function can also read directly from databases:
    
>>> from ase.io import read
>>> a = read('abc.db@42')
>>> a = read('abc.db@id=42')  # same thing
>>> b = read('abc.db@v3,abc=H')

Also the :program:`ase-gui` program can read from databases using the
same syntax.
        

Python Interface
================

.. module:: ase.db.core

First, we :func:`connect` to the database:
    
>>> from ase.db import connect
>>> c = connect('abc.db')

or

>>> import ase.db
>>> c = ase.db.connect('abc.db')

Let's do a calculation for a hydrogen molecule and write some results to a
database:
    
>>> from ase import Atoms
>>> from ase.calculators.emt import EMT
>>> h2 = Atoms('H2', [(0, 0, 0), (0, 0, 0.7)])
>>> h2.calc = EMT()
>>> h2.get_forces()
array([[ 0.        ,  0.        , -9.80290573],
       [ 0.        ,  0.        ,  9.80290573]])

Write a row to the database with a key-value pair (``'relaxed'``, ``False``):
    
>>> c.write(h2, relaxed=False)
1

The :meth:`~Database.write` method returns an integer id.

Do one more calculation and write results:
    
>>> from ase.optimize import BFGS
>>> BFGS(h2).run(fmax=0.01)
BFGS:   0  12:49:25        1.419427       9.8029
BFGS:   1  12:49:25        1.070582       0.0853
BFGS:   2  12:49:25        1.070544       0.0236
BFGS:   3  12:49:25        1.070541       0.0001
>>> c.write(h2, relaxed=True)
2

Loop over selected rows using the :meth:`~Database.select` method:
    
>>> for d in c.select(relaxed=True):
...     print d.forces[0, 2], d.relaxed
...
-9.8029057329 False
-9.2526347333e-05 True

The :meth:`~Database.select` method will generate dictionaries that one can
loop over.  The dictionaries are special in the sense that keys can be
accessed as attributes also (``d.relaxed == d['relaxed']``).

Write the energy of an isolated hydrogen atom to the database:

>>> h = Atoms('H')
>>> h.calc = EMT()
>>> h.get_potential_energy()
3.21
>>> c.write(h)
3

Select a single row with the :meth:`~Database.get` method:
    
>>> d = c.get(relaxed=1, calculator='emt')
>>> for k, v in d.items():
...     print '%-25s: %s' % (k, v)
...
user                     : jensj
key_value_pairs          : {u'relaxed': True}
energy                   : 1.07054126233
positions                : [[ ... ]]
calculator               : emt
relaxed                  : True
calculator_parameters    : {}
cell                     : [[ 1.  0.  0.] [ 0.  1.  0.] [ 0.  0.  1.]]
numbers                  : [1 1]
forces                   : [[ ... ]]
mtime                    : 14.8975933527
pbc                      : [False False False]
data                     : {u'abc': array([1, 2, 3])}
id                       : 2
unique_id                : bd6e07125cf7a46569a3ed361a2edbe8
ctime                    : 14.8975933527

Calculate the atomization energy and :meth:`~Database.update` a row in
the database:
    
>>> e2 = d.energy
>>> e1 = c.get(H=1).energy
>>> ae = 2 * e1 - e2
>>> print(ae)
5.34945873767
>>> id = c.get(relaxed=1).id
>>> c.update(id, atomization_energy=ae)
1

Delete a single row:
    
>>> del c[c.get(relaxed=0).id]

or use the :meth:`~Database.delete` method to delete several rows.


Dictionary representation of rows
---------------------------------

The first 9 keys (from "id" to "positions") are always present --- the rest
may be there:
    
=====================  =================================  ============  ======
key                    description                        datatype      shape
=====================  =================================  ============  ======
id                     Local database id                  int
unique_id              Globally unique hexadecimal id     str
ctime                  Creation time                      float
mtime                  Modification time                  float
user                   User name                          str
numbers                Atomic numbers                     int           (N,)
pbc                    Periodic boundary condition flags  bool          (3,)
cell                   Unit cell                          float         (3, 3)
positions              Atomic positions                   float         (N, 3)
initial_magmoms        Initial atomic magnetic moments    float         (N,)
initial_charges        Initial atomic charges             float         (N,)
masses                 Atomic masses                      float         (N,)
tags                   Tags                               int           (N,)
momenta                Atomic momenta                     float         (N, 3)
constraints            Constraints                        list of dict
energy                 Total energy                       float
forces                 Atomic forces                      float         (N, 3)
stress                 Stress tensor                      float         (6,)
dipole                 Electrical dipole                  float         (3,)
charges                Atomic charges                     float         (N,)
magmom                 Magnetic moment                    float
magmoms                Atomic magnetic moments            float         (N,)
calculator             Calculator name                    str
calculator_parameters  Calculator parameters              dict
=====================  =================================  ============  ======


Extracting Atoms objects from the database
------------------------------------------

If you want an Atoms object insted of a dictionary, you should use the
:meth:`~Database.get_atoms` method:

>>> h2 = c.get_atoms(H=2)

or if you want the original EMT calculator attached:
    
>>> h2 = c.get_atoms(H=2, attach_calculator=True)


Add additional data
-------------------

When you write a row to a database using the :meth:`~Database.write` method,
you can add key-value pairs where the values can be
strings, floating point numbers, integers and booleans:
    
>>> c.write(atoms, functional='LDA', distance=7.2)

More complicated data can be written like this:

>>> c.write(atoms, ..., data={'parents': [7, 34, 14], 'stuff': ...})

and accessed like this:

>>> d = c.get(...)
>>> d.data.parents
[7, 34, 14]


More details
------------

.. autofunction:: connect

.. autoclass:: ase.db.core.Database
    :members:
    :member-order: bysource
