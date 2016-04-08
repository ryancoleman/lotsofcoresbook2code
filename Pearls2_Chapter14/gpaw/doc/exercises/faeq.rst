.. _faeq:

===================================
Frequently asked exercise questions
===================================

.. contents::
    
    
.. _iso:
    
Visualizing iso-surfaces
------------------------

Iso-surfaces from data in a cube-file or a gpw-file can be visualized using
Mayavi from the command line like shown :mod:`here <ase.visualize.mlab>`.
    
It's a good idea to create a short alias like this::
    
    $ alias iso="python -m ase.visualize.mlab -C gpaw"
    
so that you can simply do::
    
    $ iso CO.cube  # plot cube file
    $ iso slab-4.gpw  # plot electron density from gpw-file
    $ iso slab-4.gpw -n 15  # plot wave function from gpw-file
    $ iso -h  # help!

    
.. _xy plot:

Making x-y plots
----------------

If you want to plot an x-y plot from data in a csv-file (comma separated
values), you can use `gnuplot <http://www.gnuplot.info/>`_::

    $ gnuplot
    gnuplot> plot "abc.csv" using 1:2
    
Alternatively, use this little :download:`Python script <xy.py>`:

.. literalinclude:: xy.py

::
    
    $ python <path-to-script>/xy.py abc.csv


Writing 3-d data to cube files
------------------------------

This can be done from Python using the :func:`ase.io.write` function::
    
    from ase.io import write
    write('abc.cube', atoms, data=data)
    
    
Square root
-----------

Square roots are calculated like this: ``2**0.5`` or ``sqrt(2)`` (the
``sqrt`` function must first be imported: ``from math import sqrt`` or
``from numpy import sqrt``).


Integer division
----------------

In python, ``/`` is used for both integer- and float
divisions. Integer division is only performed if both sides of the
operator are integers (you can always force an integer division by
using ``//``)::

  >>> 1 / 3
  0
  >>> 1 / 3.0
  0.33333333333333331

  
Why does changing one variable change another one?
--------------------------------------------------

The = operator in Python is *not* and assignment operator, it is a
*naming* operator:  It makes a new name for (reference to) the object::

  a = [1, 2, 3, 4, 5]  # Create a list
  b = a                # New name for list
  a[2] = 42
  print b              # [1, 2, 42, 4, 5]


  c = 7
  d = c
  c += 42   # d is still 7, we just did
            # c = c + 42
            # creating a new object 49 and
            # giving it the name c


Saving plots
------------

You can save plots made with matplotlib by pressing the floppy-disk
icon in the bottom of the plot, and save as a .png file.

You can save a picture of the atoms from ase-gui by choosing Save, and then
specify a .png file.

You can view .png files in the databar with the command ``eog`` ("eye
of Gnome").
