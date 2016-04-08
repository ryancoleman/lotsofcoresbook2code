.. _lattice_constants:

=========================
Finding lattice constants
=========================

.. seealso::

   * `ASE EOS tutorial
     <https://wiki.fysik.dtu.dk/ase/tutorials/eos/eos.html>`_
   * `ASE Finding lattice constants tutorial
     <https://wiki.fysik.dtu.dk/ase/tutorials/lattice_constant.html>`_

   * `ASE equation of state module
     <https://wiki.fysik.dtu.dk/ase/ase/utils.html#equation-of-state>`_


Fcc Aluminium
=============

Let's try to converge the lattice constant with respect to number of
plane-waves:

.. literalinclude:: al.py
    :lines: 1-19

.. image:: Al_conv_ecut.png

Using a plane-wave cutoff energy of 400 eV, we now check convergence
with respect to number of **k**-points:

.. literalinclude:: al.py
    :lines: 21-

.. image:: Al_conv_k.png

(see also :download:`analysis script <al.agts.py>`).
