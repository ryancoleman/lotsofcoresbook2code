.. _transport:

=========
Transport
=========

The Transport object in GPAW has been written as a calculator, 
different for other calculators, it supports open boundary
condition, it needs the information of the electrodes and 
scattering region, and can calculate the density, hamiltonian, 
total energy, forces and current as well.

Model Picture
-------------

::

        ______ ______ ____________ ______ ______ 
       |      |      |     __     |      |      |
   ....|______|______|___ /  \____|______|______|....
       |      |      |    \__/    |      |      |
       |______|______|____________|______|______|
       |      |                          |      |
       | Lead |---->   Scattering   <----| Lead |
       |      |          region          |      |



How to describe an open system
------------------------------

The total open system includes some semi-infinite electrodes
or leads connected to the scattering region we focus on. 
When it is far away enough from the scattering region, 
the lead parts are quite close to the periodical case.
So here we divide the total system into two parts: one is 
the lead parts, there is neither reconstruct nor charge 
transfer, all the information of it can be got from the 
periodical calculation, the other is the scattering region, 
we need a self-consistent procedure to get the properites
here.

Leads
-----

The influence of leads to the scattering region is absorbed
in a item named surface Green's function

.. math::

  g(E) = (E*S_l - H_l)^{-1}

`S_l` and `H_l` are the overlap and hamiltonian matrices of
leads respectively, and since they are inifinite, we need to
do some handling to get it.

We induce the concept of priciple layer, the unit cell when 
solving the surface Green's function. We assum interaction 
only exsits between two adjacent principle layers.
That means the Hamiltonian matrix is tridiagonal by the size
of a principle layer, which is necessary to get the surface
Green's function. 

The selfenergy of lead can be calculated like this

.. math::

  \Sigma _l(E) = \tau _l g(E) \tau _l^{+}

.. math::

  \tau _l = E * S_{lc} - H_{lc}

`S_{lc}` and `H_{lc}` are the coupling overlap and Hamiltonian
matrices. The detailed process can be accessed in Phys. Rev. B 28 4397.

Scattering region
-----------------

For the reason mentioned above we need to choose a relatively
big scattering region, some atoms in the leads, or generally
speaking, at least one principle layer should be included in
the scattering region.

The retarded Green's function of scattering region is written as

.. math::

  G^r(E) = ((E + i\eta)*S - H -\Sigma) ^ {-1}

`\Sigma` is the sum of the leads selfeneries.

The lesser Green's function is from the Keldysh formalism

.. math::

  G^<(E) = G^r\Sigma^<G^a

With these two Green's functions, we can get the electron
density in non-equilirbium case which will be introduced
later.

Keywords for Setup
------------------

For leads(necessary):

=================  =========    ============================
keyword            type         description
=================  =========    ============================
``pl_atoms``       list         :ref:`manual_pl_atoms`  
``pl_cells``       list         :ref:`manual_pl_cells`
``pl_kpts``        list         :ref:`manual_pl_kpts`
=================  =========    ============================

Keywords for gpaw is inherited, which is used to descibe
scattering region. Something special for Transport is:

* ``mode`` should be ``'lcao'`` always.

* if use fixed_boundary_condition, ``poissonsolver`` should be set
  like ``PoissonSolver(nn=x)``.

* ``usesymm`` does not act, Transport set a value for it automatically.
 
Other keywords (used frequently):


=========================  =====  =============  ===================================
keyword                    type   default value  description
=========================  =====  =============  ===================================
``non_sc``                 bool   False          :ref:`manual_non_sc`
``plot_eta``               float  1e-4           :ref:`manual_plot_eta`
``plot_energy_range``      list   [-5,5]         :ref:`manual_plot_energy_range`
``plot_energy_point_num``  int    201            :ref:`manual_plot_energy_point_num` 
``analysis_mode``          bool   False          :ref:`manual_analysis_mode`
=========================  =====  =============  ===================================

Usage:

Calculate transmission using hamiltonian from normal DFT

.. literalinclude:: transport_non_sc.py

Get an iv curve using NEGF:

.. literalinclude:: transport.py

A spin transport example (anti-parallel junction):

.. literalinclude:: spin_transport.py
 
Calculate transmission and DOS based on separate DFT results for electrodes and scattering region(for example
one may want to carry out a PBE+U calculation for the scattering region):

.. literalinclude:: transport_from_dft.py

Analysis:

>>> from gpaw.transport.analysor import Transport_Plotter
>>> plotter = Transport_Plotter()
>>> data = plotter.get_info(XXX, 0, 0) #information string, bias_step, ion_step 

Transport_Plotter now just get the data, users need to plot the data themselves.
XXX can be one in the list ['tc', 'dos', 'force', 'lead_fermi', 'bias', 'gate', 'nt', 'vt'].
The analysis functionality only works after a transport calculation
is done successfully and the directory analysis_data and some files in it are generated. 

Analysis Package:

Some small scripts are provided in the analysis_scripts directory to analyze the result. 
One can include them in python path, and more convenient way is to append some 
alias lines in .bashrc or cshrc. 

.. literalinclude:: alias_lines

The syntax is 
Transmission: 

>>> tc X Y

Density of States:

>>> dos X Y

X is the bias step index, Y is the ionic step index, the default for Y is 0

Average Effective Potential or Pseudo Density(average x, y directions):

>>> vt X Y

Here X support a linking symbol -, i.e. vt 5-0 plot the potential difference
for the bias step 5 and bias step 0, which can be used to judge the 
screening effects. In this plot, x axis is the realspace position in transport
direction, and the potential of one principle layer, which is from the 
electrode calculation, will be attached outside.

>>> nt X Y

Similiar to vt X Y, it plot the pseudo density.

more...

>>> vtx X Y
>>> vty X Y
>>> ntx X Y
>>> nty X Y

These commands make 2-d color plots, the data is average in one direction(x or y).

Partial Density of States:

>>> pdos X Y ZZZ

ZZZ is a description of the partial orbital,
i.e., C=40_S plot the S orbital DOS of a C atom whose index is 40;
or,   C]30[40_P plot the P orbital DOS of some C atoms whose index range from 30 to 40;
or,   C-H_S  plot the S orbital DOS of all C atoms and all H atoms;
or,   C-H_S C-H_P C-H_A plot the S orbital, P orbital, and total DOS of C and H atoms
respectively in one figure.
As a summery, the linking symbol - means plus a element, _ means the orbital type,
],[,= defines the atomic indices, A represent all the atoms or orbtials, so dos X Y is
equalivalent to pdos X Y A_A. This command support multi inputparameters seperated by
space.

IV Characteristic:

>>> iv X

It plot out the IV curve for the fisrt X bias points.

Charge:

>>> charge X Y ZZZ

Usage is similiar to pdos, it plot the charge of the partial orbitals as a function 
of bias.

Force:

>>> force X Y

It plots the force of bias step X and ionic step Y, X here also supports linking symbol -,
i.e., force 1-0 plot the force differnece for bias step 1 and bias step 0.

Optional keywords:

=====================  ===========      =============  ===============================
    keyword              type           default value          description
=====================  ===========      =============  ===============================
``bias``                 list               [0, 0]     :ref:`manual_bias`  
``gate``                 float                0        :ref:`manual_gate`
``fixed_boundary``       bool                True      :ref:`manual_fixed_boundary`
``lead_restart``         bool               False      :ref:`manual_lead_restart`
``scat_restart``         bool               False      :ref:`manual_scat_restart`
``cal_loc``              bool               False      :ref:`manual_cal_loc`
``recal_path``           bool               False      :ref:`manual_recal_path`
``use_buffer``           bool               False      :ref:`manual_use_buffer`
``buffer_atoms``         list                []        :ref:`manual_buffer_atoms`
``use_qzk_boundary``     bool               False      :ref:`manual_use_qzk_boundary`
``identical_leads``      bool               False      :ref:`manual_identical_leads`
``normalize_density``    bool               True       :ref:`manual_normalize_density`
``alpha``                float               0.0       :ref:`manual_alpha`
``gate_fun``           numpy array          None       :ref:`manual_gate_fun`
=====================  ===========      =============  ===============================
 
.. _manual_pl_atoms:


Principle Layer Atoms
---------------------

``pl_atoms`` is the index of lead atoms, whose length is the 
number of leads. For example, [[0,1,2,3],[7,8,9,10]] means there
are two leads, [0,1,2,3] is the principle layer of the first
lead and [7,8,9,10] for the second. The sequence is arbitary.

.. _manual_pl_cells:

Principle Layer Cells
---------------------

``pl_cells`` is a list of leads' cells, also has the same length
with the leads number. [[10., 10., 30], [10., 10., 30.]] for example.
For two-probe system, the lead cell should have the same size with that of 
scattering region in x, y directions.

.. _manual_pl_kpts:

Principle Layer K-points
------------------------

``pl_kpts`` is k-points sampling for leads, it is a 1*3 int sequence.
We just let all the leads have the same K number. Attention here that
the k number in the transport direction should bigger than 3, 
in principle we should have enough k points in this direction, an
experenced rule is nK * L(Å) ~ 50. L is the length of unit cell
in this direction. Note that pl_kpts should match the K sampling of 
scattering region kpts in the x and y direction and in the parallel
case, the local K sampling should match as well. So a safe way is to use
a prime number for the pl_kpts in the z axis.

.. _manual_bias:

Bias
----

``bias`` is a list of bias value for all the leads. For example,
[1.0, -1.0] means we have two leads with the bias shift 1.0 and -1.0
respectively.

.. _manual_gate:

Gate
----

``gate`` is a float number that should only make some sense with 
the fixed boundary condition. The atoms on what a constant gate
is applied is the total scattering region minus the lead's
principle layers.

.. _manual_fixed_boundary:

Fixed Boundary Condition
------------------------

``fixed_boundary`` is a bool option. If set True, we solve the
Poisson equation for the scattering region with fixed boundary
condition. It workes when ``pbc`` in the transport direction
for the scattering region is True and ``poissonsolver=PoissonSolver(nn=X)``. 
If set False, Transport object will deal with a 
regular gpaw option which depends on ``pbc``.

.. _manual_lead_restart:

Lead Calculation Restart
------------------------

``lead_restart`` is a bool option decides if restart from
some previous calculation or not. It is for leads especially.

.. _manual_scat_restart:

Scattering Region Calculation Restart
-------------------------------------

``scat_restart`` same with ``lead_restart``, this is just
for scattering region.

.. _manual_cal_loc:

Double Path Integral
--------------------

``cal_loc`` If set True, Transport will adopt a more complicate
schem for the Green's function integral, which has better
precision and costs more cpu time at the same time.

.. _manual_recal_path:

Recalculate Integral Path
-------------------------

``recal_path`` When doing the Green's function integral,
the energy points on the integral path depends on the 
hamiltonian. The default option is we get the path info with
a hamiltonian guess and then fix it, it often works fine,
saving much time for calculation the leads selfenergy, but when
the guess is not good enough, the integral result differs 
from the real value. This keyword will force to refine the 
energy points on path in each SCF iteration, it should be
a option when can not get a convergence result.

.. _manual_use_buffer:

Use Buffer Atoms
----------------

``use_buffer`` Buffer atoms are needed somtime, this part of
atoms are included when calculating hamiltonian, but will be 
neglected when calculating the Green's function, that means
the density of this part is fixed. If set True, you should
provide the information in ``buffer_atoms``.

.. _manual_buffer_atoms:

Buffer Atoms Indices
--------------------

``buffer_atoms`` is the list of buffer atoms index, just like
``pl_atoms``.

.. _manual_edge_atoms:

Edge Atoms
----------

``edge_atoms`` One needs to point which atom
is used to align the energy levels, that means in ground
state calculation, the hamiltonian matrix should have the same 
value in the orbitals of that atom. It is a list includes 
two sub lists. For example, [[0,3],[0,9]] means the atom 0 in
lead1, atom 3 in lead2 are equal to atom 0, atom 9 in scattering 
region respectively.

.. _manual_use_qzk_boundary:

Use Qzk Boundary Condition
--------------------------

``use_qzk_boundary`` is a particular keyword corresponding to
a method introduced in the paper J.Chem.Phys 194710, 127(2007)

.. _manual_identical_leads:

Identical Leads
---------------

When the two electrodes are exactly the same, including cell, atom positions,
and also stacking, set ``identical_leads`` as True can help to save
some time for electrode calculation.

.. _manual_normalize_density:

Normalize Density
-----------------

In normal DFT calculation, the density is always scaled to satisfy the charge
neutrality condition. Because the charge is conserved always by filling up
the molecular levels, the scaling factor is very close to 1. There the scaling
just helps to converge, and does not influece the calculation result. In NEGF
calculation, because the density matrix is obtain by the Green's function integral
to the fermi level, the charge neutrality is not garanteed. Without scaling, there
may be a convergence problem for the three dimensional system. So the code 
is forced to do scaling always for the first several steps. The keyword 
``normalize_density`` decides whether the scaling will be released or not at last.
In principle it should be released, but there is only very tiny diference if
the scaling factor is close to 1 at last.

.. _manual_alpha:

Alpha
-----
If scaling the electron density, sometime you can meet a trap, when the scaling 
factor oscillating far away from 1. If this situation happen, just set 
``alpha`` 0.1~0.3, then when you do scaling, there will be a small amount of
net charge after scaling, which will help the calculatoin jump out the trap and 
not get diverged.

.. _manual_gate_fun:

Gate Function
-------------
Gate Function describes the external gate potential shape in transport direction.
The numbers in this array shoule be 1 in the middle, and decays to 0 on both sides,
for example np.array([0,0,0.1,0.5,1,1,0.5,0.1,0,0]). Interpolation will be done
automatically to fit the function to the whole scaterring region.

.. _manual_non_sc:

NonSelfConsistent
-----------------

When ``non_sc`` is True, the scattering region is calculated by a normal
DFT using a periodical slab. 

.. _manual_analysis_mode:

Analysis Mode
-------------

If analysis_mode is set True, the hamiltonian in the bias data file will be
read and restart the transmission calculation with the k sampling
adjustable, this is often used to project the trammission into dense k sampling.


.. _manual_plot_eta:

Plot Eta
--------

When calculating the Green function for the transmision plot, a small 
imaginary float is added to the real energy to avoid singularity, in principle
eta should be a infinitesimal, while the bigger eta can smear the 
sharp peaks in transmission and dos.

.. _manual_plot_energy_range:

Plot Energy Range
-----------------

The energy scale in the transmission plot.

.. _manual_plot_energy_point_num:

Plot Energy Point Num
---------------------

The number of energy points used in the transmission plot which determines the
density of the plot.



