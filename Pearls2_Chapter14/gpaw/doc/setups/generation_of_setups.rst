.. _generation_of_setups:

================
Setup generation
================

The generation of setups, starts from a spin-paired atomic
all-electron calculation with spherical symmetry.


All-electron calculations for spherical atoms
=============================================

This is done by the :epydoc:`gpaw.atom.all_electron.AllElectron` class
in :epydoc:`gpaw.atom.all_electron`.  The all-electron wave functions
are defined as:

.. math::

  \phi_{n\ell m}(\mathbf{r}) =
  \frac{u_{n\ell}(r)}{r} Y_{\ell m}(\hat{\mathbf{r}}),

The `u_{n\ell}(r)` functions are stored in an attribute ``u_j`` of the
:epydoc:`gpaw.atom.all_electron.AllElectron` object.  The ``u_j``
member/attribute is an ``ndarray`` with shape ``(nj, N)``, where
``nj`` is the number of states (1s, 2s, 2p, ...) and ``N`` is the
number of radial grid points.

.. tip::

  All-electron calculations can be done with the ``gpaw-setup``
  program like this::

    $ gpaw-setup -a Cu

  Try ``gpaw-setup -h`` for more options.


Generation of setups
====================

The following parameters define a setup:

=================  =======================  =================
name               description              example
=================  =======================  =================
``core``           Froze core               ``'[Ne]'``
``rcut``           Cutoff radius/radii for  ``1.9``
                   projector functions
``extra``          Extra non-bound	    ``{0: [0.5]}``
                   projectors
``vbar``           Zero-potential	    ``('poly', 1.7)``
``filter``         Fourier-filtering	    ``(0.4, 1.75)``
                   parameters
``rcutcomp``	   Cutoff radius for	    ``1.8``
                   compensation charges
=================  =======================  =================

The default (LDA) sodium setup can be generated with the command ``gpaw-setup Na``,
which will use default parameters from the file
``gpaw/atom/generator.py``.
See :ref:`manual_xc` for other functionals.


.. _using_your_own_setups:

Using your own setups
=====================

The setups you generate must be placed in a directory which is included in
the environment variable :envvar:`GPAW_SETUP_PATH` in order for GPAW to
find them. If you want to use the setups in your local directory, add the
following lines to the beginning of your Python script::

    from gpaw import setup_paths
    setup_paths.insert(0, '.')

You can also override the environment variable :envvar:`GPAW_SETUP_PATH` so
that it lists the local directory first and the regular entries afterwards.

If you use bash, :envvar:`GPAW_SETUP_PATH` can be temporarily modified
while you run GPAW with the single command::

    GPAW_SETUP_PATH=.:$GPAW_SETUP_PATH gpaw-python script.py

or if you are using csh or tcsh, you have to first run ``setenv`` and then 
GPAW::

    setenv GPAW_SETUP_PATH .:$GPAW_SETUP_PATH&& gpaw-python script.py

