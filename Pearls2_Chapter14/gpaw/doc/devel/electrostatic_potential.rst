.. _electrostatic_potential:

===============================
Note on electrostatic potential
===============================

In the PAW formalism, the electrostatic potential from the
pseudo charge `\tilde{\rho}(\mathbf{r})` is obtained by solving a Poisson
equation:

.. math::

   \nabla^2 \tilde{v}_H(\mathbf{r})=-4\pi\tilde{\rho}(\mathbf{r}).

To get the *real* all-electron electrostatic potential, we need the
all-electron charge density:

.. math::

   \rho(\mathbf{r}) = \tilde{\rho}(\mathbf{r}) +
   \sum_a \Delta\tilde{\rho}^a(\mathbf{r} - \mathbf{R}^a),

where `\Delta\tilde{\rho}^a` is an atomic PAW correction to the pseudo
charge density:

.. math::

   \Delta\tilde{\rho}^a(\mathbf{r}) =
   n_c^a(r) - \tilde{n}_c^a(r) -
   \mathbb{Z}^a\delta(\mathbf{r}) -
   \sum_{\ell=0}^{\ell_{\text{max}}} \sum_{m=-\ell}^\ell
   Q_{\ell m}^a \hat{g}_{\ell m}^a(\mathbf{r}) +
   \sum_{\sigma i_1 i_2} D_{\sigma i_1 i_2}^a
   (\phi_{i_1}^a(\mathbf{r})\phi_{i_2}^a(\mathbf{r}) -
   \tilde{\phi}_{i_1}^a(\mathbf{r})\tilde{\phi}_{i_2}^a(\mathbf{r})).

See :ref:`here <density>` for details.

So, the all-electron potential is:

.. math::

   v_H(\mathbf{r}) = \tilde{v}_H(\mathbf{r}) +
   \sum_a \Delta\tilde{v}_H^a(\mathbf{r} - \mathbf{R}^a)

and 

.. math::

   \Delta\tilde{v}_H^a(\mathbf{r}) =
   \int d\mathbf{r}'
   \frac{\Delta\tilde{\rho}^a(\mathbf{r}')}
   {|\mathbf{r}-\mathbf{r}'|}.

Notice that the `Q_{\ell m}^a` have been choosen so that all multipole
moments of `\Delta\tilde{\rho}^a` are zero and therefore, the
potential from these correction charges (`\Delta\tilde{v}_H^a`) will
be non-zero only inside the atomic augmentation spheres.

The :meth:`~gpaw.aseinterface.GPAW.get_electrostatic_corrections`
method will return an array of integrated corrections:

.. math::

   \int d\mathbf{r} \Delta\tilde{v}_H^a(\mathbf{r})

in units of eV Å\ :sup:`3`.

