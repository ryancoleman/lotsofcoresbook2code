.. _df_theory:

========================================================
Linear dielectric response of an extended system: theory
========================================================

Introduction
============

The DF (dielectric function) object calculates optical and dielectric properties of extended systems. It computes the linear response function of an interacting many-body system from its ground state electronic structure, which are obtained from GPAW in the real-space grids scheme. 
The frequency and wave-vector dependent density response function are calculated within Time-Dependent Density-Functional Theory formalism using the projector-augmented wave method. Random phase approximation and adiabatic local density approximation are used for exchange-correlation functional. Local field, which indicate the inhomogeneity of the system, are included by constructing the density response matrix in reciprocal space. Fast Fourier Transform (FFT) are utilized to transform between real and reciprocal spaces. 


Refer to :ref:`df_tutorial` for getting started with examples. 



Non-interacting density response function
=========================================

The Non-interacting density response funtion in real space is written as, 

.. math::

  \chi^0(\mathbf{r}, \mathbf{r}^{\prime}, \omega) = \sum_{\mathbf{k}, \mathbf{q}}^{\mathrm{BZ}} \sum_{n, n^{\prime}}
  \frac{f_{n\mathbf{k}}-f_{n^{\prime} \mathbf{k} + \mathbf{q}}}{\omega + \epsilon_{n\mathbf{k}} - \epsilon_{n^{\prime} \mathbf{k} + \mathbf{q} } + i\eta} 
  \psi_{n\mathbf{k}}^{\ast}(\mathbf{r}) \psi_{n^{\prime} \mathbf{k} + \mathbf{q} }(\mathbf{r}) \psi_{n\mathbf{k}}(\mathbf{r}^{\prime}) \psi^{\ast}_{n^{\prime} \mathbf{k} + \mathbf{q} }(\mathbf{r}^{\prime}), 
 
where `\epsilon_{n \mathbf{k}}` 
and `\psi_{n \mathbf{k}}(\mathbf{r})` are the eigenvalue and eigen wavefunction, which 
is normalized to 1 in the crystal volume `\Omega (= \Omega_{\mathrm{cell}} N_k)`.
 
The sum of occupation `f_{n \mathbf{k}}` should be the total number of electrons in the crystal,  
satisfing 

.. math::

  \sum_{n \mathbf{k}} f_{n \mathbf{k}}= N_k N 
 
where `N` is the number of electrons
in a single unit cell and `N_k` is the number of unit cells (kpoints). 


For translation invariant systems,  `\chi^0` can be expanded in planewave basis as

.. math::

  \chi^0(\mathbf{r}, \mathbf{r}^{\prime},  \omega) = \frac{1}{\Omega} 
  \sum_{\mathbf{q}}^{\mathrm{BZ}} \sum_{\mathbf{G} \mathbf{G}^{\prime}}
  e^{i(\mathbf{q} + \mathbf{G}) \cdot \mathbf{r}} \chi^0_{\mathbf{G} \mathbf{G}^{\prime}}(\mathbf{q}, \omega) 
  e^{-i(\mathbf{q} + \mathbf{G}^{\prime}) \cdot \mathbf{r}^{\prime}}


where `\mathbf q` stands for the Bloch vector of the incident wave and `\mathbf G (\mathbf G^{\prime})`
are reciprocal lattice vectors.

The Fourier coefficients, derived by Adler  \ [#Adler]_ and Wiser  \ [#Wiser]_, are written as

.. math::

  \chi^0_{\mathbf{G} \mathbf{G}^{\prime}}(\mathbf{q}, \omega) = \frac{1}{\Omega} 
 \sum_{\mathbf{k}}^{\mathrm{BZ}} \sum_{n, n^{\prime}}
 \frac{f_{n\mathbf{k}}-f_{n^{\prime} \mathbf{k} + \mathbf{q} }}{\omega + \epsilon_{n\mathbf{k}} - \epsilon_{n^{\prime} \mathbf{k} + \mathbf{q} } + i\eta} 
  \langle \psi_{n \mathbf{k}} | e^{-i(\mathbf{q} + \mathbf{G}) \cdot \mathbf{r}} | \psi_{n^{\prime} \mathbf{k} + \mathbf{q} } \rangle_{\Omega_{\mathrm{cell}}} 
  \langle \psi_{n\mathbf{k}} | e^{i(\mathbf{q} + \mathbf{G}^{\prime}) \cdot \mathbf{r}^{\prime}} | \psi_{n^{\prime} \mathbf{k} + \mathbf{q} } \rangle_{\Omega_{\mathrm{cell}}} 


Interacting density response function
=====================================

The full interacting density response function is obtained by solving 
Dyson's equation:

.. math::

  \chi(\mathbf r, \mathbf{r^{\prime}}, \omega) = \chi_0(\mathbf r,  \mathbf{r^{\prime}}, \omega)
  \iint_{\Omega} d\mathbf{r}_1 d\mathbf{r}_2 \ \chi_0(\mathbf r, \mathbf{r}_1, \omega) 
  K(\mathbf{r}_1, \mathbf{r}_2) \chi(\mathbf{r}_2,  \mathbf{r^{\prime}} ,\omega),

where the kernel is the summation of coulomb and exchange-correlation interaction

.. math::

  K(\mathbf{r}_1, \mathbf{r}_2) = \frac{1}{|\mathbf{r}_1 -\mathbf{r}_2|} 
   + \frac{\partial V_{xc}[n]}{\partial n}.  


By Fourier transform, the Dyson's equation in reciprocal space becomes 

.. math::
 
  \chi_{\mathbf G \mathbf G^{\prime}}(\mathbf q, \omega)  
  = \chi^0_{\mathbf G \mathbf G^{\prime}}(\mathbf q, \omega) 
  + \sum_{\mathbf G_1 \mathbf G_2} \chi^0_{\mathbf G \mathbf G_1}(\mathbf q \omega) K_{\mathbf G_1 \mathbf G_2}(\mathbf q)
  \chi_{\mathbf G_2 \mathbf G^{\prime}}(\mathbf q, \omega). 


Note, the coulomb kernel becomes diagonal in reciprocal space

.. math::

   K^{\mathrm{Coulomb}}_{\mathbf G_1 \mathbf G_2}(\mathbf q) = 
   \frac{4\pi}{|\mathbf q+\mathbf G_1|^2} \delta_{\mathbf G_1 \mathbf G_2}


The exchange-correlation (xc) is obtained using adiabatic local density approximation (ALDA):

.. math::

   K^{xc.ALDA}_{\mathbf G_1 \mathbf G_2}(\mathbf q) = 
   \frac{1}{\Omega} \int d\mathbf{r} f_{xc}[n(\mathbf{r})] e^{-i(\mathbf{G}_1-\mathbf{G}_2)\cdot \mathbf{r}}

with 

.. math::

   f_{xc}[n(\mathbf{r})] = \left. \frac{\partial^2 E_{xc}[n]}{\partial n^2} \right|_{n_0(\mathbf{r})}. 


.. _macroscopic_dielectric_function:

Dielectric function and its relation to spectra
===============================================

The dielectric matrix is related to the density response matrix by

.. math::

  \epsilon^{-1}_{\mathbf G \mathbf G^{\prime}}(\mathbf q, \omega) 
  = \delta_{\mathbf G \mathbf G^{\prime}} + \frac{4\pi}{|\mathbf q + \mathbf G|^2} 
  \chi_{\mathbf G \mathbf G^{\prime}}(\mathbf q, \omega)

Within RPA approximation, the dielectric matrix can also be written as

.. math::

  \epsilon^{\mathrm{RPA}}_{\mathbf G \mathbf G^{\prime}}(\mathbf q, \omega)
  = \delta_{\mathbf G \mathbf G^{\prime}} - \frac{4\pi}{|\mathbf q + \mathbf G|^2} 
  \chi^0_{\mathbf G \mathbf G^{\prime}}(\mathbf q, \omega)

The macroscopic dielectric function is defined by

.. math::

  \epsilon_M(\mathbf q, \omega) = \frac{1}{\epsilon^{-1}_{00}(\mathbf q, \omega)}

Optical absorption spectrum is obtained through

.. math::

  \mathrm{ABS} = \mathrm{Im} \epsilon_M(\mathbf q \rightarrow 0, \omega)

Electron energy loss spectrum is 

.. math::

  \mathrm{EELS} = -\mathrm{Im}\frac{1}{\epsilon_M(\mathbf q, \omega)}


The f-sum rule
==============

The scalar dielectric function is related to the 
dielectric tensor by

.. math::

  \epsilon_M(\mathbf q, \omega) = \mathrm{lim}_{\mathbf q \rightarrow 0} 
  \ \hat{q}_{\alpha} \epsilon_{\alpha \beta}(\mathbf q, \omega)  
  \hat{q}_{\beta},

and the dielectric tensor  `\epsilon_{\alpha \beta}(\omega)` satify the "f-sum rule"

.. math::

  \int_0^{\infty}  d\omega \  \omega \ \mathrm{Im} \epsilon_{\alpha \beta}(\omega) 
   = \frac{2\pi^2N}{\Omega_{\mathrm{cell}}} \delta_{\alpha \beta}


where  `N` is the number of electrons in the unit cell and `\frac{N}{\Omega_{\mathrm{cell}}}`
is the electron density.


Optical limit (q -> 0)
======================

In the above sections we have derived the longitudianl dielectric function `\epsilon(\mathbf q, \omega)`. 
For external perturbation by a tranverse  electro-magnetic field, the full dielectric tensor should be 
calculated. However, in the long-wavelength limit, which is the case for light absorption, 
the dielectric tensor can be recovered by scalar or longitudinal dielectric function considering
different direction of `\hat{\mathbf q}`. 

Although  `\mathbf q` is close to zero, 
we can't use the approximation `\mathbf q = 0`
because the Coulomb kernel (`\frac{4\pi}{|\mathbf q + \mathbf G|^2}`) diverges at  `\mathbf q = \mathbf G = 0`. 
In this section we will focus on 
evaluating   `\chi_{\mathbf G \mathbf G^{\prime}}^0(\mathbf q, \omega)`
in the limit of `\mathbf q \rightarrow 0` and `\mathbf G = 0`   \ [#Louie]_. 


The dipole transition matrix  `\langle \psi_{n \mathbf k} | 
e^{-i (\mathbf q + \mathbf G) \cdot \mathbf r} | \psi_{n^{\prime} \mathbf k + \mathbf q} \rangle`
with  `\mathbf G  = 0` becomes

.. math::

  \langle \psi_{n \mathbf k} | 
  e^{-i (\mathbf q + \mathbf G) \cdot \mathbf r} | \psi_{n^{\prime} \mathbf k + \mathbf q} \rangle
  =  \langle u_{n \mathbf k} | u_{n^{\prime} \mathbf k + \mathbf q} \rangle
 

Note, `\psi_{n \mathbf k}` is all-electron wavefunction with band index `n` 
at kpoint `\mathbf k` , and  `u_{n \mathbf k}`  is 
the periodic part of the Bloch wave written as 
`\psi_{n \mathbf k}(\mathbf r) = u_{n \mathbf k}(\mathbf r) e^{i \mathbf k \cdot \mathbf r}`. 


Employing second order perturbation theory, `u_{n^{\prime} \mathbf k + \mathbf q}` 
can be expanded in terms of other orbitals written as

.. math::

  | u_{n^{\prime} \mathbf k + \mathbf q} \rangle
   =  | u_{n^{\prime} \mathbf k } \rangle
   +   \sum_{m \neq n^{\prime}} 
   \frac{ \langle u_{m \mathbf k} | \tilde V | u_{n^{\prime} \mathbf k} \rangle }{\epsilon_{n^{\prime} \mathbf k} - \epsilon_{m \mathbf k} } | u_{m \mathbf k} \rangle
  

where the perturbation `\tilde V` is obtained in the following through k.p perturbation theory. 

The k.p Hamiltonian is expressed as

.. math::
 
  H(\mathbf k) u_{n \mathbf k}(\mathbf r) = \left[ -\frac{\hbar^2}{2m}(\nabla + i\mathbf k)^2 + V(\mathbf r) \right] u_{n \mathbf k}(\mathbf r)
  = \epsilon_{n \mathbf k} u_{n \mathbf k}(\mathbf r),

where `V(\mathbf r)` is the periodic crystal potential. 

The perturbation Hamiltonian `\tilde V` is calculated by (atomic unit):

.. math::

 \tilde V = H(\mathbf k + \mathbf q) - H(\mathbf k) = -i\mathbf q \cdot (\nabla + i \mathbf k)


Substitute  `\tilde V` into the expression of  `| u_{n^{\prime} \mathbf k + \mathbf q} \rangle`, 
multiply `\langle u_{n \mathbf k} |` to the left, 
and apply the orthonormalized condition for the all-electron wavefunction 
`\langle u_{n \mathbf k} | u_{m \mathbf k} \rangle = \delta_{nm}`, we get

.. math::

    \langle \psi_{n \mathbf k} | 
  e^{-i (\mathbf q + \mathbf G) \cdot \mathbf r} | \psi_{n^{\prime} \mathbf k + \mathbf q} \rangle_{\mathbf q \rightarrow 0, \mathbf G=0}
  = -i \mathbf q \cdot \frac{ \langle u_{n \mathbf k} | \nabla + i \mathbf k |u_{n^{\prime} \mathbf k} \rangle }{\epsilon_{n^{\prime} \mathbf k} - \epsilon_{n \mathbf k}} 
  =  -i \mathbf q \cdot \frac{ \langle \psi_{n \mathbf k} | \nabla |\psi_{n^{\prime} \mathbf k} \rangle }{\epsilon_{n^{\prime} \mathbf k} - \epsilon_{n \mathbf k}} 



Hilbert Transform
=================


The non-interaction density response function 
`\chi^0_{\mathbf{G} \mathbf{G}^{\prime}}(\mathbf{q}, \omega)`  can be calculated through 
hilbert transform, written as

.. math::

   \chi^0_{\mathbf{G} \mathbf{G}^{\prime}}(\mathbf{q}, \omega)
   = \int_{-\infty}^{\infty} d\omega^{\prime}
     \frac{A_{\mathbf{G} \mathbf{G}^{\prime}}(\mathbf{q}, \omega^{\prime})}
     {\omega - \omega^{\prime}+ i\eta} 


where spectral function `A_{\mathbf{G} \mathbf{G}^{\prime}}(\mathbf{q}, \omega^{\prime})` 
is defined as

.. math::

   A_{\mathbf{G} \mathbf{G}^{\prime}}(\mathbf{q}, \omega^{\prime})
   =  \frac{1}{\Omega} 
      \sum_{\mathbf{k}}^{\mathrm{BZ}} \sum_{n, n^{\prime}}
      ( f_{n\mathbf{k}}-f_{n^{\prime} \mathbf{k} + \mathbf{q}} )
       \langle \psi_{n \mathbf{k}} | e^{-i(\mathbf{q} + \mathbf{G}) \cdot \mathbf{r}} | \psi_{n^{\prime} \mathbf{k} + \mathbf{q} } \rangle_{\Omega_{\mathrm{cell}}} 
  \langle \psi_{n\mathbf{k}} | e^{i(\mathbf{q} + \mathbf{G}^{\prime}) \cdot \mathbf{r}^{\prime}} | \psi_{n^{\prime} \mathbf{k} + \mathbf{q} } \rangle_{\Omega_{\mathrm{cell}}} 
       \times \delta( \omega^{\prime} + \epsilon_{n\mathbf{k}} - \epsilon_{n^{\prime} \mathbf{k} + \mathbf{q} }  )

Note that the integration above requires both positive and negative frequencies. 
In the following derivation, the  intergration will be reduced to only half of the frequency domain. 

In the system that possesses the time-reversal symmetry, the bloch states have the following properties

.. math::
    
   \epsilon_{n, -\mathbf{k}} = \epsilon_{n, \mathbf{k}}

   f_{n, -\mathbf{k}} = f_{n, \mathbf{k}}
   
   \psi_{n, -\mathbf{k}}(\mathbf{r}) = \psi^{\ast}_{n, \mathbf{k}}(\mathbf{r})


Change the index in `A_{\mathbf{G} \mathbf{G}^{\prime}}(\mathbf{q}, \omega^{\prime})`
as 

.. math::

   n, \mathbf{k}   \rightarrow n^{\prime}, -\mathbf{k}-\mathbf{q}

   n^{\prime}, \mathbf{k}+\mathbf{q} \rightarrow  n, -\mathbf{k}  

and employing the time-reversal symmetry, one can get

.. math::

   A_{\mathbf{G} \mathbf{G}^{\prime}}(\mathbf{q}, \omega^{\prime})
   =  - A_{\mathbf{G} \mathbf{G}^{\prime}}(\mathbf{q}, -\omega^{\prime})   

Substitute it to the integration in the beginning of this section, one get

.. math::

     \chi^0_{\mathbf{G} \mathbf{G}^{\prime}}(\mathbf{q}, \omega)
     = \int_0^{\infty} d\omega^{\prime} 
       \frac{ A_{\mathbf{G} \mathbf{G}^{\prime}}(\mathbf{q}, \omega^{\prime})}{\omega-\omega^{\prime}+i\eta}
       + \int_{-\infty}^{0}  d\omega^{\prime} 
       \frac{ A_{\mathbf{G} \mathbf{G}^{\prime}}(\mathbf{q}, \omega^{\prime})}{\omega-\omega^{\prime}+i\eta} 
     =  \int_0^{\infty} d\omega^{\prime} 
          \left[ \frac{1}{ \omega-\omega^{\prime}+i\eta } - 
                 \frac{1}{ \omega+\omega^{\prime}+i\eta }\right]
            A_{\mathbf{G} \mathbf{G}^{\prime}}(\mathbf{q}, \omega^{\prime})
         
Applying the hilbert transform can make the calculations of `\chi^0_{\mathbf{G} \mathbf{G}^{\prime}}(\mathbf{q}, \omega)` `Nw / 2` times faster, where `Nw` is the number of frequency points used. 

For the delta function, we use either a triangular function, which is described in  \ [#DeltaFunc]_ and
is normalized to 1 or a gaussian function, which is in principle normalized but in fact not due to  numerically finite frequency 
points used. We tried both and it turns out that the spectrum does not sensitively depend on the function applied.  



PAW terms
=========

The PAW terms comes in when calculating the dipole transition matrix 

.. math::

   \langle \psi_{n \mathbf k} | 
   e^{-i (\mathbf q + \mathbf G) \cdot \mathbf r} | \psi_{n^{\prime} \mathbf k + \mathbf q} \rangle 
   = \langle \tilde{\psi}_{n \mathbf k} | 
   e^{-i (\mathbf q + \mathbf G) \cdot \mathbf r} | \tilde{\psi}_{n^{\prime} \mathbf k + \mathbf q} \rangle 
   + \sum_{a,ij} 
   \langle  \tilde{\psi}_{n \mathbf k} | \tilde{p}_i^a  \rangle^{\ast}
   \langle \tilde{\psi}_{n^{\prime} \mathbf k + \mathbf q} | \tilde{p}_j^a   \rangle
   \left[ \langle \phi_i^a | e^{-i(\mathbf{q} + \mathbf{G}) \cdot \mathbf{r}} | \phi_j^a \rangle
         - \langle \tilde{\phi}_i^a | e^{-i(\mathbf{q} + \mathbf{G}) \cdot \mathbf{r}} | \tilde{\phi}_j^a \rangle
   \right]


We calculate the last term in the above equation by expanding the planewave in such a way

.. math::

   e^{i \mathbf{k} \cdot \mathbf{r}} = 4 \pi \sum_{lm} i^l j_l(kr) Y_{lm}(\hat{\mathbf{r}})  Y_{lm}(\hat{\mathbf{k}}) 

where `j_l` is spherical bessel function
and write (for simplicity, define `\mathbf{k} = \mathbf{q} + \mathbf{G}`)

.. math::

    \langle \phi_i^a | e^{-i \mathbf{k} \cdot \mathbf{r}} | \phi_j^a \rangle
         - \langle \tilde{\phi}_i^a | e^{-i \mathbf{k} \cdot \mathbf{r}} | \tilde{\phi}_j^a \rangle
    = 4 \pi e^{-i \mathbf{k} \cdot \mathbf{R}_a}  \sum_{lm} (-i)^l  Y_{lm}(\hat{\mathbf{k}}) 
        \int dr \ r^2  j_l(kr) \left[ \phi^{a}_{n_1 l_1}(r)  \phi^{a}_{n_2 l_2}(r) 
                                     -  \tilde{\phi}^{a}_{n_1 l_1}(r)  \tilde{\phi}^{a}_{n_2 l_2}(r) \right] 
        \int d\Omega \  Y_{lm} Y_{l_1 m_1} Y_{l_2 m_2}      

where `\mathbf{R}_a` are the positions of atoms in the unit cell. 


For optical limit calculation, the dipole matrix related is 

.. math::

     \langle \psi_{n \mathbf{k}} | \nabla | \psi_{n^{\prime} \mathbf{k}} \rangle
     = \langle \tilde{\psi}_{n \mathbf{k}} | \nabla | \tilde{\psi}_{n^{\prime} \mathbf{k}} \rangle
       +  \sum_{a,ij} 
   \langle  \tilde{\psi}_{n \mathbf k} | \tilde{p}_i^a  \rangle^{\ast}
   \langle \tilde{\psi}_{n^{\prime} \mathbf k} | \tilde{p}_j^a   \rangle
   \left[ \langle \phi_i^a | \nabla_{\mathbf{r}} | \phi_j^a \rangle
         - \langle \tilde{\phi}_i^a | \nabla_{\mathbf{r}} | \tilde{\phi}_j^a \rangle
   \right]

Refer to :ref:`setup_matrix_elements_nabla`
for calculation of  `\langle \phi_i^a | \nabla_{\mathbf{r}} | \phi_j^a \rangle - \langle \tilde{\phi}_i^a | \nabla_{\mathbf{r}} | \tilde{\phi}_j^a \rangle`

.. [#Adler] S. L. Adler,
            Quantum theory of the dielectric constant in real solids,
            *Phys. Rev.* **126**, 413 (1962)

.. [#Wiser] N. Wiser, 
            Dielectric constant with local field effects included, 
            *Phys. Rev.* **129**, 62 (1963).

.. [#Louie] M. S. Hybertsen and S. G. Louie, 
            Ab initio static dielectric matrices from the density-functional
            approach. I. Formulation and application to semiconductors and 
            insulators, 
            *Phys. Rev. B* **35**, 5585 (1987).

.. [#DeltaFunc] M. Shishkin and G. Kresse, 
                Implementation and performance of the frequency-dependent GW
                method within the PAW framework, 
                *Phys. Rev. B* **74**, 035101 (2006).
