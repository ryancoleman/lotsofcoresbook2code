===================================
Response and RPA correlation energy
===================================

(Jun Yan, Thomas Olsen)

Excited state properties
========================

RPA correlation energy
======================

The response part of GPAW can be used to calculate ground state correlation energies within the Adiabatic-Connection Dissipation-Fluctuation theorem (ACDF). In Particular, the Random Phase Approximation (RPA) has been implemented and works. 

A major challenge for RPA calculations is the fact that correlation energies can rarely be converged with respect to plane wave cutoff. Instead one has to calculate the correlation energy as a function of cutoff and extrapolate to infinity. When comparing energy differences between isoelectronic systems in the same unit cell one can often neglect this since the error introduced by using a low cutoff tends to cancel out. When comparing energies obtained with different unit cells (lattice constant) or different electronic structure (atomization energies) extrapolaton is essential.

Although the extrapolation is usually rather reliable, it is difficult to obtain an accuracy of a few meV which is needed for van der Waals bonded systems like graphite, BN and MoS2.

Another problem is related to the q=0 contribution. A special method is applied to handle this limit where the coulomb interaction diverges. However, for systems with many degenerate states near the Fermi level (non-noble transition metals) the method breaks down and one cannot trust the q=0 contribtion to the correlation energy. For systems with many kpoints one can simply neglect this term, but for large systems where a sparse kpoint sampling is used, the q=0 term becomes important.

Presently, the research is focussed on going beyond RPA in the ACDF formalism. In particular, including exchange-correlation kernels in the evalution of the interacting response function. More on this later.
