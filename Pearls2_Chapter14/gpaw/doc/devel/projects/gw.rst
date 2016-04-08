GW approximation
================

:Who:
    Falco

Single-shot GW calculations are implemented and in the trunk and can be used by everybody.
I've started writing documentations and tutorials, but I'm still working on some details.

Recent features are:

* non-linear frequency grid
* Coulomb cutoff for low-dimensional systems
* choice of two different methods for calculating the self energy (Hilbert transform False or True)
* customized parallelization over k-points and frequencies
* read Exact Exchange contributions from file
* spin-polarized GW
* support of grid, LCAO and planewave mode

I've also implemented static COHSEX, but it's not yet in the repository, since I still haven't found a smart way to combine COHSEX 
and GW within the same piece of code.

Right now, I'm preparing and testing self-consistent GW + COHSEX calculations by calculating off-diagonal matrix elements of the 
self energy and Kohn-Sham exchange-correlation contributions and diagonalizing the full Hamiltonian.
