.. _tools:
    
Command-line tools
==================

Finding all or some unocupied states
------------------------------------

If you have a gpw-file containing the ground-state density for a plane-wave
calculation, then you can set up the full
`H_{\mathbf{G}\mathbf{G}'}(\mathbf{k})` and
`S_{\mathbf{G}\mathbf{G}'}(\mathbf{k})` matrices in your plane-wave basis and
use direct diagonalization to find all the eigenvalues and eigenstates in one
step.

Usage::
    
    $ python -m gpaw.fulldiag [options] <gpw-file>
    
Options:
    
-h, --help            Show this help message and exit
-n BANDS, --bands=BANDS
                      Number of bands to calculate.  Defaults to all.
-s SCALAPACK, --scalapack=SCALAPACK
                      Number of cores to use for ScaLapack.  Default is one.
-d, --dry-run         Just write out size of matrices.

Typpically, you will want to run this in parallel and distrubute the matrices
using ScaLapack::
    
    $ mpiexec gpaw-python -m gpaw.fulldiag abc.gpw --scalapack=8
