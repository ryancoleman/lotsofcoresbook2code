import numpy as np

from gpaw.utilities.blas import gemm
from gpaw.utilities import unpack


def get_atomic_hamiltonian(name):
    cls = dict(dense=DenseAtomicHamiltonian,
               distributed=DistributedAtomicHamiltonian)[name]
    return cls()


class DenseAtomicHamiltonian:
    name = 'dense'
    description = 'dense with blas'

    def redistribute(self, wfs, dH_asp):
        return dH_asp
    
    def calculate(self, wfs, kpt, dH_asp, H_MM, y):
        Mstart = wfs.ksl.Mstart
        Mstop = wfs.ksl.Mstop
        dtype = wfs.dtype
        for a, P_Mi in kpt.P_aMi.items():
            dH_ii = np.asarray(unpack(dH_asp[a][kpt.s]), dtype)
            dHP_iM = np.zeros((dH_ii.shape[1], P_Mi.shape[0]), dtype)
            # (ATLAS can't handle uninitialized output array)
            gemm(1.0, P_Mi, dH_ii, 0.0, dHP_iM, 'c')
            gemm(y, dHP_iM, P_Mi[Mstart:Mstop], 1.0, H_MM)


class DistributedAtomicHamiltonian:
    name = 'distributed'
    description = 'block-sparse by atoms'
    
    def redistribute(self, wfs, dH_asp):
        def get_empty(a):
            ni = wfs.setups[a].ni
            return np.empty((wfs.ns, ni * (ni + 1) // 2))

        # just distributed over gd comm.  It's not the most aggressive
        # we can manage but we want to make band parallelization
        # a bit easier and it won't really be a problem.  I guess
        #
        # Also: This call is blocking, but we could easily do a
        # non-blocking version as we only need this stuff after
        # doing tons of real-space work.
        return wfs.atom_partition.to_even_distribution(dH_asp,
                                                       get_empty,
                                                       copy=True)
    
    def calculate(self, wfs, kpt, dH_asp, H_MM, y):
        dtype = wfs.dtype
        M_a = wfs.setups.M_a
        nM_a = np.array([setup.nao for setup in wfs.setups])
        Mstart = wfs.ksl.Mstart
        Mstop = wfs.ksl.Mstop

        # Now calculate basis-projector-basis overlap: a1 -> a3 -> a2
        #
        # specifically:
        #   < phi[a1] | p[a3] > * dH[a3] * < p[a3] | phi[a2] >
        #
        # This matrix multiplication is semi-sparse.  It works by blocks
        # of atoms, looping only over pairs that do have nonzero
        # overlaps.  But it might be even nicer with scipy sparse.
        # This we will have to check at some point.
        #
        # The projection arrays P_aaim are distributed over the grid,
        # whereas the Hamiltonian is distributed over the band comm.
        # One could choose a set of a3 to optimize the load balance.
        # Right now the load balance will be "random" and probably
        # not very good.
        #innerloops = 0

        for (a3, a1), P1_im in kpt.P_aaim.items():
            a1M1 = M_a[a1]
            nM1 = nM_a[a1]
            a1M2 = a1M1 + nM1

            if a1M1 > Mstop or a1M2 < Mstart:
                continue

            stickout1 = max(0, Mstart - a1M1)
            stickout2 = max(0, a1M2 - Mstop)
            P1_mi = np.conj(P1_im.T[stickout1:nM1 - stickout2])
            dH_ii = y * np.asarray(unpack(dH_asp[a3][kpt.s]), dtype)
            H_mM = H_MM[a1M1 + stickout1 - Mstart:a1M2 - stickout2 - Mstart]

            P1dH_mi = np.dot(P1_mi, dH_ii)

            assert len(wfs.P_neighbors_a[a3]) > 0
            for a2 in wfs.P_neighbors_a[a3]:
                # We can use symmetry somehow.  Since the entire matrix
                # is symmetrized after the Hamiltonian is constructed,
                # at least in the non-Gamma-point case, we should do
                # so conditionally somehow.  Right now let's stay out
                # of trouble.

                # Humm.  The following works with gamma point
                # but not with kpts.  XXX take a look at this.
                # also, it doesn't work with a2 < a1 for some reason.
                #if a2 > a1:
                #    continue
                a2M1 = wfs.setups.M_a[a2]
                a2M2 = a2M1 + wfs.setups[a2].nao
                P2_im = kpt.P_aaim[(a3, a2)]
                P1dHP2_mm = np.dot(P1dH_mi, P2_im)
                H_mM[:, a2M1:a2M2] += P1dHP2_mm
                #innerloops += 1
                #if wfs.world.rank == 0:
                #    print 'y', y
                #if a1 != a2:
                #    H_MM[a2M1:a2M2, a1M1:a1M2] += P1dHP2_mm.T.conj()

