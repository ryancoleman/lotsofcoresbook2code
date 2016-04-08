
import numpy as np

from ase.units import Bohr
from gpaw import debug
from gpaw.mpi import world
#from gpaw.utilities.dscftools import mpi_debug
mpi_debug = lambda x, ordered=True: None # silenced
from gpaw.overlap import Overlap
from gpaw.utilities import unpack
from gpaw.lfc import NewLocalizedFunctionsCollection as NewLFC

#from arbitrary_tci import projector_overlap_matrix2, generate_atomic_overlaps2

# -------------------------------------------------------------------

class PairOverlap:
    def __init__(self, gd, setups):
        self.gd = gd
        self.setups = setups
        self.ni_a = np.cumsum([0]+[setup.ni for setup in self.setups])

    def __len__(self):
        return self.ni_a[-1].item()

    def assign_atomic_pair_matrix(self, X_aa, a1, a2, dX_ii):
        X_aa[self.ni_a[a1]:self.ni_a[a1+1],\
             self.ni_a[a2]:self.ni_a[a2+1]] = dX_ii

    def extract_atomic_pair_matrix(self, X_aa, a1, a2):
        return X_aa[self.ni_a[a1]:self.ni_a[a1+1],\
                    self.ni_a[a2]:self.ni_a[a2+1]]

    def calculate_overlaps(self, spos_ac, lfc1, lfc2=None):
        raise RuntimeError('This is a virtual member function.')

    def calculate_atomic_pair_overlaps(self, lfs1, lfs2): #XXX Move some code here from above...
        raise RuntimeError('This is a virtual member function.')

class GridPairOverlap(PairOverlap):


    def calculate_overlaps(self, spos_ac, lfc1, lfc2=None):
        # CONDITION: The two sets of splines must belong to the same kpoint!

        if lfc2 is None:
            lfc2 = lfc1

        if isinstance(lfc1, NewLFC) and isinstance(lfc2, NewLFC):
            return self.calculate_overlaps2(spos_ac, lfc1, lfc2)

        assert not isinstance(lfc1, NewLFC) and not isinstance(lfc2, NewLFC)

        nproj = len(self)
        X_aa = np.zeros((nproj,nproj), dtype=float) # XXX always float?

        if debug:
            if world.rank == 0:
                print('DEBUG INFO')

            mpi_debug('lfc1.lfs_a.keys(): %s' % lfc1.lfs_a.keys())
            mpi_debug('lfc2.lfs_a.keys(): %s' % lfc2.lfs_a.keys())
            mpi_debug('N_c=%s, beg_c=%s, end_c=%s' % (self.gd.N_c,self.gd.beg_c,self.gd.end_c))

            assert len(lfc1.spline_aj) == len(lfc1.spos_ac) # not distributed
            assert len(lfc2.spline_aj) == len(lfc2.spos_ac) # not distributed
            #assert lfc1.lfs_a.keys() == lfc2.lfs_a.keys() # XXX must they be equal?!?

        # Both loops are over all atoms in all domains
        for a1, spline1_j in enumerate(lfc1.spline_aj):
            # We assume that all functions have the same cut-off:
            rcut1 = spline1_j[0].get_cutoff()
            if debug: mpi_debug('a1=%d, spos1_c=%s, rcut1=%g, ni1=%d' % (a1,spos_ac[a1],rcut1,self.setups[a1].ni))

            for a2, spline2_j in enumerate(lfc2.spline_aj):
                # We assume that all functions have the same cut-off:
                rcut2 = spline2_j[0].get_cutoff()
                if debug: mpi_debug('  a2=%d, spos2_c=%s, rcut2=%g, ni2=%d' % (a2,spos_ac[a2],rcut2,self.setups[a2].ni))

                X_ii = self.extract_atomic_pair_matrix(X_aa, a1, a2)

                b1 = 0
                for beg1_c, end1_c, sdisp1_c in self.gd.get_boxes(spos_ac[a1], rcut1, cut=False): # loop over lfs1.box_b instead?
                    if debug: mpi_debug('    b1=%d, beg1_c=%s, end1_c=%s, sdisp1_c=%s' % (b1,beg1_c,end1_c,sdisp1_c), ordered=False)

                    # Atom a1 has at least one piece so the LFC has LocFuncs
                    lfs1 = lfc1.lfs_a[a1]

                    # Similarly, the LocFuncs must have the piece at hand
                    box1 = lfs1.box_b[b1]

                    if debug:
                        assert lfs1.dtype == lfc1.dtype
                        assert self.setups[a1].ni == lfs1.ni, 'setups[%d].ni=%d, lfc1.lfs_a[%d].ni=%d' % (a1,self.setups[a1].ni,a1,lfs1.i)

                    b2 = 0
                    for beg2_c, end2_c, sdisp2_c in self.gd.get_boxes(spos_ac[a2], rcut2, cut=False): # loop over lfs2.box_b instead?
                        if debug: mpi_debug('      b2=%d, beg2_c=%s, end2_c=%s, sdisp2_c=%s' % (b2,beg2_c,end2_c,sdisp2_c), ordered=False)

                        # Atom a2 has at least one piece so the LFC has LocFuncs
                        lfs2 = lfc2.lfs_a[a2]

                        # Similarly, the LocFuncs must have the piece at hand
                        box2 = lfs2.box_b[b2]

                        if debug:
                            assert lfs2.dtype == lfc2.dtype
                            assert self.setups[a2].ni == lfs2.ni, 'setups[%d].ni=%d, lfc2.lfs_a[%d].ni=%d' % (a2,self.setups[a2].ni,a2,lfs2.ni)

                        # Find the intersection of the two boxes
                        beg_c = np.array([map(max, zip(beg1_c, beg2_c))]).ravel()
                        end_c = np.array([map(min, zip(end1_c, end2_c))]).ravel()

                        if debug: mpi_debug('        beg_c=%s, end_c=%s, size_c=%s' % (beg_c, end_c, tuple(end_c-beg_c)), ordered=False)

                        # Intersection is non-empty, add overlap contribution
                        if (beg_c < end_c).all():
                            bra_iB1 = box1.get_functions()
                            w1slice = [slice(None)]+[slice(b,e) for b,e in \
                                zip(beg_c-beg1_c, end_c-beg1_c)]

                            ket_iB2 = box2.get_functions()
                            w2slice = [slice(None)]+[slice(b,e) for b,e in \
                                zip(beg_c-beg2_c, end_c-beg2_c)]

                            X_ii += self.gd.dv * np.inner( \
                                bra_iB1[w1slice].reshape((lfs1.ni,-1)), \
                                ket_iB2[w2slice].reshape((lfs2.ni,-1))) #XXX phase factors for kpoints

                            del bra_iB1, ket_iB2

                        b2 += 1

                    b1 += 1

        self.gd.comm.sum(X_aa) # better to sum over X_ii?
        return X_aa


    def calculate_overlaps2(self, spos_ac, lfc1, lfc2=None):
        # CONDITION: The two sets of splines must belong to the same kpoint!

        if lfc2 is None:
            lfc2 = lfc1

        assert isinstance(lfc1, NewLFC) and isinstance(lfc2, NewLFC)

        nproj = len(self)
        X_aa = np.zeros((nproj,nproj), dtype=float) # XXX always float?

        if debug:
            if world.rank == 0:
                print('DEBUG INFO')

            mpi_debug('len(lfc1.sphere_a): %d, lfc1.atom_indices: %s' % (len(lfc1.sphere_a),lfc1.atom_indices))
            mpi_debug('len(lfc2.sphere_a): %d, lfc2.atom_indices: %s' % (len(lfc2.sphere_a),lfc2.atom_indices))
            mpi_debug('N_c=%s, beg_c=%s, end_c=%s' % (self.gd.N_c,self.gd.beg_c,self.gd.end_c))

        if debug:
            assert len(lfc1.sphere_a) == len(lfc2.sphere_a) # XXX must they be equal?!?

        # Both a-loops are over all relevant atoms which affect this domain
        for a1 in lfc1.atom_indices:
            sphere1 = lfc1.sphere_a[a1]

            # We assume that all functions have the same cut-off:
            spline1_j = sphere1.spline_j
            rcut1 = spline1_j[0].get_cutoff()
            if debug: mpi_debug('a1=%d, spos1_c=%s, rcut1=%g, ni1=%d' % (a1,spos_ac[a1],rcut1,self.setups[a1].ni), ordered=False)

            for a2 in lfc2.atom_indices:
                sphere2 = lfc2.sphere_a[a2]

                # We assume that all functions have the same cut-off:
                spline2_j = sphere2.spline_j
                rcut2 = spline2_j[0].get_cutoff()
                if debug: mpi_debug('  a2=%d, spos2_c=%s, rcut2=%g, ni2=%d' % (a2,spos_ac[a2],rcut2,self.setups[a2].ni), ordered=False)

                X_ii = self.extract_atomic_pair_matrix(X_aa, a1, a2)

                b1 = 0
                for beg1_c, end1_c, sdisp1_c in self.gd.get_boxes(spos_ac[a1], rcut1, cut=False): # loop over lfs1.box_b instead?
                    if debug: mpi_debug('    b1=%d, beg1_c=%s, end1_c=%s, sdisp1_c=%s' % (b1,beg1_c,end1_c,sdisp1_c), ordered=False)

                    b2 = 0
                    for beg2_c, end2_c, sdisp2_c in self.gd.get_boxes(spos_ac[a2], rcut2, cut=False): # loop over lfs2.box_b instead?
                        if debug: mpi_debug('      b2=%d, beg2_c=%s, end2_c=%s, sdisp2_c=%s' % (b2,beg2_c,end2_c,sdisp2_c), ordered=False)

                        # Find the intersection of the two boxes
                        beg_c = np.array([map(max, zip(beg1_c, beg2_c))]).ravel()
                        end_c = np.array([map(min, zip(end1_c, end2_c))]).ravel()

                        if debug: mpi_debug('        beg_c=%s, end_c=%s, size_c=%s' % (beg_c, end_c, tuple(end_c-beg_c)), ordered=False)

                        # Intersection is non-empty, add overlap contribution
                        if (beg_c < end_c).all():
                            i1 = 0
                            for j1, spline1 in enumerate(spline1_j):
                                bra1_mB = spline1.get_functions(self.gd, \
                                    beg_c, end_c, spos_ac[a1]-sdisp1_c)
                                nm1 = bra1_mB.shape[0]

                                i2 = 0
                                for j2, spline2 in enumerate(spline2_j):
                                    ket2_mB = spline2.get_functions(self.gd, \
                                        beg_c, end_c, spos_ac[a2]-sdisp2_c)
                                    nm2 = ket2_mB.shape[0]

                                    X_mm = X_ii[i1:i1+nm1,i2:i2+nm2]
                                    X_mm += self.gd.dv * np.inner( \
                                        bra1_mB.reshape((nm1,-1)), \
                                        ket2_mB.reshape((nm2,-1))) #XXX phase factors for kpoints

                                    del ket2_mB
                                    i2 += nm2

                                del bra1_mB
                                i1 += nm1

                        b2 += 1

                    b1 += 1

        self.gd.comm.sum(X_aa) # better to sum over X_ii?
        return X_aa


"""
def extract_projectors(gd, pt, a):
    print 'pt.lfs_a: %s' % pt.lfs_a.keys()
    lfs = pt.lfs_a[a]
    assert len(lfs.box_b) == 1
    box = lfs.box_b[0]

    pt_iB = box.get_functions()
    assert pt_iB.shape[0] == lfs.ni
    shape = pt_iB.shape[1:]

    work_iG = gd.zeros(lfs.ni, dtype=pt.dtype)
    print 'a=%d, ni=%d, pt_iB.shape' % (a,lfs.ni), pt_iB.shape

    # We assume that all functions have the same cut-off:
    rcut = pt.spline_aj[a][0].get_cutoff()
    gdbox_b = gd.get_boxes(pt.spos_ac[a], rcut, cut=False)
    assert len(gdbox_b) == 1, 'What?'
    beg_c, end_c, sdisp_c = gdbox_b[0]
    work_iG[:,beg_c[0]:end_c[0],beg_c[1]:end_c[1],beg_c[2]:end_c[2]] = pt_iB

    #for i,spline in enumerate(pt.spline_aj[a]): #XXX wrong, loop over pt_iB!!!
    #    rcut = spline.get_cutoff()
    #    gdbox_b = gd.get_boxes(pt.spos_ac[a], rcut, cut=False)
    #    assert len(gdbox_b) == 1, 'What?'
    #    beg_c, end_c, sdisp_c = gdbox_b[0]
    #    print 'a=%d, pt_iB[%d].shape' % (a,i), pt_iB[i].shape
    #    print 'a=%d, work_iG[%d].shape' % (a,i), work_iG[i,beg_c[0]:end_c[0],beg_c[1]:end_c[1],beg_c[2]:end_c[2]].shape
    #    work_iG[i,beg_c[0]:end_c[0],beg_c[1]:end_c[1],beg_c[2]:end_c[2]] = pt_iB[i]

    return work_iG

def overlap_projectors(gd, pt, setups):
    work_aiG = {}
    for a in range(len(setups)):
        work_aiG[a] = extract_projectors(gd,pt,a)

    ni_a = np.cumsum([0]+[setup.ni for setup in setups])
    nproj = ni_a[-1]
    B_aa = np.zeros((nproj,nproj), dtype=float)
    for a1, work1_iG in work_aiG.items():
        for a2, work2_iG in work_aiG.items():
            B_ii = np.zeros((setups[a1].ni,setups[a2].ni), dtype=float)
            for i1, work1_G in enumerate(work1_iG):
                for i2, work2_G in enumerate(work2_iG):
                    B_ii[i1,i2] = np.dot(work1_G.flat, work2_G.flat)*gd.dv
            B_aa[ni_a[a1]:ni_a[a1+1], ni_a[a2]:ni_a[a2+1]] = B_ii
    return B_aa
"""

class ProjectorPairOverlap(Overlap, GridPairOverlap):
    """
    TODO
    """
    def __init__(self, wfs, atoms):
        """TODO

        Attributes:

        ============  ======================================================
        ``B_aa``      < p_i^a | p_i'^a' >
        ``xO_aa``     TODO
        ``dC_aa``     TODO
        ``xC_aa``     TODO
        ============  ======================================================
        """

        Overlap.__init__(self, wfs.orthoksl, wfs.timer)
        GridPairOverlap.__init__(self, wfs.gd, wfs.setups)
        self.natoms = len(atoms)
        if debug:
            assert len(self.setups) == self.natoms
        self.update(wfs, atoms)

    def update(self, wfs, atoms):
        self.timer.start('Update two-center overlap')

        nproj = len(self)

        """
        self.B_aa = np.zeros((nproj, nproj), dtype=float) #always float?
        for a1,setup1 in enumerate(self.setups):
            for a2 in wfs.pt.my_atom_indices:
                setup2 = self.setups[a2]
                R = (atoms[a1].get_position() - atoms[a2].get_position()) / Bohr

                if a1 == a2:
                    B_ii = setup1.B_ii
                else:
                    B_ii = projector_overlap_matrix2(setup1, setup2, R)
                #if a1 < a2:
                #    B_ii = projector_overlap_matrix2(setup1, setup2, R)
                #elif a1 == a2:
                #    B_ii = setup1.B_ii
                #else:
                #    B_ii = self.B_aa[ni_a[a2]:ni_a[a2+1], ni_a[a1]:ni_a[a1+1]].T

                #self.B_aa[self.ni_a[a1]:self.ni_a[a1+1], \
                #          self.ni_a[a2]:self.ni_a[a2+1]] = B_ii
                self.assign_atomic_pair_matrix(self.B_aa, a1, a2, B_ii)
        self.gd.comm.sum(self.B_aa) #TODO too heavy?
        """
        #self.B_aa = overlap_projectors(wfs.gd, wfs.pt, wfs.setups)

        #spos_ac = wfs.pt.spos_ac # not in NewLFC
        spos_ac = atoms.get_scaled_positions() % 1.0
        self.B_aa = self.calculate_overlaps(spos_ac, wfs.pt)

        # Create two-center (block-diagonal) coefficients for overlap operator
        dO_aa = np.zeros((nproj, nproj), dtype=float) #always float?
        for a,setup in enumerate(self.setups):
            self.assign_atomic_pair_matrix(dO_aa, a, a, setup.dO_ii)

        # Calculate two-center rotation matrix for overlap projections
        self.xO_aa = self.get_rotated_coefficients(dO_aa)

        # Calculate two-center coefficients for inverse overlap operator
        lhs_aa = np.eye(nproj) + self.xO_aa
        rhs_aa = -dO_aa
        self.dC_aa = np.linalg.solve(lhs_aa.T, rhs_aa.T).T #TODO parallel

        # Calculate two-center rotation matrix for inverse overlap projections
        self.xC_aa = self.get_rotated_coefficients(self.dC_aa)

        self.timer.stop('Update two-center overlap')

    def get_rotated_coefficients(self, X_aa):
        """Rotate two-center projector expansion coefficients with
        the projector-projector overlap integrals as basis.

        Performs the following operation and returns the result::

                    ---
            a1,a3   \       a1    a2     a2,a3
           Y      =  )  <  p   | p   >  X
            i1,i3   /       i1    i2     i2,i3
                    ---
                   a2,i2
        """
        return np.dot(self.B_aa, X_aa)

    def apply_to_atomic_matrices(self, dI_asp, P_axi, wfs, kpt, shape=()):

        self.timer.start('Update two-center projections')

        nproj = len(self)
        dI_aa = np.zeros((nproj, nproj), dtype=float) #always float?

        for a, dI_sp in dI_asp.items():
            dI_p = dI_sp[kpt.s]
            dI_ii = unpack(dI_p)
            self.assign_atomic_pair_matrix(dI_aa, a, a, dI_ii)
        self.gd.comm.sum(dI_aa) #TODO too heavy?

        dM_aa = self.get_rotated_coefficients(dI_aa)
        Q_axi = wfs.pt.dict(shape, zero=True)
        for a1 in range(self.natoms):
            if a1 in Q_axi.keys():
                Q_xi = Q_axi[a1]
            else:
                # Atom a1 is not in domain so allocate a temporary buffer
                Q_xi = np.zeros(shape+(self.setups[a1].ni,), dtype=wfs.pt.dtype) #TODO
            for a2, P_xi in P_axi.items():
                dM_ii = self.extract_atomic_pair_matrix(dM_aa, a1, a2)
                Q_xi += np.dot(P_xi, dM_ii.T) #sum over a2 and last i in dM_ii
            self.gd.comm.sum(Q_xi)

        self.timer.stop('Update two-center projections')

        return Q_axi

    def apply(self, a_xG, b_xG, wfs, kpt, calculate_P_ani=True,
              extrapolate_P_ani=False):
        """Apply the overlap operator to a set of vectors.

        Parameters
        ==========
        a_nG: ndarray
            Set of vectors to which the overlap operator is applied.
        b_nG: ndarray, output
            Resulting S times a_nG vectors.
        kpt: KPoint object
            k-point object defined in kpoint.py.
        calculate_P_ani: bool
            When True, the integrals of projector times vectors
            P_ni = <p_i | a_nG> are calculated.
            When False, existing P_ani are used
        extrapolate_P_ani: bool
            When True, the integrals of projector times vectors#XXX TODO
            P_ni = <p_i | a_nG> are calculated.
            When False, existing P_ani are used

        """

        self.timer.start('Apply overlap')
        b_xG[:] = a_xG
        shape = a_xG.shape[:-3]
        P_axi = wfs.pt.dict(shape)

        if calculate_P_ani:
            wfs.pt.integrate(a_xG, P_axi, kpt.q)
        else:
            for a,P_ni in kpt.P_ani.items():
                P_axi[a][:] = P_ni

        Q_axi = wfs.pt.dict(shape)
        for a, Q_xi in Q_axi.items():
            Q_xi[:] = np.dot(P_axi[a], self.setups[a].dO_ii)

        wfs.pt.add(b_xG, Q_axi, kpt.q)
        self.timer.stop('Apply overlap')

        if extrapolate_P_ani:
            for a1 in range(self.natoms):
                if a1 in Q_axi.keys():
                    Q_xi = Q_axi[a1]
                    Q_xi[:] = P_axi[a1]
                else:
                    # Atom a1 is not in domain so allocate a temporary buffer
                    Q_xi = np.zeros(shape+(self.setups[a1].ni,), dtype=wfs.pt.dtype) #TODO
                for a2, P_xi in P_axi.items():
                    # xO_aa are the overlap extrapolators across atomic pairs
                    xO_ii = self.extract_atomic_pair_matrix(self.xO_aa, a1, a2)
                    Q_xi += np.dot(P_xi, xO_ii.T) #sum over a2 and last i in xO_ii
                self.gd.comm.sum(Q_xi)

            return Q_axi
        else:
            return P_axi

    def apply_inverse(self, a_xG, b_xG, wfs, kpt, calculate_P_ani=True,
                      extrapolate_P_ani=False):

        self.timer.start('Apply inverse overlap')
        b_xG[:] = a_xG
        shape = a_xG.shape[:-3]
        P_axi = wfs.pt.dict(shape)

        if calculate_P_ani:
            wfs.pt.integrate(a_xG, P_axi, kpt.q)
        else:
            for a,P_ni in kpt.P_ani.items():
                P_axi[a][:] = P_ni

        Q_axi = wfs.pt.dict(shape, zero=True)
        for a1 in range(self.natoms):
            if a1 in Q_axi.keys():
                Q_xi = Q_axi[a1]
            else:
                # Atom a1 is not in domain so allocate a temporary buffer
                Q_xi = np.zeros(shape+(self.setups[a1].ni,), dtype=wfs.pt.dtype) #TODO
            for a2, P_xi in P_axi.items():
                # dC_aa are the inverse coefficients across atomic pairs
                dC_ii = self.extract_atomic_pair_matrix(self.dC_aa, a1, a2)
                Q_xi += np.dot(P_xi, dC_ii.T) #sum over a2 and last i in dC_ii
            self.gd.comm.sum(Q_xi)

        wfs.pt.add(b_xG, Q_axi, kpt.q)
        self.timer.stop('Apply inverse overlap')

        if extrapolate_P_ani:
            for a1 in range(self.natoms):
                if a1 in Q_axi.keys():
                    Q_xi = Q_axi[a1]
                    Q_xi[:] = P_axi[a1]
                else:
                    # Atom a1 is not in domain so allocate a temporary buffer
                    Q_xi = np.zeros(shape+(self.setups[a1].ni,), dtype=wfs.pt.dtype) #TODO
                for a2, P_xi in P_axi.items():
                    # xC_aa are the inverse extrapolators across atomic pairs
                    xC_ii = self.extract_atomic_pair_matrix(self.xC_aa, a1, a2)
                    Q_xi += np.dot(P_xi, xC_ii.T) #sum over a2 and last i in xC_ii
                self.gd.comm.sum(Q_xi)

            return Q_axi
        else:
            return P_axi





