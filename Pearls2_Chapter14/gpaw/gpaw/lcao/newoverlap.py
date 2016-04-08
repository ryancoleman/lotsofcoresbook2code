import numpy as np
from ase import Atoms
from ase.calculators.neighborlist import NeighborList

from gpaw.lcao.overlap import AtomicDisplacement, TwoCenterIntegralCalculator
from gpaw.utilities.partition import EvenPartitioning

class DistsAndOffsets:
    def __init__(self, nl, spos_ac, cell_cv):
        #assert nl.bothways
        r_and_offset_aao = {}

        def add(a1, a2, R_c, offset):
            r_and_offset_aao.setdefault((a1, a2), []).append((R_c, offset))

        for a1, spos1_c in enumerate(spos_ac):
            a2_a, offsets = nl.get_neighbors(a1)
            for a2, offset in zip(a2_a, offsets):
                spos2_c = spos_ac[a2] + offset

                R_c = np.dot(spos2_c - spos1_c, cell_cv)
                add(a1, a2, R_c, offset)
                if a1 != a2 or offset.any():
                    add(a2, a1, -R_c, -offset)
        self.r_and_offset_aao = r_and_offset_aao

    def get(self, a1, a2):
        R_ca_and_offset_a = self.r_and_offset_aao.get((a1, a2))
        return R_ca_and_offset_a


def newoverlap(wfs, spos_ac):
    assert wfs.ksl.block_comm.size == wfs.gd.comm.size * wfs.bd.comm.size
    even_part = EvenPartitioning(wfs.gd.comm, #wfs.ksl.block_comm,
                                 len(wfs.atom_partition.rank_a))
    atom_partition = even_part.as_atom_partition()
    
    tci = wfs.tci

    gd = wfs.gd
    kd = wfs.kd
    nq = len(kd.ibzk_qc)

    # New neighbor list because we want it "both ways", heh.  Or do we?
    neighbors = NeighborList(tci.cutoff_a, skin=0,
                             sorted=True, self_interaction=True, bothways=False)
    atoms = Atoms('X%d' % len(tci.cutoff_a), cell=gd.cell_cv, pbc=gd.pbc_c)
    atoms.set_scaled_positions(spos_ac)
    neighbors.update(atoms)
    
    # XXX
    pcutoff_a = []
    phicutoff_a = []
    for setup in wfs.setups:
        if setup.pt_j:
            pcutoff = max([pt.get_cutoff() for pt in setup.pt_j])
        else:
            pcutoff = 0.0
        if setup.phit_j:
            phicutoff = max([phit.get_cutoff() for phit in setup.phit_j])
        else:
            phicutoff = 0.0
        pcutoff_a.append(pcutoff)
        phicutoff_a.append(phicutoff)

    # Calculate the projector--basis function overlaps:
    #
    #    a1        ~a1
    #   P      = < p   | Phi   > ,
    #    i mu       i       mu
    #
    # i.e. projector is on a1 and basis function is on what we will call a2.

    overlapcalc = TwoCenterIntegralCalculator(wfs.kd.ibzk_qc,
                                              derivative=False)
    
    P_aaqim = {} # keys: (a1, a2).  Values: matrix blocks
    dists_and_offsets = DistsAndOffsets(neighbors, spos_ac, gd.cell_cv)

    #ng = 2**extra_parameters.get('log2ng', 10)
    #transformer = FourierTransformer(rcmax, ng)
    #tsoc = TwoSiteOverlapCalculator(transformer)
    #msoc = ManySiteOverlapCalculator(tsoc, I_a, I_a)

    msoc = wfs.tci.msoc

    phit_Ij = [setup.phit_j for setup in tci.setups_I]
    l_Ij = []
    for phit_j in phit_Ij:
        l_Ij.append([phit.get_angular_momentum_number()
                     for phit in phit_j])

    pt_l_Ij = [setup.l_j for setup in tci.setups_I]        
    pt_Ij = [setup.pt_j for setup in tci.setups_I]
    phit_Ijq = msoc.transform(phit_Ij)
    pt_Ijq = msoc.transform(pt_Ij)

    #self.Theta_expansions = msoc.calculate_expansions(l_Ij, phit_Ijq,
    #                                                  l_Ij, phit_Ijq)
    #self.T_expansions = msoc.calculate_kinetic_expansions(l_Ij, phit_Ijq)
    P_expansions = msoc.calculate_expansions(pt_l_Ij, pt_Ijq,
                                             l_Ij, phit_Ijq)
    P_neighbors_a = {}
    
    for a1 in atom_partition.my_indices:
        for a2 in range(len(wfs.setups)):
            R_ca_and_offset_a = dists_and_offsets.get(a1, a2)
            if R_ca_and_offset_a is None: # No overlap between a1 and a2
                continue

            maxdistance = pcutoff_a[a1] + phicutoff_a[a2]
            expansion = P_expansions.get(a1, a2)
            P_qim = expansion.zeros((nq,), dtype=wfs.dtype)
            disp = None
            for R_c, offset in R_ca_and_offset_a:
                r = np.linalg.norm(R_c)
                if r > maxdistance:
                    continue
                
                # Below lines are meant to make use of symmetry.  Will not
                # be relevant for P.
                #remainder = (a1 + a2) % 2
                #if a1 < a2 and not remainder:
                #    continue
                # if a1 > a2 and remainder:
                #    continue
                
                phases = overlapcalc.phaseclass(overlapcalc.ibzk_qc, offset)
                disp = AtomicDisplacement(None, a1, a2, R_c, offset, phases)
                disp.evaluate_overlap(expansion, P_qim)
            
            if disp is not None: # there was at least one non-zero overlap
                assert (a1, a2) not in P_aaqim
                P_aaqim[(a1, a2)] = P_qim
                P_neighbors_a.setdefault(a1, []).append(a2)
    
    Pkeys = P_aaqim.keys()
    Pkeys.sort()

    def get_M1M2(a):
        M1 = wfs.setups.M_a[a]
        M2 = M1 + wfs.setups[a].nao
        return M1, M2
    
    oldstyle_P_aqMi = None
    if 0:#wfs.world.size == 1:
        oldstyle_P_aqMi = {}
        for a, setup in enumerate(wfs.setups):
            oldstyle_P_aqMi[a] = np.zeros((nq, wfs.setups.nao, setup.ni),
                                          dtype=wfs.dtype)
        print([(s.ni, s.nao) for s in wfs.setups])
        for a1, a2 in Pkeys:
            M1, M2 = get_M1M2(a2)
            Pconj_qmi = P_aaqim[(a1, a2)].transpose(0, 2, 1).conjugate()
            oldstyle_P_aqMi[a1][:, M1:M2, :] = Pconj_qmi
        
    # XXX mind distribution

    return P_neighbors_a, P_aaqim, oldstyle_P_aqMi
