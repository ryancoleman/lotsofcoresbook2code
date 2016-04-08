# -*- coding: utf-8 -*-
# Copyright (C) 2003  CAMP
# Please see the accompanying LICENSE file for further information.

"""This module defines a Hamiltonian."""

from math import pi, sqrt

import numpy as np

from gpaw.poisson import PoissonSolver
from gpaw.transformers import Transformer
from gpaw.lfc import LFC
from gpaw.utilities import pack2, unpack, unpack2
from gpaw.io import read_atomic_matrices
from gpaw.utilities.partition import AtomPartition, AtomicMatrixDistributor


class Hamiltonian:
    """Hamiltonian object.

    Attributes:
     =============== =====================================================
     ``xc``          ``XC3DGrid`` object.
     ``poisson``     ``PoissonSolver``.
     ``gd``          Grid descriptor for coarse grids.
     ``finegd``      Grid descriptor for fine grids.
     ``restrict``    Function for restricting the effective potential.
     =============== =====================================================

    Soft and smooth pseudo functions on uniform 3D grids:
     ========== =========================================
     ``vHt_g``  Hartree potential on the fine grid.
     ``vt_sG``  Effective potential on the coarse grid.
     ``vt_sg``  Effective potential on the fine grid.
     ========== =========================================

    Energy contributions and forces:

    =========== ==========================================
                Description
    =========== ==========================================
    ``Ekin``    Kinetic energy.
    ``Epot``    Potential energy.
    ``Etot``    Total energy.
    ``Exc``     Exchange-Correlation energy.
    ``Eext``    Energy of external potential
    ``Eref``    Reference energy for all-electron atoms.
    ``S``       Entropy.
    ``Ebar``    Should be close to zero!
    =========== ==========================================

    """

    def __init__(self, gd, finegd, nspins, setups, timer, xc,
                 world, kptband_comm, vext=None, collinear=True):
        """Create the Hamiltonian."""
        self.gd = gd
        self.finegd = finegd
        self.nspins = nspins
        self.setups = setups
        self.timer = timer
        self.xc = xc
        self.collinear = collinear
        self.ncomp = 2 - int(collinear)
        self.ns = self.nspins * self.ncomp**2
        self.world = world
        self.kptband_comm = kptband_comm
        
        self.dH_asp = None

        # The external potential
        self.vext = vext

        self.vt_sG = None
        self.vHt_g = None
        self.vt_sg = None

        self.rank_a = None
        self.atom_partition = None

        self.Ekin0 = None
        self.Ekin = None
        self.Epot = None
        self.Ebar = None
        self.Eext = None
        self.Exc = None
        self.Etot = None
        self.S = None

        self.ref_vt_sG = None
        self.ref_dH_asp = None


    def summary(self, fd):
        fd.write('XC and Coulomb potentials evaluated on a %d*%d*%d grid\n' %
                 tuple(self.finegd.N_c))

    def set_positions(self, spos_ac, rank_a):
        atom_partition = AtomPartition(self.gd.comm, rank_a)
        
        self.spos_ac = spos_ac
        self.vbar.set_positions(spos_ac)
        self.xc.set_positions(spos_ac)
        
        # If both old and new atomic ranks are present, start a blank dict if
        # it previously didn't exist but it will needed for the new atoms.
        # XXX what purpose does this serve?  In what case does it happen?
        # How would one even go about figuring it out?  Why does it all have
        # to be so unreadable? -Ask
        #
        if (self.rank_a is not None and
            self.dH_asp is None and (rank_a == self.gd.comm.rank).any()):
            self.dH_asp = {}
            
        if self.rank_a is not None and self.dH_asp is not None:
            self.timer.start('Redistribute')
            def get_empty(a):
                ni = self.setups[a].ni
                return np.empty((self.ns, ni * (ni + 1) // 2))
            self.atom_partition.redistribute(atom_partition, self.dH_asp,
                                             get_empty)
            self.timer.stop('Redistribute')

        self.rank_a = rank_a
        self.atom_partition = atom_partition
        self.dh_distributor = AtomicMatrixDistributor(atom_partition,
                                                      self.setups,
                                                      self.kptband_comm,
                                                      self.ns)


    def aoom(self, DM, a, l, scale=1):
        """Atomic Orbital Occupation Matrix.
        
        Determine the Atomic Orbital Occupation Matrix (aoom) for a
        given l-quantum number.
        
        This operation, takes the density matrix (DM), which for
        example is given by unpack2(D_asq[i][spin]), and corrects for
        the overlap between the selected orbitals (l) upon which the
        the density is expanded (ex <p|p*>,<p|p>,<p*|p*> ).

        Returned is only the "corrected" part of the density matrix,
        which represents the orbital occupation matrix for l=2 this is
        a 5x5 matrix.
        """
        S=self.setups[a]
        l_j = S.l_j
        lq  = S.lq
        nl  = np.where(np.equal(l_j, l))[0]
        V = np.zeros(np.shape(DM))
        if len(nl) == 2:
            aa = (nl[0])*len(l_j)-((nl[0]-1)*(nl[0])/2)
            bb = (nl[1])*len(l_j)-((nl[1]-1)*(nl[1])/2)
            ab = aa+nl[1]-nl[0]
            
            if not scale:
                lq_a  = lq[aa]
                lq_ab = lq[ab]
                lq_b  = lq[bb]
            else:
                lq_a  = 1
                lq_ab = lq[ab]/lq[aa]
                lq_b  = lq[bb]/lq[aa]
 
            # and the correct entrances in the DM
            nn = (2*np.array(l_j)+1)[0:nl[0]].sum()
            mm = (2*np.array(l_j)+1)[0:nl[1]].sum()
            
            # finally correct and add the four submatrices of NC_DM
            A = DM[nn:nn+2*l+1,nn:nn+2*l+1]*(lq_a)
            B = DM[nn:nn+2*l+1,mm:mm+2*l+1]*(lq_ab)
            C = DM[mm:mm+2*l+1,nn:nn+2*l+1]*(lq_ab)
            D = DM[mm:mm+2*l+1,mm:mm+2*l+1]*(lq_b)
            
            V[nn:nn+2*l+1,nn:nn+2*l+1]=+(lq_a)
            V[nn:nn+2*l+1,mm:mm+2*l+1]=+(lq_ab)
            V[mm:mm+2*l+1,nn:nn+2*l+1]=+(lq_ab)
            V[mm:mm+2*l+1,mm:mm+2*l+1]=+(lq_b)
 
            return  A+B+C+D, V
        else:
            nn =(2*np.array(l_j)+1)[0:nl[0]].sum()
            A=DM[nn:nn+2*l+1,nn:nn+2*l+1]*lq[-1]
            V[nn:nn+2*l+1,nn:nn+2*l+1]=+lq[-1]
            return A,V

    def update(self, density):
        """Calculate effective potential.

        The XC-potential and the Hartree potential are evaluated on
        the fine grid, and the sum is then restricted to the coarse
        grid."""

        self.timer.start('Hamiltonian')

        if self.vt_sg is None:
            self.timer.start('Initialize Hamiltonian')
            self.vt_sg = self.finegd.empty(self.ns)
            self.vHt_g = self.finegd.zeros()
            self.vt_sG = self.gd.empty(self.ns)
            self.poisson.initialize()
            self.timer.stop('Initialize Hamiltonian')

        Ekin, Epot, Ebar, Eext, Exc, W_aL = \
            self.update_pseudo_potential(density)

        self.timer.start('Atomic')
        self.dH_asp = None # XXXX

        dH_asp = {}
        for a, D_sp in density.D_asp.items():
            W_L = W_aL[a]
            setup = self.setups[a]

            D_p = D_sp[:self.nspins].sum(0)
            dH_p = (setup.K_p + setup.M_p +
                    setup.MB_p + 2.0 * np.dot(setup.M_pp, D_p) +
                    np.dot(setup.Delta_pL, W_L))
            Ekin += np.dot(setup.K_p, D_p) + setup.Kc
            Ebar += setup.MB + np.dot(setup.MB_p, D_p)
            Epot += setup.M + np.dot(D_p, (setup.M_p +
                                           np.dot(setup.M_pp, D_p)))

            if self.vext is not None:
                vext = self.vext.get_taylor(spos_c=self.spos_ac[a, :])
                # Tailor expansion to the zeroth order
                Eext += vext[0][0] * (sqrt(4 * pi) * density.Q_aL[a][0]
                                      + setup.Z)
                dH_p += vext[0][0] * sqrt(4 * pi) * setup.Delta_pL[:, 0]
                if len(vext) > 1:
                    # Tailor expansion to the first order
                    Eext += sqrt(4 * pi / 3) * np.dot(vext[1],
                                                      density.Q_aL[a][1:4])
                    # there must be a better way XXXX
                    Delta_p1 = np.array([setup.Delta_pL[:, 1],
                                          setup.Delta_pL[:, 2],
                                          setup.Delta_pL[:, 3]])
                    dH_p += sqrt(4 * pi / 3) * np.dot(vext[1], Delta_p1)

            dH_asp[a] = dH_sp = np.zeros_like(D_sp)

            if setup.HubU is not None:
                assert self.collinear
                nspins = len(D_sp)
                
                l_j = setup.l_j
                l   = setup.Hubl
                scale = setup.Hubs
                nl  = np.where(np.equal(l_j,l))[0]
                nn  = (2*np.array(l_j)+1)[0:nl[0]].sum()
                
                for D_p, H_p in zip(D_sp, dH_asp[a]):
                    [N_mm,V] =self.aoom(unpack2(D_p),a,l, scale)
                    N_mm = N_mm / 2 * nspins
                     
                    Eorb = setup.HubU / 2. * (N_mm - np.dot(N_mm,N_mm)).trace()
                    Vorb = setup.HubU * (0.5 * np.eye(2*l+1) - N_mm)
                    Exc += Eorb
                    if nspins == 1:
                        # add contribution of other spin manyfold
                        Exc += Eorb
                    
                    if len(nl)==2:
                        mm  = (2*np.array(l_j)+1)[0:nl[1]].sum()
                        
                        V[nn:nn+2*l+1,nn:nn+2*l+1] *= Vorb
                        V[mm:mm+2*l+1,nn:nn+2*l+1] *= Vorb
                        V[nn:nn+2*l+1,mm:mm+2*l+1] *= Vorb
                        V[mm:mm+2*l+1,mm:mm+2*l+1] *= Vorb
                    else:
                        V[nn:nn+2*l+1,nn:nn+2*l+1] *= Vorb
                    
                    Htemp = unpack(H_p)
                    Htemp += V
                    H_p[:] = pack2(Htemp)

            dH_sp[:self.nspins] += dH_p
            if self.ref_dH_asp:
                dH_sp += self.ref_dH_asp[a]
            # We are not yet done with dH_sp; still need XC correction below

        Ddist_asp = self.dh_distributor.distribute(density.D_asp)
        
        dHdist_asp = {}
        Exca = 0.0
        self.timer.start('XC Correction')
        for a, D_sp in Ddist_asp.items():
            setup = self.setups[a]
            dH_sp = np.zeros_like(D_sp)
            Exca += self.xc.calculate_paw_correction(setup, D_sp, dH_sp, a=a)
            # XXX Exc are added on the "wrong" distribution; sum only works
            # when gd.comm and distribution comm are the same
            dHdist_asp[a] = dH_sp
        self.timer.stop('XC Correction')
        
        dHdist_asp = self.dh_distributor.collect(dHdist_asp)

        # Exca has contributions from all cores so modify it so it is
        # parallel in the same way as the other energies.
        Exca = self.world.sum(Exca)
        if self.gd.comm.rank == 0:
            Exc += Exca
        
        assert len(dHdist_asp) == len(self.atom_partition.my_indices)

        for a, D_sp in density.D_asp.items():
            dH_sp = dH_asp[a]
            dH_sp += dHdist_asp[a]
            Ekin -= (D_sp * dH_sp).sum()  # NCXXX
        self.dH_asp = dH_asp
        self.timer.stop('Atomic')

        # Make corrections due to non-local xc:
        #xcfunc = self.xc.xcfunc
        self.Enlxc = 0.0  # XXXxcfunc.get_non_local_energy()
        Ekin += self.xc.get_kinetic_energy_correction() / self.gd.comm.size
        
        energies = np.array([Ekin, Epot, Ebar, Eext, Exc])
        self.timer.start('Communicate energies')
        self.gd.comm.sum(energies)
        # Make sure that all CPUs have the same energies
        self.world.broadcast(energies, 0)
        self.timer.stop('Communicate energies')
        (self.Ekin0, self.Epot, self.Ebar, self.Eext, self.Exc) = energies

        #self.Exc += self.Enlxc
        #self.Ekin0 += self.Enlkin

        self.timer.stop('Hamiltonian')

    def get_energy(self, occupations):
        self.Ekin = self.Ekin0 + occupations.e_band
        self.S = occupations.e_entropy

        # Total free energy:
        self.Etot = (self.Ekin + self.Epot + self.Eext +
                     self.Ebar + self.Exc - self.S)

        return self.Etot

    def linearize_to_xc(self, new_xc, density):
        # Store old hamiltonian
        ref_vt_sG = self.vt_sG.copy()
        ref_dH_asp = {}
        for a, dH_sp in self.dH_asp.items():
            ref_dH_asp[a] = dH_sp.copy()
        self.xc = new_xc
        self.xc.set_positions(self.spos_ac)
        self.update(density)

        ref_vt_sG -= self.vt_sG
        for a, dH_sp in self.dH_asp.items():
            ref_dH_asp[a] -= dH_sp
        self.ref_vt_sG = ref_vt_sG
        self.ref_dH_asp = ref_dH_asp



    def calculate_forces(self, dens, F_av):
        ghat_aLv = dens.ghat.dict(derivative=True)
        nct_av = dens.nct.dict(derivative=True)
        vbar_av = self.vbar.dict(derivative=True)

        self.calculate_forces2(dens, ghat_aLv, nct_av, vbar_av)

        # Force from compensation charges:
        for a, dF_Lv in ghat_aLv.items():
            F_av[a] += np.dot(dens.Q_aL[a], dF_Lv)

        # Force from smooth core charge:
        for a, dF_v in nct_av.items():
            F_av[a] += dF_v[0]

        # Force from zero potential:
        for a, dF_v in vbar_av.items():
            F_av[a] += dF_v[0]

        self.xc.add_forces(F_av)
        self.gd.comm.sum(F_av, 0)

    def apply_local_potential(self, psit_nG, Htpsit_nG, s):
        """Apply the Hamiltonian operator to a set of vectors.

        XXX Parameter description is deprecated!
        
        Parameters:

        a_nG: ndarray
            Set of vectors to which the overlap operator is applied.
        b_nG: ndarray, output
            Resulting H times a_nG vectors.
        kpt: KPoint object
            k-point object defined in kpoint.py.
        calculate_projections: bool
            When True, the integrals of projector times vectors
            P_ni = <p_i | a_nG> are calculated.
            When False, existing P_uni are used
        local_part_only: bool
            When True, the non-local atomic parts of the Hamiltonian
            are not applied and calculate_projections is ignored.
        
        """
        vt_G = self.vt_sG[s]
        if psit_nG.ndim == 3:
            Htpsit_nG += psit_nG * vt_G
        else:
            for psit_G, Htpsit_G in zip(psit_nG, Htpsit_nG):
                Htpsit_G += psit_G * vt_G

    def apply(self, a_xG, b_xG, wfs, kpt, calculate_P_ani=True):
        """Apply the Hamiltonian operator to a set of vectors.

        Parameters:

        a_nG: ndarray
            Set of vectors to which the overlap operator is applied.
        b_nG: ndarray, output
            Resulting S times a_nG vectors.
        wfs: WaveFunctions
            Wave-function object defined in wavefunctions.py
        kpt: KPoint object
            k-point object defined in kpoint.py.
        calculate_P_ani: bool
            When True, the integrals of projector times vectors
            P_ni = <p_i | a_nG> are calculated.
            When False, existing P_ani are used
        
        """

        wfs.kin.apply(a_xG, b_xG, kpt.phase_cd)
        self.apply_local_potential(a_xG, b_xG, kpt.s)
        shape = a_xG.shape[:-3]
        P_axi = wfs.pt.dict(shape)

        if calculate_P_ani:  # TODO calculate_P_ani=False is experimental
            wfs.pt.integrate(a_xG, P_axi, kpt.q)
        else:
            for a, P_ni in kpt.P_ani.items():
                P_axi[a][:] = P_ni

        for a, P_xi in P_axi.items():
            dH_ii = unpack(self.dH_asp[a][kpt.s])
            P_axi[a] = np.dot(P_xi, dH_ii)
        wfs.pt.add(b_xG, P_axi, kpt.q)

    def get_xc_difference(self, xc, density):
        """Calculate non-selfconsistent XC-energy difference."""
        if density.nt_sg is None:
            density.interpolate_pseudo_density()
        nt_sg = density.nt_sg
        if hasattr(xc, 'hybrid'):
            xc.calculate_exx()
        Exc = xc.calculate(density.finegd, nt_sg) / self.gd.comm.size
        for a, D_sp in density.D_asp.items():
            setup = self.setups[a]
            Exc += xc.calculate_paw_correction(setup, D_sp)
        Exc = self.gd.comm.sum(Exc)
        return Exc - self.Exc

    def estimate_memory(self, mem):
        nbytes = self.gd.bytecount()
        nfinebytes = self.finegd.bytecount()
        arrays = mem.subnode('Arrays', 0)
        arrays.subnode('vHt_g', nfinebytes)
        arrays.subnode('vt_sG', self.nspins * nbytes)
        arrays.subnode('vt_sg', self.nspins * nfinebytes)
        self.xc.estimate_memory(mem.subnode('XC'))
        self.poisson.estimate_memory(mem.subnode('Poisson'))
        self.vbar.estimate_memory(mem.subnode('vbar'))

    def read(self, reader, parallel):
        self.Ekin = reader['Ekin']
        self.Epot = reader['Epot']
        self.Ebar = reader['Ebar']
        try:
            self.Eext = reader['Eext']
        except (AttributeError, KeyError):
            self.Eext = 0.0
        self.Exc = reader['Exc']
        self.S = reader['S']
        self.Etot = reader.get('PotentialEnergy',
                               broadcast=True) - 0.5 * self.S

        if not reader.has_array('PseudoPotential'):
            return

        hdf5 = hasattr(reader, 'hdf5')
        version = reader['version']

        # Read pseudo potential on the coarse grid
        # and broadcast on kpt/band comm:
        if version > 0.3:
            self.vt_sG = self.gd.empty(self.nspins)
            if hdf5:
                indices = [slice(0, self.nspins), ] + self.gd.get_slice()
                do_read = (self.kptband_comm.rank == 0)
                reader.get('PseudoPotential', out=self.vt_sG,
                           parallel=parallel,
                           read=do_read, *indices)  # XXX read=?
                self.kptband_comm.broadcast(self.vt_sG, 0)
            else:
                for s in range(self.nspins):
                    self.gd.distribute(reader.get('PseudoPotential', s),
                                       self.vt_sG[s])

        # Read non-local part of hamiltonian
        self.dH_asp = {}
        natoms = len(self.setups)
        self.rank_a = np.zeros(natoms, int)
        if version > 0.3:
            all_H_sp = reader.get('NonLocalPartOfHamiltonian', broadcast=True)

        if self.gd.comm.rank == 0 and version > 0.3:
            self.dH_asp = read_atomic_matrices(all_H_sp, self.setups)


class RealSpaceHamiltonian(Hamiltonian):
    def __init__(self, gd, finegd, nspins, setups, timer, xc, world,
                 kptband_comm, vext=None, collinear=True, psolver=None,
                 stencil=3):
        Hamiltonian.__init__(self, gd, finegd, nspins, setups, timer, xc,
                             world, kptband_comm, vext, collinear)

        # Solver for the Poisson equation:
        if psolver is None:
            psolver = PoissonSolver(nn=3, relax='J')
        self.poisson = psolver
        self.poisson.set_grid_descriptor(finegd)

        # Restrictor function for the potential:
        self.restrictor = Transformer(self.finegd, self.gd, stencil)
        self.restrict = self.restrictor.apply

        self.vbar = LFC(self.finegd, [[setup.vbar] for setup in setups],
                        forces=True)
        self.vbar_g = None

    def summary(self, fd):
        Hamiltonian.summary(self, fd)

        degree = self.restrictor.nn * 2 - 1
        name = ['linear', 'cubic', 'quintic', 'heptic'][degree // 2]
        fd.write('Interpolation: tri-%s ' % name +
                 '(%d. degree polynomial)\n' % degree)

        fd.write('Poisson solver: %s\n' % self.poisson.get_description())

    def set_positions(self, spos_ac, rank_a):
        Hamiltonian.set_positions(self, spos_ac, rank_a)
        if self.vbar_g is None:
            self.vbar_g = self.finegd.empty()
        self.vbar_g[:] = 0.0
        self.vbar.add(self.vbar_g)

    def update_pseudo_potential(self, density):
        self.timer.start('vbar')
        Ebar = self.finegd.integrate(self.vbar_g, density.nt_g,
                                     global_integral=False)

        vt_g = self.vt_sg[0]
        vt_g[:] = self.vbar_g
        self.timer.stop('vbar')

        Eext = 0.0
        if self.vext is not None:
            assert self.collinear
            vt_g += self.vext.get_potential(self.finegd)
            Eext = self.finegd.integrate(vt_g, density.nt_g,
                                         global_integral=False) - Ebar

        self.vt_sg[1:self.nspins] = vt_g

        self.vt_sg[self.nspins:] = 0.0
            
        self.timer.start('XC 3D grid')
        Exc = self.xc.calculate(self.finegd, density.nt_sg, self.vt_sg)
        Exc /= self.gd.comm.size
        self.timer.stop('XC 3D grid')

        self.timer.start('Poisson')
        # npoisson is the number of iterations:
        self.npoisson = self.poisson.solve(self.vHt_g, density.rhot_g,
                                           charge=-density.charge)
        self.timer.stop('Poisson')

        self.timer.start('Hartree integrate/restrict')
        Epot = 0.5 * self.finegd.integrate(self.vHt_g, density.rhot_g,
                                           global_integral=False)

        Ekin = 0.0
        s = 0
        for s, (vt_g, vt_G, nt_G) in enumerate(zip(self.vt_sg, self.vt_sG, density.nt_sG)):
            if s < self.nspins:
                vt_g += self.vHt_g

            self.restrict(vt_g, vt_G)
            if self.ref_vt_sG is not None:
                vt_G += self.ref_vt_sG[s]

            if s < self.nspins:
                Ekin -= self.gd.integrate(vt_G, nt_G - density.nct_G,
                                          global_integral=False)
            else:
                Ekin -= self.gd.integrate(vt_G, nt_G, global_integral=False)
            s += 1

        self.timer.stop('Hartree integrate/restrict')
            
        # Calculate atomic hamiltonians:
        W_aL = {}
        for a in density.D_asp:
            W_aL[a] = np.empty((self.setups[a].lmax + 1)**2)
        density.ghat.integrate(self.vHt_g, W_aL)

        return Ekin, Epot, Ebar, Eext, Exc, W_aL

    def calculate_forces2(self, dens, ghat_aLv, nct_av, vbar_av):
        if self.nspins == 2:
            vt_G = self.vt_sG.mean(0)
        else:
            vt_G = self.vt_sG[0]

        dens.ghat.derivative(self.vHt_g, ghat_aLv)
        dens.nct.derivative(vt_G, nct_av)
        self.vbar.derivative(dens.nt_g, vbar_av)
