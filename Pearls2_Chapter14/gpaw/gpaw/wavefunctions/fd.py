import numpy as np

from gpaw.kpoint import KPoint
from gpaw.mpi import serial_comm
from gpaw.io import FileReference
from gpaw.utilities.blas import axpy
from gpaw.transformers import Transformer
from gpaw.hs_operators import MatrixOperator
from gpaw.preconditioner import Preconditioner
from gpaw.fd_operators import Laplace, Gradient
from gpaw.kpt_descriptor import KPointDescriptor
from gpaw.wavefunctions.fdpw import FDPWWaveFunctions
from gpaw.lfc import LocalizedFunctionsCollection as LFC
from gpaw import use_mic
from gpaw.mic import stream


class FD:
    name = 'fd'

    def __init__(self):
        pass  # Could hold parameters like FD stencils

    def __call__(self, *args, **kwargs):
        return FDWaveFunctions(*args, **kwargs)


class FDWaveFunctions(FDPWWaveFunctions):
    mode = 'fd'

    def __init__(self, stencil, diagksl, orthoksl, initksl,
                 gd, nvalence, setups, bd,
                 dtype, world, kd, kptband_comm, timer=None):
        FDPWWaveFunctions.__init__(self, diagksl, orthoksl, initksl,
                                   gd, nvalence, setups, bd,
                                   dtype, world, kd, kptband_comm, timer)

        # Kinetic energy operator:
        self.kin = Laplace(self.gd, -0.5, stencil, self.dtype)

        self.matrixoperator = MatrixOperator(self.orthoksl)

        self.taugrad_v = None  # initialized by MGGA functional

    def empty(self, n=(), global_array=False, realspace=False, q=-1):
        return self.gd.empty(n, self.dtype, global_array)

    def integrate(self, a_xg, b_yg=None, global_integral=True):
        return self.gd.integrate(a_xg, b_yg, global_integral)

    def bytes_per_wave_function(self):
        return self.gd.bytecount(self.dtype)

    def set_setups(self, setups):
        self.pt = LFC(self.gd, [setup.pt_j for setup in setups],
                      self.kd, dtype=self.dtype, forces=True)
        FDPWWaveFunctions.set_setups(self, setups)

    def set_positions(self, spos_ac):
        FDPWWaveFunctions.set_positions(self, spos_ac)

    def summary(self, fd):
        fd.write('Wave functions: Uniform real-space grid\n')
        fd.write('Kinetic energy operator: %s\n' % self.kin.description)
        
    def make_preconditioner(self, block=1):
        return Preconditioner(self.gd, self.kin, self.dtype, block)
    
    def apply_pseudo_hamiltonian(self, kpt, hamiltonian, psit_xG, Htpsit_xG):
        self.timer.start('Apply hamiltonian')
        self.kin.apply(psit_xG, Htpsit_xG, kpt.phase_cd)
        hamiltonian.apply_local_potential(psit_xG, Htpsit_xG, kpt.s)
        self.timer.stop('Apply hamiltonian')

    def add_orbital_density(self, nt_G, kpt, n):
        if self.dtype == float:
            axpy(1.0, kpt.psit_nG[n]**2, nt_G)
        else:
            axpy(1.0, kpt.psit_nG[n].real**2, nt_G)
            axpy(1.0, kpt.psit_nG[n].imag**2, nt_G)

    def add_to_density_from_k_point_with_occupation(self, nt_sG, kpt, f_n):
        # Used in calculation of response part of GLLB-potential
        nt_G = nt_sG[kpt.s]
        if self.dtype == float:
            for f, psit_G in zip(f_n, kpt.psit_nG):
                axpy(f, psit_G**2, nt_G)
        else:
            for f, psit_G in zip(f_n, kpt.psit_nG):
                axpy(f, psit_G.real**2, nt_G)
                axpy(f, psit_G.imag**2, nt_G)

        # Hack used in delta-scf calculations:
        if hasattr(kpt, 'c_on'):
            assert self.bd.comm.size == 1
            d_nn = np.zeros((self.bd.mynbands, self.bd.mynbands),
                            dtype=complex)
            for ne, c_n in zip(kpt.ne_o, kpt.c_on):
                d_nn += ne * np.outer(c_n.conj(), c_n)
            for d_n, psi0_G in zip(d_nn, kpt.psit_nG):
                for d, psi_G in zip(d_n, kpt.psit_nG):
                    if abs(d) > 1.e-12:
                        nt_G += (psi0_G.conj() * d * psi_G).real

    def calculate_kinetic_energy_density(self):
        if self.taugrad_v is None:
            self.taugrad_v = [
                Gradient(self.gd, v, n=3, dtype=self.dtype).apply
                for v in range(3)]
            
        assert not hasattr(self.kpt_u[0], 'c_on')
        if self.kpt_u[0].psit_nG is None:
            raise RuntimeError('No wavefunctions yet')
        if isinstance(self.kpt_u[0].psit_nG, FileReference):
            # XXX initialize
            raise RuntimeError('Wavefunctions have not been initialized.')

        taut_sG = self.gd.zeros(self.nspins)
        dpsit_G = self.gd.empty(dtype=self.dtype)
        for kpt in self.kpt_u:
            for f, psit_G in zip(kpt.f_n, kpt.psit_nG):
                for v in range(3):
                    self.taugrad_v[v](psit_G, dpsit_G, kpt.phase_cd)
                    axpy(0.5 * f, abs(dpsit_G)**2, taut_sG[kpt.s])

        self.kd.comm.sum(taut_sG)
        self.band_comm.sum(taut_sG)
        return taut_sG
        
    def apply_mgga_orbital_dependent_hamiltonian(self, kpt, psit_xG,
                                                 Htpsit_xG, dH_asp,
                                                 dedtaut_G):
        a_G = self.gd.empty(dtype=psit_xG.dtype)
        for psit_G, Htpsit_G in zip(psit_xG, Htpsit_xG):
            for v in range(3):
                self.taugrad_v[v](psit_G, a_G, kpt.phase_cd)
                self.taugrad_v[v](dedtaut_G * a_G, a_G, kpt.phase_cd)
                axpy(-0.5, a_G, Htpsit_G)

    def ibz2bz(self, atoms):
        """Transform wave functions in IBZ to the full BZ."""

        assert self.kd.comm.size == 1

        # New k-point descriptor for full BZ:
        kd = KPointDescriptor(self.kd.bzk_kc, nspins=self.nspins)
        #kd.set_symmetry(atoms, self.setups, enabled=False)
        kd.set_communicator(serial_comm)

        self.pt = LFC(self.gd, [setup.pt_j for setup in self.setups],
                      kd, dtype=self.dtype)
        self.pt.set_positions(atoms.get_scaled_positions())

        self.initialize_wave_functions_from_restart_file()

        weight = 2.0 / kd.nspins / kd.nbzkpts
        
        # Build new list of k-points:
        kpt_u = []
        for s in range(self.nspins):
            for k in range(kd.nbzkpts):
                # Index of symmetry related point in the IBZ
                ik = self.kd.bz2ibz_k[k]
                r, u = self.kd.get_rank_and_index(s, ik)
                assert r == 0
                kpt = self.kpt_u[u]
            
                phase_cd = np.exp(2j * np.pi * self.gd.sdisp_cd *
                                  kd.bzk_kc[k, :, np.newaxis])

                # New k-point:
                kpt2 = KPoint(weight, s, k, k, phase_cd)
                kpt2.f_n = kpt.f_n / kpt.weight / kd.nbzkpts * 2 / self.nspins
                kpt2.eps_n = kpt.eps_n.copy()
                
                # Transform wave functions using symmetry operation:
                Psit_nG = self.gd.collect(kpt.psit_nG)
                if Psit_nG is not None:
                    Psit_nG = Psit_nG.copy()
                    for Psit_G in Psit_nG:
                        Psit_G[:] = self.kd.transform_wave_function(Psit_G, k)
                kpt2.psit_nG = self.gd.empty(self.bd.nbands, dtype=self.dtype)
                self.gd.distribute(Psit_nG, kpt2.psit_nG)

                # Calculate PAW projections:
                kpt2.P_ani = self.pt.dict(len(kpt.psit_nG))
                self.pt.integrate(kpt2.psit_nG, kpt2.P_ani, k)
                
                kpt_u.append(kpt2)

        self.kd = kd
        self.kpt_u = kpt_u

    def write(self, writer, write_wave_functions=False):
        writer['Mode'] = 'fd'

        if not write_wave_functions:
            return

        writer.add('PseudoWaveFunctions',
                   ('nspins', 'nibzkpts', 'nbands',
                    'ngptsx', 'ngptsy', 'ngptsz'),
                   dtype=self.dtype)

        if hasattr(writer, 'hdf5'):
            parallel = (self.world.size > 1)
            for kpt in self.kpt_u:
                indices = [kpt.s, kpt.k]
                indices.append(self.bd.get_slice())
                indices += self.gd.get_slice()
                writer.fill(kpt.psit_nG, parallel=parallel, *indices)
        else:
            for s in range(self.nspins):
                for k in range(self.kd.nibzkpts):
                    for n in range(self.bd.nbands):
                        psit_G = self.get_wave_function_array(n, k, s)
                        writer.fill(psit_G, s, k, n)

    def read(self, reader, hdf5):
        if ((not hdf5 and self.bd.comm.size == 1) or
            (hdf5 and self.world.size == 1)):
            # We may not be able to keep all the wave
            # functions in memory - so psit_nG will be a special type of
            # array that is really just a reference to a file:
            for kpt in self.kpt_u:
                kpt.psit_nG = reader.get_reference('PseudoWaveFunctions',
                                                   (kpt.s, kpt.k))
        else:
            for kpt in self.kpt_u:
                kpt.psit_nG = self.empty(self.bd.mynbands)
                if hdf5:
                    indices = [kpt.s, kpt.k]
                    indices.append(self.bd.get_slice())
                    indices += self.gd.get_slice()
                    reader.get('PseudoWaveFunctions', out=kpt.psit_nG,
                               parallel=(self.world.size > 1), *indices)
                else:
                    # Read band by band to save memory
                    for myn, psit_G in enumerate(kpt.psit_nG):
                        n = self.bd.global_index(myn)
                        if self.gd.comm.rank == 0:
                            big_psit_G = np.array(
                                reader.get('PseudoWaveFunctions',
                                           kpt.s, kpt.k, n),
                                self.dtype)
                        else:
                            big_psit_G = None
                        self.gd.distribute(big_psit_G, psit_G)
        
    def initialize_from_lcao_coefficients(self, basis_functions, mynbands):
        for kpt in self.kpt_u:
            kpt.psit_nG = self.gd.zeros(self.bd.mynbands, self.dtype)
            basis_functions.lcao_to_grid(kpt.C_nM,
                                         kpt.psit_nG[:mynbands], kpt.q)
            kpt.C_nM = None
            if use_mic:
                kpt.psit_nG_mic = stream.bind(kpt.psit_nG)
                stream.sync()

    def random_wave_functions(self, nao):
        """Generate random wave functions."""

        gpts = self.gd.N_c[0] * self.gd.N_c[1] * self.gd.N_c[2]
        
        if self.bd.nbands < gpts / 64:
            gd1 = self.gd.coarsen()
            gd2 = gd1.coarsen()

            psit_G1 = gd1.empty(dtype=self.dtype)
            psit_G2 = gd2.empty(dtype=self.dtype)

            interpolate2 = Transformer(gd2, gd1, 1, self.dtype).apply
            interpolate1 = Transformer(gd1, self.gd, 1, self.dtype).apply

            shape = tuple(gd2.n_c)
            scale = np.sqrt(12 / abs(np.linalg.det(gd2.cell_cv)))

            old_state = np.random.get_state()

            np.random.seed(4 + self.world.rank)

            for kpt in self.kpt_u:
                for psit_G in kpt.psit_nG[nao:]:
                    if self.dtype == float:
                        psit_G2[:] = (np.random.random(shape) - 0.5) * scale
                    else:
                        psit_G2.real = (np.random.random(shape) - 0.5) * scale
                        psit_G2.imag = (np.random.random(shape) - 0.5) * scale

                    interpolate2(psit_G2, psit_G1, kpt.phase_cd)
                    interpolate1(psit_G1, psit_G, kpt.phase_cd)
            np.random.set_state(old_state)
        
        elif gpts / 64 <= self.bd.nbands < gpts / 8:
            gd1 = self.gd.coarsen()

            psit_G1 = gd1.empty(dtype=self.dtype)

            interpolate1 = Transformer(gd1, self.gd, 1, self.dtype).apply

            shape = tuple(gd1.n_c)
            scale = np.sqrt(12 / abs(np.linalg.det(gd1.cell_cv)))

            old_state = np.random.get_state()

            np.random.seed(4 + self.world.rank)

            for kpt in self.kpt_u:
                for psit_G in kpt.psit_nG[nao:]:
                    if self.dtype == float:
                        psit_G1[:] = (np.random.random(shape) - 0.5) * scale
                    else:
                        psit_G1.real = (np.random.random(shape) - 0.5) * scale
                        psit_G1.imag = (np.random.random(shape) - 0.5) * scale

                    interpolate1(psit_G1, psit_G, kpt.phase_cd)
            np.random.set_state(old_state)
               
        else:
            shape = tuple(self.gd.n_c)
            scale = np.sqrt(12 / abs(np.linalg.det(self.gd.cell_cv)))

            old_state = np.random.get_state()

            np.random.seed(4 + self.world.rank)

            for kpt in self.kpt_u:
                for psit_G in kpt.psit_nG[nao:]:
                    if self.dtype == float:
                        psit_G[:] = (np.random.random(shape) - 0.5) * scale
                    else:
                        psit_G.real = (np.random.random(shape) - 0.5) * scale
                        psit_G.imag = (np.random.random(shape) - 0.5) * scale

            np.random.set_state(old_state)

    def estimate_memory(self, mem):
        FDPWWaveFunctions.estimate_memory(self, mem)
