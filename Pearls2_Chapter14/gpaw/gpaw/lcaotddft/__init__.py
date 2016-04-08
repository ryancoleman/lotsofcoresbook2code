from __future__ import print_function
from gpaw import GPAW
from gpaw.external_potential import ConstantElectricField
import numpy as np
from math import sqrt
from gpaw.utilities.blas import gemm
from numpy.linalg import inv
from numpy import dot
from gpaw.mixer import DummyMixer
from math import pi, log
from ase.parallel import paropen  # XXX function refers directly to world
from gpaw.mpi import world  # XXXXXXXX let's get rid of world and
# use correct communicator, i.e., calc.wfs.world or similar. -askhl
from gpaw.tddft.units import attosec_to_autime, autime_to_attosec
from numpy.linalg import solve
from gpaw.xc import XC

from gpaw.utilities.scalapack import pblas_simple_hemm, pblas_simple_gemm, \
                                     scalapack_inverse, scalapack_solve, \
                                     scalapack_zero, pblas_tran, scalapack_set
                                     
from gpaw.utilities.tools import tri2full

import sys
from time import localtime


def print_matrix(M, file=None, rank=0):
    # XXX Debugging stuff. Remove.
    if world.rank == 0:
        if file is not None:
            f = open(file, 'w')
        else:
            f = sys.stdout
        a, b = M.shape
        for i in range(a):
            for j in range(b):
                print('%.7f' % M[i][j].real, '%.7f' % M[i][j].imag, end=' ', file=f)
            print(file=f)
        if file is not None:
            f.close()


def verify(data1, data2, id, uplo='B'):
    # Debugging stuff. Remove
    if uplo == 'B':
        err = sum(abs(data1 - data2).ravel()**2)
    else:
        err = 0
        N, M = data1.shape
        for i in range(N):
            for j in range(i, M):
                if uplo == 'L':
                    if i >= j:
                        err += abs(data1[i][j] - data2[i][j])**2
                if uplo == 'U':
                    if i <= j:
                        err += abs(data1[i][j] - data2[i][j])**2
    if err > 1e-7:
        print('verify err', err)
    if err > 1e-5:
        print('Parallel assert failed: ', id, ' norm: ', err)
        print('Data from proc ', world.rank)
        print('First', data1)
        print('Second', data2)
        print('Diff', data1-data2)
        assert False


def mpiverify(data, id):
    # Do some debugging when running on two procs XXX REMOVE
    if world.size == 2:
        if world.rank == 0:
            temp = -data.copy()
        else:
            temp = data.copy()
        world.sum(temp)
        err = sum(abs(temp).ravel()**2)
        if err > 1e-10:
            if world.rank == 0:
                print('Parallel assert failed: ', id, ' norm: ',
                      sum(temp.ravel()**2))
            print('Data from proc ', world.rank)
            print(data)
            assert False

            
class KickHamiltonian:
    def __init__(self, calc, ext):
        self.ext = ext
        self.vt_sG = [ext.get_potential(gd=calc.density.gd)]
        self.dH_asp = {}

        # This code is copy-paste from hamiltonian.update
        for a, D_sp in calc.density.D_asp.items():
            setup = calc.hamiltonian.setups[a]
            vext = ext.get_taylor(spos_c=calc.hamiltonian.spos_ac[a, :])
            # Taylor expansion to the zeroth order
            self.dH_asp[a] = [vext[0][0] * sqrt(4 * pi) * setup.Delta_pL[:, 0]]
            if len(vext) > 1:
                # Taylor expansion to the first order
                Delta_p1 = np.array([setup.Delta_pL[:, 1],
                                     setup.Delta_pL[:, 2],
                                     setup.Delta_pL[:, 3]])
                self.dH_asp[a] += sqrt(4 * pi / 3) * np.dot(vext[1], Delta_p1)


class LCAOTDDFT(GPAW):
    def __init__(self, filename=None, propagator_debug=False,
                 propagator='cn', fxc=None, **kwargs):
        self.time = 0.0
        self.niter = 0
        self.kick_strength = [0.0, 0.0, 0.0]
        GPAW.__init__(self, filename, **kwargs)
        self.propagator_debug = propagator_debug
        self.tddft_initialized = False
        self.fxc = fxc
        self.propagator = propagator

        # Restarting from a file
        if filename is not None:
            self.initialize()
            self.set_positions()

    def propagate_wfs(self, sourceC_nm, targetC_nm, S_MM, H_MM, dt):
        if self.propagator == 'cn':
            return self.linear_propagator(sourceC_nm, targetC_nm, S_MM, H_MM, dt)
        raise NotImplementedError

    def linear_propagator(self, sourceC_nM, targetC_nM, S_MM, H_MM, dt):
        self.timer.start('Linear solve')
        # XXX Debugging stuff. Remove
        if self.propagator_debug:
            if self.blacs:
                globalH_MM = self.blacs_mm_to_global(H_MM)
                globalS_MM = self.blacs_mm_to_global(S_MM)
                if world.rank == 0:
                    tri2full(globalS_MM, 'L')
                    tri2full(globalH_MM, 'L')
                    U_MM = dot(inv(globalS_MM-0.5j*globalH_MM*dt), globalS_MM+0.5j*globalH_MM*dt)
                    debugC_nM = dot(sourceC_nM, U_MM.T.conjugate())
                    #print 'PASS PROPAGATOR'
                    #debugC_nM = sourceC_nM.copy()
            else:
                if world.rank == 0:
                    U_MM = dot(inv(S_MM-0.5j*H_MM*dt), S_MM+0.5j*H_MM*dt)
                    debugC_nM = dot(sourceC_nM, U_MM.T.conjugate())
                #print 'PASS PROPAGATOR'
                #debugC_nM = sourceC_nM.copy()

        if self.blacs:
            target_blockC_nm = self.Cnm_block_descriptor.empty(dtype=complex) # XXX, Preallocate
            temp_blockC_nm = self.Cnm_block_descriptor.empty(dtype=complex) # XXX, Preallocate
            temp_block_mm = self.mm_block_descriptor.empty(dtype=complex)
            if self.density.gd.comm.rank != 0:
                # XXX Fake blacks nbands, nao, nbands, nao grid because some weird asserts
                # (these are 0,x or x,0 arrays)
                sourceC_nM = self.CnM_unique_descriptor.zeros(dtype=complex)

            # 1. target = (S+0.5j*H*dt) * source
            # Wave functions to target
            self.CnM2nm.redistribute(sourceC_nM, temp_blockC_nm)

            # XXX It can't be this f'n hard to symmetrize a matrix (tri2full)
            scalapack_zero(self.mm_block_descriptor, H_MM, 'U') # Remove upper diagonal
            temp_block_mm[:] = S_MM - (0.5j*dt) * H_MM  # Lower diagonal matrix
            scalapack_set(self.mm_block_descriptor, temp_block_mm, 0, 0, 'U')
            # Note it's stricly lower diagonal matrix
            pblas_tran(-0.5j*dt, H_MM, 1.0, temp_block_mm, self.mm_block_descriptor, self.mm_block_descriptor) # Add transpose of H
            pblas_tran(1.0, S_MM, 1.0, temp_block_mm, self.mm_block_descriptor, self.mm_block_descriptor) # Add transpose of S

            pblas_simple_gemm(self.Cnm_block_descriptor,
                              self.mm_block_descriptor,
                              self.Cnm_block_descriptor,
                              temp_blockC_nm,
                              temp_block_mm,
                              target_blockC_nm)
            # 2. target = (S-0.5j*H*dt)^-1 * target
            #temp_block_mm[:] = S_MM + (0.5j*dt) * H_MM
            # XXX It can't be this f'n hard to symmetrize a matrix (tri2full)
            temp_block_mm[:] = S_MM + (0.5j*dt) * H_MM  # Lower diagonal matrix
            scalapack_set(self.mm_block_descriptor, temp_block_mm, 0, 0, 'U') # Not it's stricly lower diagonal matrix           
            pblas_tran(+0.5j*dt, H_MM, 1.0, temp_block_mm, self.mm_block_descriptor, self.mm_block_descriptor) # Add transpose of H
            pblas_tran(1.0, S_MM, 1.0, temp_block_mm, self.mm_block_descriptor, self.mm_block_descriptor) # Add transpose of S

            scalapack_solve(self.mm_block_descriptor, 
                            self.Cnm_block_descriptor, 
                            temp_block_mm,
                            target_blockC_nm)

            if self.density.gd.comm.rank != 0: # XXX is this correct?
                # XXX Fake blacks nbands, nao, nbands, nao grid because some weird asserts
                # (these are 0,x or x,0 arrays)
                target = self.CnM_unique_descriptor.zeros(dtype=complex)
            else:
                target = targetC_nM
            self.Cnm2nM.redistribute(target_blockC_nm, target)
            self.density.gd.comm.broadcast(targetC_nM, 0) # Is this required?
        else:
            # Note: The full equation is conjugated (therefore -+, not +-)
            targetC_nM[:] = solve(S_MM-0.5j*H_MM*dt, np.dot(S_MM+0.5j*H_MM*dt, sourceC_nM.T.conjugate())).T.conjugate()
        
        # XXX Debugging stuff. Remove
        if self.propagator_debug:
            if world.rank == 0:
                verify(targetC_nM, debugC_nM,
                       'Linear solver propagator vs. reference')

        self.timer.stop('Linear solve')

    def taylor_propagator(self, sourceC_nM, targetC_nM, S_MM, H_MM, dt):
        self.timer.start('Taylor propagator')
        # XXX Debugging stuff. Remove
        if self.propagator_debug:
            if self.blacs:
                globalH_MM = self.blacs_mm_to_global(H_MM)
                globalS_MM = self.blacs_mm_to_global(S_MM) 
                if world.rank == 0:
                    tri2full(globalS_MM, 'L')
                    tri2full(globalH_MM, 'L')
                    U_MM = dot(inv(globalS_MM-0.5j*globalH_MM*dt), globalS_MM+0.5j*globalH_MM*dt)
                    debugC_nM = dot(sourceC_nM, U_MM.T.conjugate())
                    #print 'PASS PROPAGATOR'
                    #debugC_nM = sourceC_nM.copy()
            else:
                if world.rank == 0:
                    U_MM = dot(inv(S_MM - 0.5j * H_MM * dt),
                               S_MM + 0.5j * H_MM * dt)
                    debugC_nM = dot(sourceC_nM, U_MM.T.conjugate())
                #print 'PASS PROPAGATOR'
                #debugC_nM = sourceC_nM.copy()

        if self.blacs:
            target_blockC_nm = self.Cnm_block_descriptor.empty(dtype=complex) # XXX, Preallocate
            if self.density.gd.comm.rank != 0: 
                # XXX Fake blacks nbands, nao, nbands, nao grid because some weird asserts
                # (these are 0,x or x,0 arrays)
                sourceC_nM = self.CnM_unique_descriptor.zeros(dtype=complex)

            # Zeroth order taylor to target
            self.CnM2nm.redistribute(sourceC_nM, target_blockC_nm) 

            # XXX, preallocate, optimize use of temporal arrays
            temp_blockC_nm = target_blockC_nm.copy()
            temp2_blockC_nm = target_blockC_nm.copy()

            order = 4
            assert self.wfs.kd.comm.size == 1
            for n in range(order):
                # Multiply with hamiltonian
                pblas_simple_hemm(self.mm_block_descriptor, 
                                  self.Cnm_block_descriptor, 
                                  self.Cnm_block_descriptor, 
                                  H_MM, 
                                  temp_blockC_nm, 
                                  temp2_blockC_nm, side='R') 
                # XXX: replace with not simple gemm
                temp2_blockC_nm *= -1j*dt/(n+1) 
                # Multiply with inverse overlap
                pblas_simple_hemm(self.mm_block_descriptor, 
                                  self.Cnm_block_descriptor,
                                  self.Cnm_block_descriptor, 
                                  self.wfs.kpt_u[0].invS_MM, # XXX
                                  temp2_blockC_nm, 
                                  temp_blockC_nm, side='R')
                target_blockC_nm += temp_blockC_nm
            if self.density.gd.comm.rank != 0: # Todo: Change to gd.rank
                # XXX Fake blacks nbands, nao, nbands, nao grid because some weird asserts
                # (these are 0,x or x,0 arrays)
                target = self.CnM_unique_descriptor.zeros(dtype=complex)
            else:
                target = targetC_nM
            self.Cnm2nM.redistribute(target_blockC_nm, target)

            self.density.gd.comm.broadcast(targetC_nM, 0)
        else:
            assert self.wfs.kd.comm.size == 1
            if self.density.gd.comm.rank == 0:
                targetC_nM[:] = sourceC_nM[:]
                tempC_nM = sourceC_nM.copy()
                order = 4
                for n in range(order):
                    tempC_nM[:] = np.dot(self.wfs.kpt_u[0].invS, np.dot(H_MM, 1j*dt/(n+1)*tempC_nM.T.conjugate())).T.conjugate()
                    targetC_nM += tempC_nM
            self.density.gd.comm.broadcast(targetC_nM, 0)
                
        if self.propagator_debug:
            if world.rank == 0:
                verify(targetC_nM, debugC_nM,
                       'Linear solver propagator vs. reference')

        self.timer.stop('Taylor propagator')

    def kick(self, strength):
        self.tddft_init()
        self.timer.start('Kick')
        self.kick_strength = strength

        # magnitude
        magnitude = np.sqrt(strength[0]*strength[0] 
                             + strength[1]*strength[1] 
                             + strength[2]*strength[2])

        # normalize
        direction = strength / magnitude

        self.text('Applying absorbtion kick')
        self.text('Magnitude: %.8f ' % magnitude)
        self.text('Direction: %.4f %.4f %.4f' % tuple(direction))

        # Create hamiltonian object for absorbtion kick
        kick_hamiltonian = KickHamiltonian(self, ConstantElectricField(magnitude, direction=direction))
        for k, kpt in enumerate(self.wfs.kpt_u):
            Vkick_MM = self.wfs.eigensolver.calculate_hamiltonian_matrix(kick_hamiltonian, self.wfs, kpt, add_kinetic=False, root=-1)
            for i in range(10):
                self.propagate_wfs(kpt.C_nM, kpt.C_nM, kpt.S_MM, Vkick_MM, 0.1)
        self.timer.stop('Kick')

    def blacs_mm_to_global(self, H_mm):
        target = self.MM_descriptor.empty(dtype=complex)
        self.mm2MM.redistribute(H_mm, target)
        world.barrier()
        return target

    def blacs_nm_to_global(self, C_nm):
        target = self.CnM_unique_descriptor.empty(dtype=complex)
        self.Cnm2nM.redistribute(C_nm, target)
        world.barrier()
        return target

    def tddft_init(self):
        if not self.tddft_initialized:
            if world.rank == 0:
                print('Initializing real time LCAO TD-DFT calculation.')
                print('XXX Warning: Array use not optimal for memory.')
                print('XXX Taylor propagator probably doesn\'t work')
                print('XXX ...and no arrays are listed in memory estimate yet.')
            self.blacs = self.wfs.ksl.using_blacs
            if self.blacs:
                self.ksl = ksl = self.wfs.ksl    
                nao = ksl.nao
                nbands = ksl.bd.nbands
                mynbands = ksl.bd.mynbands
                blocksize = ksl.blocksize

                from gpaw.blacs import Redistributor
                if world.rank == 0:
                    print('BLACS Parallelization')

                # Parallel grid descriptors
                self.MM_descriptor = ksl.blockgrid.new_descriptor(nao, nao, nao, nao) # FOR DEBUG
                self.mm_block_descriptor = ksl.blockgrid.new_descriptor(nao, nao, blocksize, blocksize)
                self.Cnm_block_descriptor = ksl.blockgrid.new_descriptor(nbands, nao, blocksize, blocksize)
                #self.CnM_descriptor = ksl.blockgrid.new_descriptor(nbands, nao, mynbands, nao)
                self.mM_column_descriptor = ksl.single_column_grid.new_descriptor(nao, nao, ksl.naoblocksize, nao)
                self.CnM_unique_descriptor = ksl.single_column_grid.new_descriptor(nbands, nao, mynbands, nao)

                # Redistributors
                self.mm2MM = Redistributor(ksl.block_comm,
                                           self.mm_block_descriptor,
                                           self.MM_descriptor) # XXX FOR DEBUG
                self.MM2mm = Redistributor(ksl.block_comm,
                                           self.MM_descriptor,
                                           self.mm_block_descriptor) # XXX FOR DEBUG
                self.Cnm2nM = Redistributor(ksl.block_comm,
                                            self.Cnm_block_descriptor,
                                            self.CnM_unique_descriptor) 
                self.CnM2nm = Redistributor(ksl.block_comm,
                                            self.CnM_unique_descriptor,
                                            self.Cnm_block_descriptor) 
                self.mM2mm =  Redistributor(ksl.block_comm,
                                            self.mM_column_descriptor,
                                            self.mm_block_descriptor)

                for kpt in self.wfs.kpt_u:
                    scalapack_zero(self.mm_block_descriptor, kpt.S_MM,'U')
                    scalapack_zero(self.mm_block_descriptor, kpt.T_MM,'U')

                # XXX to propagator class
                if self.propagator == 'taylor' and self.blacs:  
                    # cholS_mm = self.mm_block_descriptor.empty(dtype=complex)
                    for kpt in self.wfs.kpt_u:
                        kpt.invS_MM = kpt.S_MM.copy()
                        scalapack_inverse(self.mm_block_descriptor,
                                          kpt.invS_MM, 'L')
                    if self.propagator_debug:
                        if world.rank == 0:
                            print('XXX Doing serial inversion of overlap matrix.')
                        self.timer.start('Invert overlap (serial)')
                        invS2_MM = self.MM_descriptor.empty(dtype=complex)
                        for kpt in self.wfs.kpt_u:
                            #kpt.S_MM[:] = 128.0*(2**world.rank)
                            self.mm2MM.redistribute(self.wfs.S_qMM[kpt.q], invS2_MM)
                            world.barrier()
                            if world.rank == 0:
                                tri2full(invS2_MM,'L')
                                invS2_MM[:] = inv(invS2_MM.copy())
                                self.invS2_MM = invS2_MM
                            kpt.invS2_MM = ksl.mmdescriptor.empty(dtype=complex)
                            self.MM2mm.redistribute(invS2_MM, kpt.invS2_MM)
                            verify(kpt.invS_MM, kpt.invS2_MM, 'overlap par. vs. serial', 'L')
                        self.timer.stop('Invert overlap (serial)')
                        if world.rank == 0:
                            print('XXX Overlap inverted.')
                if self.propagator == 'taylor' and not self.blacs:
                    tmp = inv(self.wfs.kpt_u[0].S_MM)
                    self.wfs.kpt_u[0].invS = tmp

            # Reset the density mixer
            self.density.mixer = DummyMixer()    
            self.tddft_initialized = True
            for k, kpt in enumerate(self.wfs.kpt_u):
                kpt.C2_nM = kpt.C_nM.copy()
                #kpt.firstC_nM = kpt.C_nM.copy()

    def update_projectors(self):
        self.timer.start('LCAO update projectors') 
        # Loop over all k-points
        for k, kpt in enumerate(self.wfs.kpt_u):
            for a, P_ni in kpt.P_ani.items():
                print('Update projector: Rank:', world.rank, 'a', a)
                P_ni.fill(117)
                gemm(1.0, kpt.P_aMi[a], kpt.C_nM, 0.0, P_ni, 'n')
        self.timer.stop('LCAO update projectors') 

    def save_wfs(self):
        for k, kpt in enumerate(self.wfs.kpt_u):
            kpt.C2_nM[:] = kpt.C_nM

    def update_hamiltonian(self):
        self.update_projectors()
        self.density.update(self.wfs)
        self.hamiltonian.update(self.density)

    def propagate(self, time_step=10, iterations=2000, out='lcao.dm',
                  dump_interval=50):
        assert self.wfs.dtype == complex
        time_step *= attosec_to_autime
        self.time_step = time_step
        self.dump_interval = dump_interval
        maxiter = self.niter + iterations

        if self.time < self.time_step:
            self.dm_file = paropen(out,'w') # XXXX
            # Bug: will fail if world != self.wfs.world.  -askhl
            header = '# Kick = [%22.12le, %22.12le, %22.12le]\n' \
                   % (self.kick_strength[0], self.kick_strength[1], \
                      self.kick_strength[2])
            header += '# %15s %15s %22s %22s %22s\n' \
                   % ('time', 'norm', 'dmx', 'dmy', 'dmz')
            self.dm_file.write(header)
            self.dm_file.flush()
            self.text('About to do %d propagation steps.' % iterations)
        else:
            self.dm_file = paropen(out,'a') # XXXX
            self.text('About to continue from iteration %d and do %d propagation steps' % (self.niter, maxiter)) 
        self.tddft_init()

        dm0 = None # Initial dipole moment
        self.timer.start('Propagate')
        while self.niter < maxiter:
            dm = self.density.finegd.calculate_dipole_moment(self.density.rhot_g)
            if dm0 is None:
                dm0 = dm
            norm = self.density.finegd.integrate(self.density.rhot_g)
            line = '%20.8lf %20.8le %22.12le %22.12le %22.12le' % (self.time, norm, dm[0], dm[1], dm[2])
            T = localtime()
            if world.rank == 0:
                print(line, file=self.dm_file)

            if world.rank == 0 and self.niter%10==0:
                print('iter: %3d  %02d:%02d:%02d %11.2f   %9.1f %12.8f' % (self.niter,
                                                                           T[3], T[4], T[5],
                                                                           self.time * autime_to_attosec,
                                                                           log(abs(norm)+1e-16)/log(10),
                                                                           np.sqrt(dm[0]**2+dm[1]**2+dm[2]**2)))
                self.dm_file.flush()

            # ---------------------------------------------------------------------------- 
            # Predictor step
            # ----------------------------------------------------------------------------
            # 1. Calculate H(t)
            self.save_wfs() # kpt.C2_nM = kpt.C_nM
            # 2. H_MM(t) = <M|H(t)|H>
            #    Solve Psi(t+dt) from (S_MM - 0.5j*H_MM(t)*dt) Psi(t+dt) = (S_MM + 0.5j*H_MM(t)*dt) Psi(t)

            for k, kpt in enumerate(self.wfs.kpt_u):
                if self.fxc is not None:
                    if self.time == 0.0:
                        kpt.deltaXC_H_MM = self.wfs.eigensolver.calculate_hamiltonian_matrix(\
                            self.hamiltonian, self.wfs, kpt, root=-1)
                        self.hamiltonian.xc = XC(self.fxc)
                        self.update_hamiltonian()
                        assert len(self.wfs.kpt_u) == 1
                        kpt.deltaXC_H_MM -= self.wfs.eigensolver.calculate_hamiltonian_matrix(\
                            self.hamiltonian, self.wfs, kpt, root=-1)

            self.update_hamiltonian()

            for k, kpt in enumerate(self.wfs.kpt_u):
                kpt.H0_MM = self.wfs.eigensolver.calculate_hamiltonian_matrix(self.hamiltonian, self.wfs, kpt, root=-1)
                if self.fxc is not None:
                    kpt.H0_MM += kpt.deltaXC_H_MM
                self.propagate_wfs(kpt.C_nM, kpt.C_nM, kpt.S_MM, kpt.H0_MM, self.time_step)
            # ----------------------------------------------------------------------------
            # Propagator step
            # ----------------------------------------------------------------------------
            # 1. Calculate H(t+dt)
            self.update_hamiltonian()
            # 2. Estimate H(t+0.5*dt) ~ H(t) + H(t+dT)
            for k, kpt in enumerate(self.wfs.kpt_u):
                kpt.H0_MM *= 0.5
                if self.fxc is not None:
                    #  Store this to H0_MM and maybe save one extra H_MM of memory?
                    kpt.H0_MM += 0.5 * (self.wfs.eigensolver.calculate_hamiltonian_matrix( \
                                             self.hamiltonian, self.wfs, kpt, root=-1) +\
                                             kpt.deltaXC_H_MM)
                else:
                    #  Store this to H0_MM and maybe save one extra H_MM of memory?
                    kpt.H0_MM += 0.5 * self.wfs.eigensolver.calculate_hamiltonian_matrix( \
                                             self.hamiltonian, self.wfs, kpt, root=-1)

                # 3. Solve Psi(t+dt) from (S_MM - 0.5j*H_MM(t+0.5*dt)*dt) Psi(t+dt) = (S_MM + 0.5j*H_MM(t+0.5*dt)*dt) Psi(t)
                self.propagate_wfs(kpt.C2_nM, kpt.C_nM, kpt.S_MM, kpt.H0_MM, self.time_step)

            self.niter += 1
            self.time += self.time_step
            
            # Call registered callback functions
            self.call_observers(self.niter)

        self.call_observers(self.niter, final=True)
        self.dm_file.close()
        self.timer.stop('Propagate')
