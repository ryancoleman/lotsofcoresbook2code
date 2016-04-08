# Copyright (C) 2010  CAMd
# Copyright (C) 2010  Argonne National Laboratory
# Please see the accompanying LICENSE file for further information.
import numpy as np

from gpaw.mpi import SerialCommunicator, serial_comm
from gpaw.matrix_descriptor import MatrixDescriptor, \
    BandMatrixDescriptor, \
    BlacsBandMatrixDescriptor
from blacs import BlacsGrid, Redistributor
from gpaw.utilities import uncamelcase
from gpaw.utilities.blas import gemm, r2k, gemmdot
from gpaw.utilities.lapack import diagonalize, general_diagonalize, \
    inverse_cholesky
from gpaw.utilities.scalapack import pblas_simple_gemm, pblas_tran
from gpaw.utilities.tools import tri2full
from gpaw.utilities.timing import nulltimer


def get_KohnSham_layouts(sl, mode, gd, bd, block_comm, dtype, **kwargs):
    """Create Kohn-Sham layouts object."""
    # Not needed for AtomPAW special mode, as usual we just provide whatever
    # happens to make the code not crash
    if not isinstance(mode, str):
        return None  #XXX
    name = {'fd': 'BandLayouts', 'lcao': 'OrbitalLayouts'}[mode]
    args = (gd, bd, block_comm, dtype)
    if sl is not None:
        name = 'Blacs' + name
        assert len(sl) == 3
        args += tuple(sl)
    ksl = {'BandLayouts':         BandLayouts,
           'BlacsBandLayouts':    BlacsBandLayouts,
           'BlacsOrbitalLayouts': BlacsOrbitalLayouts,
           'OrbitalLayouts':      OrbitalLayouts,
            }[name](*args, **kwargs)
    if 0: #XXX debug
        print(('USING KSL: %s' % repr(ksl)))
    assert isinstance(ksl, KohnShamLayouts)
    assert isinstance(ksl, BlacsLayouts) == (sl is not None)
    return ksl


class KohnShamLayouts:
    using_blacs = False  # This is only used by a regression test
    matrix_descriptor_class = None

    def __init__(self, gd, bd, block_comm, dtype, timer=nulltimer):
        assert gd.comm.parent is bd.comm.parent  # must have same parent comm
        self.world = bd.comm.parent
        self.gd = gd
        self.bd = bd
        self.dtype = dtype
        self.block_comm = block_comm
        self.timer = timer
        self._kwargs = {'timer': timer}

        if gd.comm.rank == 0:
            self.column_comm = bd.comm
        else:
            self.column_comm = None

    def get_keywords(self):
        return self._kwargs.copy()  # just a shallow copy...

    def diagonalize(self, *args, **kwargs):
        raise RuntimeError('Virtual member function should not be called.')

    def inverse_cholesky(self, *args, **kwargs):
        raise RuntimeError('Virtual member function should not be called.')

    def new_descriptor(self):
        return self.matrix_descriptor_class(self.bd, self.gd, self)

    def __repr__(self):
        return uncamelcase(self.__class__.__name__)

    def get_description(self):
        """Description of this object in prose, e.g. for logging.

        Subclasses are expected to override this with something useful."""
        return repr(self)


class BlacsLayouts(KohnShamLayouts):
    using_blacs = True  # This is only used by a regression test

    def __init__(self, gd, bd, block_comm, dtype, mcpus, ncpus,
                 blocksize, timer=nulltimer):
        KohnShamLayouts.__init__(self, gd, bd, block_comm, dtype,
                                 timer)
        # WARNING: Do not create the BlacsGrid on a communicator which does not
        # contain block_comm.rank = 0. This will break BlacsBandLayouts which
        # assume eps_M will be broadcast over block_comm.
        self.blocksize = blocksize
        self.blockgrid = BlacsGrid(self.block_comm, mcpus, ncpus)

    def get_description(self):
        title = 'BLACS'
        template = '%d x %d grid with %d x %d blocksize'
        return (title, template)


class BandLayouts(KohnShamLayouts):
    matrix_descriptor_class = BandMatrixDescriptor

    def __init__(self, gd, bd, block_comm, dtype,
                 buffer_size=None, timer=nulltimer):
        KohnShamLayouts.__init__(self, gd, bd, block_comm, dtype,
                                 timer)
        self.buffer_size = buffer_size

    def diagonalize(self, H_NN, eps_n):
        """Serial diagonalizer must handle two cases:
        1. Parallelization over domains only.
        2. Simultaneous parallelization over domains and bands.
        """
        nbands = self.bd.nbands
        mynbands = self.bd.mynbands
        eps_N = np.empty(nbands)
        self.timer.start('Diagonalize')
        # Broadcast on block_comm since result
        # is k-point and spin-dependent only
        self.block_comm.broadcast(H_NN, 0)
        # The result on different processor is not necessarily bit-wise
        # identical, so only domain master performs diagonalization
        if self.gd.comm.rank == 0:
            self._diagonalize(H_NN, eps_N)
        self.gd.comm.broadcast(H_NN, 0)
        self.gd.comm.broadcast(eps_N, 0)
        self.timer.stop('Diagonalize')

        self.timer.start('Distribute results')
        # Copy the portion that belongs to my band group
        eps_n[:] = eps_N[self.bd.get_slice()]
        self.timer.stop('Distribute results')

    def _diagonalize(self, H_NN, eps_N):
        """Serial diagonalize via LAPACK."""
        # This is replicated computation but ultimately avoids
        # additional communication.
        diagonalize(H_NN, eps_N)

    def inverse_cholesky(self, S_NN):
        """Serial inverse Cholesky must handle two cases:
        1. Parallelization over domains only.
        2. Simultaneous parallelization over domains and bands.
        """
        self.timer.start('Inverse Cholesky')
        # Broadcast on block_comm since result
        # is k-point and spin-dependent only
        self.block_comm.broadcast(S_NN, 0)
        # The result on different processor is not necessarily bit-wise
        # identical, so only domain master performs computation
        if self.gd.comm.rank == 0:
            self._inverse_cholesky(S_NN)
        self.gd.comm.broadcast(S_NN, 0)
        self.timer.stop('Inverse Cholesky')

    def _inverse_cholesky(self, S_NN):
        """Serial inverse Cholesky via LAPACK."""
        # This is replicated computation but ultimately avoids
        # additional communication.
        inverse_cholesky(S_NN)

    def get_description(self):
        return 'Serial LAPACK'


class BlacsBandLayouts(BlacsLayouts):  #XXX should derive from BandLayouts too!
    """ScaLAPACK Dense Linear Algebra.

    This class is instantiated in the real-space code.  Not for
    casual use, at least for now.
    
    Requires two distributors and three descriptors for initialization
    as well as grid descriptors and band descriptors. Distributors are
    for cols2blocks (1D -> 2D BLACS grid) and blocks2rows (2D -> 1D
    BLACS grid). ScaLAPACK operations must occur on a 2D BLACS grid for
    performance and scalability. Redistribute of 1D *column* layout
    matrix will operate only on lower half of H or S. Redistribute of
    2D block will operate on entire matrix for U, but only lower half
    of C.

    inverse_cholesky is "hard-coded" for real-space code.
    Expects overlap matrix (S) and the coefficient matrix (C) to be a
    replicated data structures and *not* created by the BLACS descriptor class.
    This is due to the MPI_Reduce and MPI_Broadcast that will occur
    in the parallel matrix multiply. Input matrices should be:
    S = np.empty((nbands, mybands), dtype)
    C = np.empty((mybands, nbands), dtype)

    
    _standard_diagonalize is "hard-coded" for the real-space code.
    Expects both hamiltonian (H) and eigenvector matrix (U) to be a
    replicated data structures and not created by the BLACS descriptor class.
    This is due to the MPI_Reduce and MPI_Broadcast that will occur
    in the parallel matrix multiply. Input matrices should be:
    H = np.empty((nbands, mynbands), dtype)
    U = np.empty((mynbands, nbands), dtype)
    eps_n = np.empty(mynbands, dtype = float)
    """ #XXX rewrite this docstring a bit!

    matrix_descriptor_class = BlacsBandMatrixDescriptor

    # This class 'describes' all the realspace Blacs-related layouts
    def __init__(self, gd, bd, block_comm, dtype, mcpus, ncpus,
                 blocksize, buffer_size=None, timer=nulltimer):
        BlacsLayouts.__init__(self, gd, bd, block_comm, dtype,
                              mcpus, ncpus, blocksize, timer)
        self.buffer_size = buffer_size
        nbands = bd.nbands
        self.mynbands = mynbands = bd.mynbands
        self.blocksize = blocksize

        # 1D layout - columns
        self.columngrid = BlacsGrid(self.column_comm, 1, bd.comm.size)
        self.Nndescriptor = self.columngrid.new_descriptor(nbands, nbands,
                                                           nbands, mynbands)

        # 2D layout
        self.nndescriptor = self.blockgrid.new_descriptor(nbands, nbands,
                                                          blocksize, blocksize)

        # 1D layout - rows
        self.rowgrid = BlacsGrid(self.column_comm, bd.comm.size, 1)
        self.nNdescriptor = self.rowgrid.new_descriptor(nbands, nbands,
                                                        mynbands, nbands)

        # Only redistribute filled out half for Hermitian matrices
        self.Nn2nn = Redistributor(self.block_comm, self.Nndescriptor,
                                   self.nndescriptor)
        #self.Nn2nn = Redistributor(self.block_comm, self.Nndescriptor,
        #                           self.nndescriptor, 'L') #XXX faster but...

        # Resulting matrix will be used in dgemm which is symmetry obvlious
        self.nn2nN = Redistributor(self.block_comm, self.nndescriptor,
                                   self.nNdescriptor)
        
    def diagonalize(self, H_nn, eps_n):
        nbands = self.bd.nbands
        eps_N = np.empty(nbands)
        self.timer.start('Diagonalize')
        self._diagonalize(H_nn, eps_N)
        self.timer.stop('Diagonalize')

        self.timer.start('Distribute results')
        # eps_N is already on block_comm.rank = 0
        # easier to broadcast eps_N to all and
        # get the correct slice afterward.
        self.block_comm.broadcast(eps_N, 0)
        eps_n[:] = eps_N[self.bd.get_slice()]
        self.timer.stop('Distribute results')

    def _diagonalize(self, H_nn, eps_N):
        """Parallel diagonalizer."""
        self.nndescriptor.diagonalize_dc(H_nn.copy(), H_nn, eps_N, 'L')
        
    def inverse_cholesky(self, S_nn):
        self.timer.start('Inverse Cholesky')
        self._inverse_cholesky(S_nn)
        self.block_comm.barrier() # removing barrier may lead to race condition
        self.timer.stop('Inverse Cholesky')
        
    def _inverse_cholesky(self, S_nn):
        self.nndescriptor.inverse_cholesky(S_nn, 'L')

    def get_description(self):
        (title, template) = BlacsLayouts.get_description(self)
        bg = self.blockgrid
        desc = self.nndescriptor
        s = template % (bg.nprow, bg.npcol, desc.mb, desc.nb)
        return ' '.join([title, s])


class BlacsOrbitalLayouts(BlacsLayouts):
    """ScaLAPACK Dense Linear Algebra.

    This class is instantiated in LCAO.  Not for casual use, at least for now.
    
    Requires two distributors and three descriptors for initialization
    as well as grid descriptors and band descriptors. Distributors are
    for cols2blocks (1D -> 2D BLACS grid) and blocks2cols (2D -> 1D
    BLACS grid). ScaLAPACK operations must occur on 2D BLACS grid for
    performance and scalability.

    _general_diagonalize is "hard-coded" for LCAO.
    Expects both Hamiltonian and Overlap matrix to be on the 2D BLACS grid.
    This is done early on to save memory.
    """
    # XXX rewrite this docstring a bit!

    # This class 'describes' all the LCAO Blacs-related layouts
    def __init__(self, gd, bd, block_comm, dtype, mcpus, ncpus,
                 blocksize, nao, timer=nulltimer):
        BlacsLayouts.__init__(self, gd, bd, block_comm, dtype,
                              mcpus, ncpus, blocksize, timer)
        nbands = bd.nbands
        self.blocksize = blocksize
        self.mynbands = mynbands = bd.mynbands
        
        self.orbital_comm = self.bd.comm
        self.naoblocksize = naoblocksize = -((-nao) // self.orbital_comm.size)
        self.nao = nao

        # Range of basis functions for BLACS distribution of matrices:
        self.Mmax = nao
        self.Mstart = bd.comm.rank * naoblocksize
        self.Mstop = min(self.Mstart + naoblocksize, self.Mmax)
        self.mynao = self.Mstop - self.Mstart

        # Column layout for one matrix per band rank:
        self.columngrid = BlacsGrid(bd.comm, bd.comm.size, 1)
        self.mMdescriptor = self.columngrid.new_descriptor(nao, nao,
                                                           naoblocksize, nao)
        self.nMdescriptor = self.columngrid.new_descriptor(nbands, nao,
                                                           mynbands, nao)

        #parallelprint(world, (mynao, self.mMdescriptor.shape))

        # Column layout for one matrix in total (only on grid masters):
        self.single_column_grid = BlacsGrid(self.column_comm, bd.comm.size, 1)
        self.mM_unique_descriptor = self.single_column_grid.new_descriptor( \
            nao, nao, naoblocksize, nao)

        # nM_unique_descriptor is meant to hold the coefficients after
        # diagonalization.  BLACS requires it to be nao-by-nao, but
        # we only fill meaningful data into the first nbands columns.
        #
        # The array will then be trimmed and broadcast across
        # the grid descriptor's communicator.
        self.nM_unique_descriptor = self.single_column_grid.new_descriptor( \
            nbands, nao, mynbands, nao)

        # Fully blocked grid for diagonalization with many CPUs:
        self.mmdescriptor = self.blockgrid.new_descriptor(nao, nao, blocksize,
                                                          blocksize)

        #self.nMdescriptor = nMdescriptor
        self.mM2mm = Redistributor(self.block_comm, self.mM_unique_descriptor,
                                   self.mmdescriptor)
        self.mm2nM = Redistributor(self.block_comm, self.mmdescriptor,
                                   self.nM_unique_descriptor)

    def diagonalize(self, H_mm, C_nM, eps_n, S_mm):
        # C_nM needs to be simultaneously compatible with:
        # 1. outdescriptor
        # 2. broadcast with gd.comm
        # We will does this with a dummy buffer C2_nM
        indescriptor = self.mM2mm.srcdescriptor  # cols2blocks
        outdescriptor = self.mm2nM.dstdescriptor  # blocks2cols
        blockdescriptor = self.mM2mm.dstdescriptor  # cols2blocks

        dtype = S_mm.dtype
        eps_M = np.empty(C_nM.shape[-1])  # empty helps us debug
        subM, subN = outdescriptor.gshape
        
        C_mm = blockdescriptor.zeros(dtype=dtype)
        self.timer.start('General diagonalize')
        # general_diagonalize_ex may have a buffer overflow, so
        # we no longer use it
        #blockdescriptor.general_diagonalize_ex(H_mm, S_mm.copy(), C_mm, eps_M,
        #                                       UL='L', iu=self.bd.nbands)
        blockdescriptor.general_diagonalize_dc(H_mm, S_mm.copy(), C_mm, eps_M,
                                               UL='L')
        self.timer.stop('General diagonalize')
 
        # Make C_nM compatible with the redistributor
        self.timer.start('Redistribute coefs')
        if outdescriptor:
            C2_nM = C_nM
        else:
            C2_nM = outdescriptor.empty(dtype=dtype)
        assert outdescriptor.check(C2_nM)
        self.mm2nM.redistribute(C_mm, C2_nM, subM, subN)  # blocks2cols
        self.timer.stop('Redistribute coefs')

        self.timer.start('Send coefs to domains')
        # eps_M is already on block_comm.rank = 0
        # easier to broadcast eps_M to all and
        # get the correct slice afterward.
        self.block_comm.broadcast(eps_M, 0)
        eps_n[:] = eps_M[self.bd.get_slice()]
        self.gd.comm.broadcast(C_nM, 0)
        self.timer.stop('Send coefs to domains')

    def distribute_overlap_matrix(self, S_qmM, root=0,
                                  add_hermitian_conjugate=False):
        # Some MPI implementations need a lot of memory to do large
        # reductions.  To avoid trouble, we do comm.sum on smaller blocks
        # of S (this code is also safe for arrays smaller than blocksize)
        Sflat_x = S_qmM.ravel()
        blocksize = 2**23 // Sflat_x.itemsize  # 8 MiB
        nblocks = -(-len(Sflat_x) // blocksize)
        Mstart = 0
        for i in range(nblocks):
            self.gd.comm.sum(Sflat_x[Mstart:Mstart + blocksize], root=root)
            Mstart += blocksize
        assert Mstart + blocksize >= len(Sflat_x)

        xshape = S_qmM.shape[:-2]
        nm, nM = S_qmM.shape[-2:]
        S_qmM = S_qmM.reshape(-1, nm, nM)
        
        blockdesc = self.mmdescriptor
        coldesc = self.mM_unique_descriptor
        S_qmm = blockdesc.zeros(len(S_qmM), S_qmM.dtype)

        if not coldesc:  # XXX ugly way to sort out inactive ranks
            S_qmM = coldesc.zeros(len(S_qmM), S_qmM.dtype)
        
        self.timer.start('Distribute overlap matrix')
        for S_mM, S_mm in zip(S_qmM, S_qmm):
            self.mM2mm.redistribute(S_mM, S_mm)
            if add_hermitian_conjugate:
                if blockdesc.active:
                    pblas_tran(1.0, S_mm.copy(), 1.0, S_mm,
                               blockdesc, blockdesc)
                
        self.timer.stop('Distribute overlap matrix')
        return S_qmm.reshape(xshape + blockdesc.shape)

    def get_overlap_matrix_shape(self):
        return self.mmdescriptor.shape

    def calculate_blocked_density_matrix(self, f_n, C_nM):
        nbands = self.bd.nbands
        mynbands = self.bd.mynbands
        nao = self.nao
        dtype = C_nM.dtype
        
        self.nMdescriptor.checkassert(C_nM)
        if self.gd.rank == 0:
            Cf_nM = (C_nM * f_n[:, None]).conj()
        else:
            C_nM = self.nM_unique_descriptor.zeros(dtype=dtype)
            Cf_nM = self.nM_unique_descriptor.zeros(dtype=dtype)

        r = Redistributor(self.block_comm, self.nM_unique_descriptor,
                          self.mmdescriptor)

        Cf_mm = self.mmdescriptor.zeros(dtype=dtype)
        r.redistribute(Cf_nM, Cf_mm, nbands, nao)
        del Cf_nM
        
        C_mm = self.mmdescriptor.zeros(dtype=dtype)
        r.redistribute(C_nM, C_mm, nbands, nao)
        # no use to delete C_nM as it's in the input...

        rho_mm = self.mmdescriptor.zeros(dtype=dtype)
        
        pblas_simple_gemm(self.mmdescriptor,
                          self.mmdescriptor,
                          self.mmdescriptor,
                          Cf_mm, C_mm, rho_mm, transa='T')
        return rho_mm

    def calculate_density_matrix(self, f_n, C_nM, rho_mM=None):
        """Calculate density matrix from occupations and coefficients.

        Presently this function performs the usual scalapack 3-step trick:
        redistribute-numbercrunching-backdistribute.
        
        
        Notes on future performance improvement.
        
        As per the current framework, C_nM exists as copies on each
        domain, i.e. this is not parallel over domains.  We'd like to
        correct this and have an efficient distribution using e.g. the
        block communicator.

        The diagonalization routine and other parts of the code should
        however be changed to accommodate the following scheme:
        
        Keep coefficients in C_mm form after the diagonalization.
        rho_mm can then be directly calculated from C_mm without
        redistribution, after which we only need to redistribute
        rho_mm across domains.
        
        """
        dtype = C_nM.dtype
        rho_mm = self.calculate_blocked_density_matrix(f_n, C_nM)
        rback = Redistributor(self.block_comm, self.mmdescriptor,
                              self.mM_unique_descriptor)
        rho1_mM = self.mM_unique_descriptor.zeros(dtype=dtype)
        rback.redistribute(rho_mm, rho1_mM)
        del rho_mm

        if rho_mM is None:
            if self.gd.rank == 0:
                rho_mM = rho1_mM
            else:
                rho_mM = self.mMdescriptor.zeros(dtype=dtype)

        self.gd.comm.broadcast(rho_mM, 0)
        return rho_mM

    def distribute_to_columns(self, rho_mm, srcdescriptor):
        redistributor = Redistributor(self.block_comm, # XXX
                                      srcdescriptor,
                                      self.mM_unique_descriptor)
        rho_mM = redistributor.redistribute(rho_mm)
        if self.gd.rank != 0:
            rho_mM = self.mMdescriptor.zeros(dtype=rho_mm.dtype)
        self.gd.comm.broadcast(rho_mM, 0)
        return rho_mM

    def oldcalculate_density_matrix(self, f_n, C_nM, rho_mM=None):
        # This version is parallel over the band descriptor only.
        # This is inefficient, but let's keep it for a while in case
        # there's trouble with the more efficient version
        nbands = self.bd.nbands
        mynbands = self.bd.mynbands
        nao = self.nao
        
        if rho_mM is None:
            rho_mM = self.mMdescriptor.zeros(dtype=C_nM.dtype)
        
        Cf_nM = (C_nM * f_n[:, None]).conj()
        pblas_simple_gemm(self.nMdescriptor, self.nMdescriptor,
                          self.mMdescriptor, Cf_nM, C_nM, rho_mM, transa='T')
        return rho_mM

    def get_transposed_density_matrix(self, f_n, C_nM, rho_mM=None):
        return self.calculate_density_matrix(f_n, C_nM, rho_mM).conj()

    def get_description(self):
        (title, template) = BlacsLayouts.get_description(self)
        bg = self.blockgrid
        desc = self.mmdescriptor
        s = template % (bg.nprow, bg.npcol, desc.mb, desc.nb)
        return ' '.join([title, s])


class OrbitalLayouts(KohnShamLayouts):
    def __init__(self, gd, bd, block_comm, dtype, nao,
                 timer=nulltimer):
        KohnShamLayouts.__init__(self, gd, bd, block_comm, dtype,
                                 timer)
        self.mMdescriptor = MatrixDescriptor(nao, nao)
        self.nMdescriptor = MatrixDescriptor(bd.mynbands, nao)
        
        self.Mstart = 0
        self.Mstop = nao
        self.Mmax = nao
        self.mynao = nao
        self.nao = nao
        self.orbital_comm = bd.comm

    def diagonalize(self, H_MM, C_nM, eps_n, S_MM):
        eps_M = np.empty(C_nM.shape[-1])
        self.block_comm.broadcast(H_MM, 0)
        self.block_comm.broadcast(S_MM, 0)
        # The result on different processor is not necessarily bit-wise
        # identical, so only domain master performs computation
        if self.gd.comm.rank == 0:
            self._diagonalize(H_MM, S_MM.copy(), eps_M)
        self.gd.comm.broadcast(H_MM, 0)
        self.gd.comm.broadcast(eps_M, 0)
        eps_n[:] = eps_M[self.bd.get_slice()]
        C_nM[:] = H_MM[self.bd.get_slice()]
    
    def _diagonalize(self, H_MM, S_MM, eps_M):
        """Serial diagonalize via LAPACK."""
        # This is replicated computation but ultimately avoids
        # additional communication
        general_diagonalize(H_MM, eps_M, S_MM)

    def estimate_memory(self, mem, dtype):
        nao = self.setups.nao
        itemsize = mem.itemsize[dtype]
        mem.subnode('eps [M]', self.nao * mem.floatsize)
        mem.subnode('H [MM]', self.nao * self.nao * itemsize)

    def distribute_overlap_matrix(self, S_qMM, root=0,
                                  add_hermitian_conjugate=False):
        self.gd.comm.sum(S_qMM, root)
        if add_hermitian_conjugate:
            S_qMM += S_qMM.swapaxes(-1, -2).conj()
        return S_qMM

    def get_overlap_matrix_shape(self):
        return self.nao, self.nao

    def calculate_density_matrix(self, f_n, C_nM, rho_MM=None, C2_nM=None):
        # Only a madman would use a non-transposed density matrix.
        # Maybe we should use the get_transposed_density_matrix instead
        if rho_MM is None:
            rho_MM = np.zeros((self.mynao, self.nao), dtype=C_nM.dtype)
        # XXX Should not conjugate, but call gemm(..., 'c')
        # Although that requires knowing C_Mn and not C_nM.
        # that also conforms better to the usual conventions in literature
        if C2_nM is None:
            C2_nM = C_nM
        Cf_Mn = np.ascontiguousarray(C2_nM.T.conj() * f_n)
        gemm(1.0, C_nM, Cf_Mn, 0.0, rho_MM, 'n')
        return rho_MM

    def get_transposed_density_matrix(self, f_n, C_nM, rho_MM=None):
        return self.calculate_density_matrix(f_n, C_nM, rho_MM).T.copy()

        #if rho_MM is None:
        #    rho_MM = np.zeros((self.mynao, self.nao), dtype=C_nM.dtype)
        #C_Mn = C_nM.T.copy()
        #gemm(1.0, C_Mn, f_n[np.newaxis, :] * C_Mn, 0.0, rho_MM, 'c')
        #self.bd.comm.sum(rho_MM)
        #return rho_MM

    def alternative_calculate_density_matrix(self, f_n, C_nM, rho_MM=None):
        if rho_MM is None:
            rho_MM = np.zeros((self.mynao, self.nao), dtype=C_nM.dtype)
        # Alternative suggestion. Might be faster. Someone should test this
        C_Mn = C_nM.T.copy()
        r2k(0.5, C_Mn, f_n * C_Mn, 0.0, rho_MM)
        tri2full(rho_MM)
        return rho_MM

    def get_description(self):
        return 'Serial LAPACK'

    def calculate_density_matrix_delta(self, d_nn, C_nM, rho_MM=None):
        # Only a madman would use a non-transposed density matrix.
        # Maybe we should use the get_transposed_density_matrix instead
        if rho_MM is None:
            rho_MM = np.zeros((self.mynao, self.nao), dtype=C_nM.dtype)
        Cd_Mn = np.zeros((self.nao, self.bd.mynbands), dtype=C_nM.dtype)
        # XXX Should not conjugate, but call gemm(..., 'c')
        # Although that requires knowing C_Mn and not C_nM.
        # that also conforms better to the usual conventions in literature
        C_Mn = C_nM.T.conj().copy()
        gemm(1.0, d_nn, C_Mn, 0.0, Cd_Mn, 'n')
        gemm(1.0, C_nM, Cd_Mn, 0.0, rho_MM, 'n')
        self.bd.comm.sum(rho_MM)
        return rho_MM

    def get_transposed_density_matrix_delta(self, d_nn, C_nM, rho_MM=None):
        return self.calculate_density_matrix_delta(d_nn, C_nM, rho_MM).T.copy()
