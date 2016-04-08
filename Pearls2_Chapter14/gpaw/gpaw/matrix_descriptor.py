
import numpy as np

from gpaw import debug
from gpaw.utilities.lapack import general_diagonalize


class MatrixDescriptor:
    """Class representing a 2D matrix shape.  Base class for parallel
    matrix descriptor with BLACS."""
    
    def __init__(self, M, N):
        self.shape = (M, N)
    
    def __nonzero__(self):
        return self.shape[0] != 0 and self.shape[1] != 0

    def zeros(self, n=(), dtype=float):
        """Return array of zeroes with the correct size on all CPUs.

        The last two dimensions will be equal to the shape of this
        descriptor.  If specified as a tuple, can have any preceding
        dimension."""
        return self._new_array(np.zeros, n, dtype)

    def empty(self, n=(), dtype=float):
        """Return array of zeros with the correct size on all CPUs.

        See zeros()."""
        return self._new_array(np.empty, n, dtype)

    def _new_array(self, func, n, dtype):
        if isinstance(n, int):
            n = n,
        shape = n + self.shape
        return func(shape, dtype)

    def check(self, a_mn):
        """Check that specified array is compatible with this descriptor."""
        return a_mn.shape == self.shape and a_mn.flags.contiguous

    def checkassert(self, a_mn):
        ok = self.check(a_mn)
        if not ok:
            if not a_mn.flags.contiguous:
                msg = 'Matrix with shape %s is not contiguous' % (a_mn.shape,)
            else:
                msg = ('%s-descriptor incompatible with %s-matrix' %
                       (self.shape, a_mn.shape))
            raise AssertionError(msg)

    def general_diagonalize_dc(self, H_mm, S_mm, C_mm, eps_M,
                               UL='L', iu=None):
        general_diagonalize(H_mm, eps_M, S_mm, iu=iu)
        C_mm[:] = H_mm

    def my_blocks(self, array_mn):
        yield (0, self.shape[0], 0, self.shape[1], array_mn)

    def estimate_memory(self, mem, dtype):
        """Handled by subclass."""
        pass
# -------------------------------------------------------------------

class BandMatrixDescriptor(MatrixDescriptor):
    """Descriptor-class for square matrices of bands times bands."""

    def __init__(self, bd, gd, ksl):
        MatrixDescriptor.__init__(self, bd.nbands, bd.nbands)
        self.bd = bd
        self.gd = gd #XXX used?
        self.ksl = ksl # not really used...

    def assemble_blocks(self, A_qnn, A_NN, hermitian):
        """Assign all distributed sub-blocks pertaining from various rank to
        the relevant parts of a Hermitian or non-Hermitian matrix A_NN.

        Parameters:

        A_qnn: ndarray
            Sub-blocks belonging to the specified rank.
        A_NN: ndarray
            Full matrix in which to write contributions from sub-blocks.
        hermitian: bool
            Indicates whether A_NN is a Hermitian matrix, in which
            case only the lower triangular part is assigned to.

        Note that the sub-block buffers are used for communicating across the
        band communicator, hence A_qnn will be altered during the assembly.
        """
        if self.bd.comm.size == 1:
            if hermitian:
                self.triangular_blockwise_assign(A_qnn, A_NN, 0)
            else:
                self.full_blockwise_assign(A_qnn, A_NN, 0)
            return

        if self.bd.comm.rank == 0:
            for band_rank in range(self.bd.comm.size):
                if band_rank > 0:
                    self.bd.comm.receive(A_qnn, band_rank, 13)
                if hermitian:
                    self.triangular_blockwise_assign(A_qnn, A_NN, band_rank)
                else:
                    self.full_blockwise_assign(A_qnn, A_NN, band_rank)
        else:
            self.bd.comm.send(A_qnn, 0, 13)

    def triangular_blockwise_assign(self, A_qnn, A_NN, band_rank):
        """Assign the sub-blocks pertaining from a given rank to the lower
        triangular part of a Hermitian matrix A_NN. This subroutine is used
        for matrix assembly.

        Parameters:

        A_qnn: ndarray
            Sub-blocks belonging to the specified rank.
        A_NN: ndarray
            Full matrix in which to write contributions from sub-blocks.
        band_rank: int
            Communicator rank to which the sub-blocks belongs.

        Note that a Hermitian matrix requires Q=B//2+1 blocks of M x M
        elements where B is the communicator size and M=N//B for N bands.
        """
        N = self.bd.mynbands
        B = self.bd.comm.size
        assert band_rank in xrange(B)

        if B == 1:
            # Only fill in the lower part
            mask = np.tri(N).astype(bool)
            A_NN[mask] = A_qnn.reshape((N,N))[mask]
            return

        # A_qnn[q2,myn1,myn2] on rank q1 is the q2'th overlap calculated
        # between <psi_n1| and A|psit_n2> where n1 <-> (q1,myn1) and 
        # n2 <-> ((q1+q2)%B,myn2) since we've sent/recieved q2 times.
        q1 = band_rank
        Q = B // 2 + 1
        if debug:
            assert A_qnn.shape == (Q,N,N)

        # Note that for integer inequalities, these relations are useful (X>0):
        #     A*X > B   <=>   A > B//X   ^   A*X <= B   <=>   A <= B//X

        if self.bd.strided:
            A_nbnb = A_NN.reshape((N, B, N, B))
            mask = np.empty((N,N), dtype=bool)
            for q2 in range(Q):
                # n1 = (q1+q2)%B + myn1*B   ^   n2 = q1 + myn2*B
                #
                # We seek the lower triangular part i.e. n1 >= n2
                #   <=>   (myn2-myn1)*B <= (q1+q2)%B-q1
                #   <=>   myn2-myn1 <= dq//B
                dq = (q1+q2)%B-q1 # within ]-B; Q[ so dq//B is -1 or 0

                # Create mask for lower part of current block
                mask[:] = np.tri(N, N, dq//B)
                if debug:
                    m1,m2 = np.indices((N,N))
                    assert dq in xrange(-B+1,Q)
                    assert (mask == (m1 >= m2 - dq//B)).all()

                # Copy lower part of A_qnn[q2] to its rightfull place
                A_nbnb[:, (q1+q2)%B, :, q1][mask] = A_qnn[q2][mask]

                # Negate the transposed mask to get complementary mask
                mask = ~mask.T

                # Copy upper part of Hermitian conjugate of A_qnn[q2]
                A_nbnb[:, q1, :, (q1+q2)%B][mask] = A_qnn[q2].T.conj()[mask]
        else:
            A_bnbn = A_NN.reshape((B, N, B, N))

            # Optimization for the first block
            if q1 == 0:
                A_bnbn[:Q, :, 0] = A_qnn
                return

            for q2 in range(Q):
                # n1 = ((q1+q2)%B)*N + myn1   ^   n2 = q1*N + myn2
                #
                # We seek the lower triangular part i.e. n1 >= n2
                #   <=>   ((q1+q2)%B-q1)*N >= myn2-myn1
                #   <=>   myn2-myn1 <= dq*N
                #   <=>   entire block if dq > 0,
                #   ...   myn2 <= myn1 if dq == 0,
                #   ...   copy nothing if dq < 0
                if q1 + q2 < B:
                    A_bnbn[q1 + q2, :, q1] = A_qnn[q2]
                else:
                    A_bnbn[q1, :, q1 + q2 - B] = A_qnn[q2].T.conj()

    def full_blockwise_assign(self, A_qnn, A_NN, band_rank):
        """Assign the sub-blocks pertaining from a given rank to the full
        non-Hermitian matrix A_NN. This subroutine is used for matrix assembly.

        Parameters:

        A_qnn: ndarray
            Sub-blocks belonging to the specified rank.
        A_NN: ndarray
            Full matrix, in which to write contributions from sub-blocks.
        band_rank: int
            Communicator rank to which the sub-blocks belongs.

        Note that a non-Hermitian matrix requires Q=B blocks of M x M
        elements where B is the communicator size and M=N//B for N bands.
        """
        N = self.bd.mynbands
        B = self.bd.comm.size
        assert band_rank in xrange(B)

        if B == 1:
            A_NN[:] = A_qnn.reshape((N,N))
            return

        # A_qnn[q2,myn1,myn2] on rank q1 is the q2'th overlap calculated
        # between <psi_n1| and A|psit_n2> where n1 <-> (q1,myn1) and 
        # n2 <-> ((q1+q2)%B,myn2) since we've sent/recieved q2 times.
        q1 = band_rank
        Q = B
        if debug:
            assert A_qnn.shape == (Q,N,N)

        if self.bd.strided:
            A_nbnb = A_NN.reshape((N, B, N, B))
            for q2 in range(Q):
                A_nbnb[:, (q1+q2)%B, :, q1] = A_qnn[q2]
        else:
            A_bnbn = A_NN.reshape((B, N, B, N))

            # Optimization for the first block
            if q1 == 0:
                A_bnbn[:Q, :, 0] = A_qnn
                return

            for q2 in range(Q):
                A_bnbn[(q1+q2)%B, :, q1] = A_qnn[q2]

    def extract_block(self, A_NN, q1, q2):
        """Extract the sub-block pertaining from a given pair of ranks within
        the full matrix A_NN. Extraction may result in copies to assure unit
        stride, thus one should not utilize this routine for altering A_NN.

        Parameters:

        A_NN: ndarray
            Full matrix, from which to read the requested sub-block.
        q1: int
            Communicator rank to which the sub-block belongs (row index).
        q2: int
            Communicator rank the sub-block originated from (column index).

        Note that a Hermitian matrix requires just Q=B//2+1 blocks of M x M
        elements where B is the communicator size and M=N//B for N bands.
        Therefor, care should be taken to only request q1,q2 pairs which
        are connected by Q shifts or less if A_NN is lower triangular.
        """
        N = self.bd.mynbands
        B = self.bd.comm.size

        if B == 1:
            return A_NN

        if self.bd.strided:
            A_nbnb = A_NN.reshape((N, B, N, B))
            return A_nbnb[:, q2, :, q1].copy() # last dim must have unit stride
        else:
            A_bnbn = A_NN.reshape((B, N, B, N))
            return A_bnbn[q2, :, q1]

    def redistribute_input(self, A_NN): # do nothing
        if debug:
            self.checkassert(A_NN)
        return A_NN

    def redistribute_output(self, A_NN): # do nothing
        if debug:
            self.checkassert(A_NN)
        return A_NN

    def estimate_memory(self, mem, dtype):
        # Temporary work arrays included in estimate #
        nbands = self.bd.nbands
        itemsize = mem.itemsize[dtype]
        mem.subnode('A_NN', nbands*nbands*itemsize)

# -------------------------------------------------------------------

#from gpaw.blacs import BlacsDescriptor #TODO XXX derive from BlacsDescriptor
#from gpaw.blacs import BlacsBandLayouts

class BlacsBandMatrixDescriptor(MatrixDescriptor):#, BlacsBandLayouts):
    """Descriptor-class for square BLACS matrices of bands times bands."""

    def __init__(self, bd, gd, ksl): #mcpus, ncpus, blocksize):
        MatrixDescriptor.__init__(self, bd.nbands, bd.mynbands) #XXX a hack...
        #BlacsBandLayouts.__init__(self, gd, bd, mcpus, ncpus, blocksize)
        #BlacsDescriptor.__init__(self, blacsgrid, M, N, mb, nb, rsrc, csrc)
        self.bd = bd
        self.gd = gd #XXX used?
        self.ksl = ksl

    def assemble_blocks(self, A_qnn, A_Nn, hermitian):
        """Assign all distributed sub-blocks pertaining from various rank to
        the relevant parts of a Hermitian or non-Hermitian matrix A_NN.

        Parameters:

        A_qnn: ndarray
            Sub-blocks belonging to the specified rank.
        A_Nn: ndarray
            Full column vector in which to write contributions from sub-blocks.
        hermitian: bool
            Indicates whether A_Nn represents a Hermitian matrix, in which
            case only the lower triangular part is assigned to.

        Note that the sub-block buffers are used for communicating across the
        band communicator, hence A_qnn will be altered during the assembly.
        """

        band_rank = self.bd.comm.rank
        if hermitian:
            self.triangular_columnwise_assign(A_qnn, A_Nn, band_rank)
        else:
            self.full_columnwise_assign(A_qnn, A_Nn, band_rank)

    def triangular_columnwise_assign(self, A_qnn, A_Nn, band_rank):
        """Assign the sub-blocks pertaining from a given rank to the lower
        triangular part of a Hermitian matrix A_NN. This subroutine is used
        for matrix assembly.

        Parameters:

        A_qnn: ndarray
            Sub-blocks belonging to the specified rank.
        A_Nn: ndarray
            Full column vector in which to write contributions from sub-blocks.
        band_rank: int
            Communicator rank to which the sub-blocks belongs.

        Note that a Hermitian matrix requires Q=B//2+1 blocks of M x M
        elements where B is the communicator size and M=N//B for N bands.
        """
        N = self.bd.mynbands
        B = self.bd.comm.size
        assert band_rank in xrange(B)

        if B == 1:
            # Only fill in the lower part
            mask = np.tri(N).astype(bool)
            A_Nn[mask] = A_qnn.reshape((N,N))[mask]
            return

        # A_qnn[q2,myn1,myn2] on rank q1 is the q2'th overlap calculated
        # between <psi_n1| and A|psit_n2> where n1 <-> (q1,myn1) and 
        # n2 <-> ((q1+q2)%B,myn2) since we've sent/recieved q2 times.
        q1 = band_rank
        Q = B // 2 + 1
        if debug:
            assert A_qnn.shape == (Q,N,N)

        # Note that for integer inequalities, these relations are useful (X>0):
        #     A*X > B   <=>   A > B//X   ^   A*X <= B   <=>   A <= B//X

        if self.bd.strided:
            raise NotImplementedError

            """
            A_nbn = A_NN.reshape((N, B, N))
            mask = np.empty((N,N), dtype=bool)
            for q2 in range(Q):
                # n1 = (q1+q2)%B + myn1*B   ^   n2 = q1 + myn2*B
                #
                # We seek the lower triangular part i.e. n1 >= n2
                #   <=>   (myn2-myn1)*B <= (q1+q2)%B-q1
                #   <=>   myn2-myn1 <= dq//B
                dq = (q1+q2)%B-q1 # within ]-B; Q[ so dq//B is -1 or 0

                # Create mask for lower part of current block
                mask[:] = np.tri(N, N, dq//B)
                if debug:
                    m1,m2 = np.indices((N,N))
                    assert dq in xrange(-B+1,Q)
                    assert (mask == (m1 >= m2 - dq//B)).all()

                # Copy lower part of A_qnn[q2] to its rightfull place
                A_nbn[:, (q1+q2)%B][mask] = A_qnn[q2][mask]

                # Negate the transposed mask to get complementary mask
                mask = ~mask.T

                # Copy upper part of Hermitian conjugate of A_qnn[q2]
                A_nbn[:, q1][mask] = A_qnn[q2].T.conj()[mask] #XXX on rank (q1+q2)%B
            """
        else:
            A_bnn = A_Nn.reshape((B, N, N))

            for q2 in range(Q):
                # n1 = ((q1+q2)%B)*N + myn1   ^   n2 = q1*N + myn2
                #
                # We seek the lower triangular part i.e. n1 >= n2
                #   <=>   ((q1+q2)%B-q1)*N >= myn2-myn1
                #   <=>   myn2-myn1 <= dq*N
                #   <=>   entire block if dq > 0,
                #   ...   myn2 <= myn1 if dq == 0,
                #   ...   copy nothing if dq < 0

                if q1 + q2 < B:
                    A_bnn[q1 + q2] = A_qnn[q2]

            reqs = []                
            if q1 < Q-1: # receive from ranks >= Q
                if debug:
                    Q2 = np.arange(Q,B-q1)
                    print('q1=%d, q2: %12s | recv from q1+q2:%12s -> A_bnn%s' % (q1,Q2.tolist(),(q1+Q2).tolist(),(q1+Q2).tolist()))

                for q2 in range(Q, B-q1):
                    rrank = q1 + q2
                    A_nn = A_bnn[q1 + q2]
                    reqs.append(self.bd.comm.receive(A_nn, rrank, block=False))
            elif q1 >= Q: # send to ranks < Q-1
                if debug:
                    Q2 = np.arange(B-q1,B-Q+1)[::-1]
                    print('q1=%d, q2: %12s | send to q1+q2-B:%12s <- A_qnn%s.T.conj()' % (q1,Q2.tolist(),(q1+Q2-B).tolist(),Q2.tolist()))

                for q2 in reversed(range(B-q1, B-Q+1)): # symmetrize comm.
                    srank = q1 + q2 - B
                    sbuf_nn = np.ascontiguousarray(np.conjugate(A_qnn[q2].T)) # always a copy!
                    reqs.append(self.bd.comm.send(sbuf_nn, srank, block=False))
            else:
                if debug:
                    print('q1=%d, do nothing...' % q1)

            self.bd.comm.waitall(reqs)

    def full_columnwise_assign(self, A_qnn, A_Nn, band_rank):
        """Assign the sub-blocks pertaining from a given rank to the columns of
        non-Hermitian matrix A_Nn. This subroutine is used for column assembly.

        Parameters:
        
        A_qnn: ndarray
            Sub-blocks belonging to the specified rank.
        A_Nn: ndarray
            Full column vector, in which to write contributions from sub-blocks.
        band_rank: int
            Communicator rank to which the sub-blocks belongs.

        Note that a non-Hermitian matrix requires Q=B blocks of M x M
        elements where B is the communicator size and M=N//B for N bands.
        """
        N = self.bd.mynbands
        B = self.bd.comm.size
        assert band_rank in xrange(B)

        if B == 1:
            A_Nn[:] = A_qnn.reshape((N,N))
            return

        # A_qnn[q2,myn1,myn2] on rank q1 is the q2'th overlap calculated
        # between <psi_n1| and A|psit_n2> where n1 <-> (q1,myn1) and
        # n2 <-> ((q1+q2)%B,myn2) since we've sent/recieved q2 times.
        q1 = band_rank
        Q = B

        if self.bd.strided:
            A_nbn = A_Nn.reshape((N, B, N))
            for q2 in range(Q):
                A_nbn[:, (q1+q2)%B] = A_qnn[q2]
        else:
            A_bnn = A_Nn.reshape((B, N, N))

            # Optimization for the first block
            if q1 == 0:
                A_bnn[:Q] = A_qnn
                return

            for q2 in range(Q):
                A_bnn[(q1+q2)%B] = A_qnn[q2]

    def extract_block_from_column(self, A_Nn, q1, q2):
        """Extract the sub-block pertaining from a given pair of ranks within
        the full matrix A_NN. Extraction may result in copies to assure unit
        stride, thus one should not utilize this routine for altering A_NN.

        Parameters:

        A_Nn: ndarray
            Full column vector, from which to read the requested sub-block.
        q1: int
            Communicator rank to which the sub-block belongs (row index).
        q2: int
            Communicator rank the sub-block originated from (column index).

        Note that a Hermitian matrix requires just Q=B//2+1 blocks of M x M
        elements where B is the communicator size and M=N//B for N bands.
        Therefor, care should be taken to only request q1,q2 pairs which
        are connected by Q shifts or less if A_NN is lower triangular.
        """
        N = self.bd.mynbands
        B = self.bd.comm.size

        if B == 1:
            return A_Nn

        if q1 == self.bd.comm.rank: #XXX must evaluate the same on all ranks!
            if self.bd.strided:
                A_nbn = A_Nn.reshape((N, B, N))
                return A_nbn[:, q2, :].copy() # block must be contiguous
            else:
                A_bnn = A_Nn.reshape((B, N, N))
                return A_bnn[q2]

        # Time for us to put the cards on the table
        Qs = np.empty((self.bd.comm.size,2), dtype=int)
        self.bd.comm.all_gather(np.array([q1,q2]), Qs)
        Q1, Q2 = Qs.T

        # Block exchange. What do we want and who should we be sending to?
        rrank = q1
        srank = np.argwhere(Q1 == self.bd.comm.rank).ravel().item()
        sq2 = Q2[srank]

        if debug:
            S = np.empty(self.bd.comm.size, dtype=int)
            self.bd.comm.all_gather(np.array([srank]), S)
            if self.bd.comm.rank == 0:
                #print 'Q1: %s, Q2: %s' % (Q1.tolist(),Q2.tolist())
                print('recv(Q1): %s, send(Q1==rank): %s' % (Q1.tolist(),S.tolist()))

        if self.bd.strided:
            A_nbn = A_Nn.reshape((N, B, N))
            sbuf_nn = A_nbn[:, sq2, :].copy() # block must be contiguous
        else:
            A_bnn = A_Nn.reshape((B, N, N))
            sbuf_nn = A_bnn[sq2]

        A_nn = np.empty_like(sbuf_nn)
        self.bd.comm.sendreceive(sbuf_nn, srank, A_nn, rrank)
        return A_nn

    def extract_block_from_row(self, A_nN, q1, q2):
        """Extract the sub-block pertaining from a given pair of ranks within
        the full matrix A_NN. Extraction may result in copies to assure unit
        stride, thus one should not utilize this routine for altering A_NN.

        Parameters:

        A_nN: ndarray
            Full row vector, from which to read the requested sub-block.
        q1: int
            Communicator rank to which the sub-block belongs (row index).
        q2: int
            Communicator rank the sub-block originated from (column index).

        Note that a Hermitian matrix requires just Q=B//2+1 blocks of M x M
        elements where B is the communicator size and M=N//B for N bands.
        Therefor, care should be taken to only request q1,q2 pairs which
        are connected by Q shifts or less if A_NN is lower triangular.
        """
        N = self.bd.mynbands
        B = self.bd.comm.size

        if B == 1:
            return A_nN

        if q2 == self.bd.comm.rank: #XXX must evaluate the same on all ranks!
            if self.bd.strided:
                A_nnb = A_nN.reshape((N, N, B))
                return A_nbn[..., q1].copy() # block must be contiguous
            else:
                A_nbn = A_nN.reshape((N, B, N))
                return A_nbn[:, q1, :].copy() # block must be contiguous
        else:
            raise NotImplementedError

    extract_block = extract_block_from_row #XXX ugly but works

    def redistribute_input(self, A_nn, A_nN=None): # 2D -> 1D row layout
        if A_nN is None:
            A_nN = self.ksl.nNdescriptor.empty(dtype=A_nn.dtype)
        self.ksl.nn2nN.redistribute(A_nn, A_nN)
        if not self.ksl.nNdescriptor.blacsgrid.is_active():
            assert A_nN.shape == (0,0)
            A_nN = np.empty((self.bd.mynbands, self.bd.nbands), dtype=A_nN.dtype)
        self.gd.comm.broadcast(A_nN, 0)
        return A_nN

    def redistribute_output(self, A_Nn, A_nn=None): # 1D column -> 2D layout
        if not self.ksl.Nndescriptor.blacsgrid.is_active():
            A_Nn = np.empty((0,0), dtype=A_Nn.dtype)
        if A_nn is None:
            A_nn = self.ksl.nndescriptor.empty(dtype=A_Nn.dtype)
        self.ksl.Nn2nn.redistribute(A_Nn, A_nn)
        return A_nn
    
    def estimate_memory(self, mem, dtype):
        # Temporary work arrays included in estimate #
        mynbands = self.bd.mynbands
        nbands = self.bd.nbands
        itemsize = mem.itemsize[dtype]
        mem.subnode('2 A_nN', 2*mynbands*nbands*itemsize)
        mem.subnode('2 A_nn', 2*nbands*nbands/self.ksl.blockgrid.ncpus*itemsize)

    #def redistribute_input(self, A_NN): # 2D -> 1D row layout
    #    # XXX instead of a BLACS-distribute from 2D, we disassemble the full matrix
    #    A_nN = np.empty((self.bd.mynbands,self.bd.nbands), dtype=A_NN.dtype)
    #    self.bd.distribute(A_NN, A_nN)
    #    return A_nN

    #def redistribute_input(self, A_NN): # 2D -> 1D column layout
    #    # XXX instead of a BLACS-distribute from 2D, we disassemble the full matrix
    #    A_nN = np.empty((self.bd.mynbands,self.bd.nbands), dtype=A_NN.dtype)
    #    self.bd.distribute(A_NN.T.copy(), A_nN)
    #    return A_nN.T.copy()

    #def redistribute_output(self, A_Nn): # 1D column -> 2D layout
    #    # XXX instead of a BLACS-distribute to 2D, we assemble the full matrix
    #    A_NN = self.bd.collect(A_Nn.T.copy())
    #    if self.bd.comm.rank == 0:
    #        return A_NN.T.copy()
    #    else:
    #        return np.empty((self.bd.nbands,self.bd.nbands), dtype=A_Nn.dtype)

