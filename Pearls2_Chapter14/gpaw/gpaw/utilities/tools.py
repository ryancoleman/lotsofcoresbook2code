import numpy as np

from ase.units import Hartree, Bohr

def L_to_lm(L):
    """Convert L index to (l, m) index."""
    l = int(np.sqrt(L))
    m = L - l**2 - l
    return l, m


def lm_to_L(l, m):
    """Convert (l, m) index to L index."""
    return l**2 + l + m


def split_formula(formula):
    """Count elements in a chemical formula.

    E.g. split_formula('C2H3Mg') -> ['C', 'C', 'H', 'H', 'H', 'Mg']
    """
    res = []
    for c in formula:
        if c.isupper():
            res.append(c)
        elif c.islower():
            res[-1] += c
        else:
            res.extend([res[-1],] * (eval(c) - 1))
    return res


def construct_reciprocal(gd, q_c=None):
    """Construct the reciprocal lattice from ``GridDescriptor`` instance.

    The generated reciprocal lattice has lattice vectors correspoding to the
    real-space lattice defined in the input grid. Note that it is the squared
    length of the reciprocal lattice vectors that are returned.

    The ordering of the reciprocal lattice agrees with the one typically used
    in fft algorithms, i.e. positive k-values followed by negative.

    Note that the G=(0,0,0) entrance is set to one instead of zero. This
    bit should probably be moved somewhere else ...
    
    Parameters
    ----------
    q_c: ndarray
        Offset for the reciprocal lattice vectors (in scaled coordinates of the
        reciprocal lattice vectors, i.e. array with index ``c``). When
        specified, the returned array contains the values of (q+G)^2 where G
        denotes the reciprocal lattice vectors.
    
    """
    
    assert gd.pbc_c.all(), 'Works only with periodic boundary conditions!'

    # Check q_c
    if q_c is not None:
        assert q_c.shape in [(3,), (3,1)]
        q_c = q_c.reshape((3,1))
        
    # Calculate reciprocal lattice vectors
    N_c1 = gd.N_c[:, np.newaxis]
    i_cq = np.indices(gd.N_c).reshape((3, -1))
    i_cq += N_c1 // 2
    i_cq %= N_c1
    i_cq -= N_c1 // 2

    if q_c is not None:
        i_cq = np.array(i_cq, dtype=float)
        i_cq += q_c

    # Convert from scaled to absolute coordinates
    B_vc = 2.0 * np.pi * gd.icell_cv.T
    k_vq = np.dot(B_vc, i_cq)

    k_vq *= k_vq
    k2_Q = k_vq.sum(axis=0).reshape(gd.N_c)

    if q_c is None:
        k2_Q[0, 0, 0] = 1.0
    elif abs(q_c).sum() < 1e-8:
        k2_Q[0, 0, 0] = 1.0
        
    # Determine N^3
    N3 = gd.n_c[0] * gd.n_c[1] * gd.n_c[2]

    return k2_Q, N3

def coordinates(gd, origin=None, tiny=1e-12):
    """Constructs and returns matrices containing cartesian coordinates,
       and the square of the distance from the origin.

       The origin can be given explicitely (in Bohr units, not Anstroms). 
       Otherwise the origin is placed in the center of the box described 
       by the given grid-descriptor 'gd'.
    """
    
    if origin is None:
        origin = 0.5 * gd.cell_cv.sum(0)
    r0_v = np.array(origin)

    r_vG = gd.get_grid_point_distance_vectors(r0_v)
    r2_G = np.sum(r_vG**2, axis=0)
    # Remove singularity at origin and replace with small number
    r2_G = np.where(r2_G < tiny, tiny, r2_G)

    # Return r^2 matrix
    return r_vG, r2_G

def pick(a_ix, i):
    """Take integer index of a, or a linear combination of the elements of a"""
    if isinstance(i, int):
        return a_ix[i]
    shape = a_ix.shape
    a_x = np.dot(i, a_ix[:].reshape(shape[0], -1))
    return a_x.reshape(shape[1:])


def dagger(a, copy=True):
    """Return Hermitian conjugate of input

    If copy is False, the original array might be overwritten. This is faster,
    but use with care.
    """
    if copy:
        return np.conj(a.T)
    else:
        a = a.T
        if a.dtype == complex:
            a.imag *= -1
        return a


def project(a, b):
    """Return the projection of b onto a."""
    return a * (np.dot(a.conj(), b) / np.linalg.norm(a))


def normalize(U):
    """Normalize columns of U."""
    for col in U.T:
        col /= np.linalg.norm(col)


def get_matrix_index(ind1, ind2=None):
    if ind2 is None:
        dim1 = len(ind1)
        return np.resize(ind1, (dim1, dim1))
    else:
        dim1 = len(ind1)
        dim2 = len(ind2)
    return np.resize(ind1, (dim2, dim1)).T, np.resize(ind2, (dim1, dim2))


def gram_schmidt(U):
    """Orthonormalize columns of U according to the Gram-Schmidt procedure."""
    for i, col in enumerate(U.T):
        for col2 in U.T[:i]:
            col -= col2 * np.dot(col2.conj(), col)
        col /= np.linalg.norm(col)


def lowdin(U, S=None):
    """Orthonormalize columns of U according to the Lowdin procedure.

    If the overlap matrix is know, it can be specified in S.
    """
    if S is None:
        S = np.dot(dagger(U), U)
    eig, rot = np.linalg.eigh(S)
    rot = np.dot(rot / np.sqrt(eig), dagger(rot))
    U[:] = np.dot(U, rot)


def lowdin_svd(U):
    """Orthogonalize according to the Lowdin procedure
       using singular value decomposition.

       U is an N x M matrix containing M vectors as its columns.
    """
    Z, D, V = np.linalg.svd(U, full_matrices=0)
    return np.dot(Z, V)


def symmetrize(matrix):
    """Symmetrize input matrix."""
    np.add(dagger(matrix), matrix, matrix)
    np.multiply(.5, matrix, matrix)
    return matrix


def tri2full(H_nn, UL='L', map=np.conj):
    """Fill in values of hermitian or symmetric matrix.

    Fill values in lower or upper triangle of H_nn based on the opposite
    triangle, such that the resulting matrix is symmetric/hermitian.

    UL='U' will copy (conjugated) values from upper triangle into the
    lower triangle.

    UL='L' will copy (conjugated) values from lower triangle into the
    upper triangle.

    The map parameter can be used to specify a different operation than
    conjugation, which should work on 1D arrays.  Example::

      def antihermitian(src, dst):
            np.conj(-src, dst)

      tri2full(H_nn, map=antihermitian)

    """
    N, tmp = H_nn.shape
    assert N == tmp, 'Matrix must be square'
    #assert np.isreal(H_nn.diagonal()).all(), 'Diagonal should be real'
    if UL != 'L':
        H_nn = H_nn.T

    for n in range(N - 1):
        map(H_nn[n + 1:, n], H_nn[n, n + 1:])


def apply_subspace_mask(H_nn, f_n):
    """Uncouple occupied and unoccupied subspaces.

    This method forces the H_nn matrix into a block-diagonal form
    in the occupied and unoccupied states respectively.
    """
    occ = 0
    nbands = len(f_n)
    while occ < nbands and f_n[occ] > 1e-3: occ += 1
    H_nn[occ:, :occ] = H_nn[:occ, occ:] = 0


def cutoff2gridspacing(E):
    """Convert planewave energy cutoff to a real-space gridspacing."""
    return np.pi / np.sqrt(2 * E / Hartree) * Bohr


def gridspacing2cutoff(h):
    """Convert real-space gridspacing to planewave energy cutoff."""
    # In Hartree units, E = k^2 / 2, where k_max is approx. given by pi / h
    # See PRB, Vol 54, 14362 (1996)
    return 0.5 * (np.pi * Bohr / h)**2 * Hartree


def tridiag(a, b, c, r, u):
    """Solve linear system with tridiagonal coefficient matrix.

    a is the lower band, b is the diagonal, c is the upper band, and
    r is the right hand side.
    The solution is returned in u.


    [b1 c1  0  ...            ] [u1]   [r1]
    [a1 b2 c2 0 ...           ] [ :]   [ :]
    [ 0 a2 b3 c3 0 ...        ] [  ] = [  ]
    [                         ] [  ]   [  ]
    [     ... 0 an-2 bn-1 cn-1] [ :]   [ :]
    [          ... 0 an-1 bn  ] [un]   [rn]
    """
    n = len(b)
    tmp = np.zeros(n-1) # necessary temporary array
    if b[0] == 0:
        raise RuntimeError, 'System is effectively order N-1'

    beta = b[0]
    u[0] = r[0] / beta
    for i in range(1, n):
        # Decompose and forward substitution
        tmp[i-1] = c[i-1] / beta
        beta = b[i] - a[i-1] * tmp[i-1]
        if beta == 0:
            raise RuntimeError, 'Method failure'
        u[i] = (r[i] - a[i-1] * u[i-1]) / beta

    for i in range(n-1, 0, -1):
        # Backward substitution
        u[i-1] -= tmp[i-1] * u[i]


def signtrim(data, decimals=None):
    """Trim off the sign of potential zeros, usually occuring after round.

    data is the ndarray holding NumPy data to round and trim.
    decimals is an integer specifying how many decimals to round to.
    """
    if decimals is not None:
        data = data.round(decimals) #np.round is buggy because -0 != 0

    shape = data.shape
    data = data.reshape(-1)

    if data.dtype == complex:
        i = np.argwhere(np.sign(data.real) == 0).ravel()
        j = np.argwhere(np.sign(data.imag) == 0).ravel()
        data.real[i] = 0
        data.imag[j] = 0
    else:
        i = np.argwhere(np.sign(data)==0).ravel()
        data[i] = 0

    return data.reshape(shape)


try:
    from hashlib import md5 as md5_new
    import hashlib as md5
except ImportError:
    from md5 import new as md5_new
    import md5


def md5_array(data, numeric=False):
    """Create MD5 hex digest from NumPy array.

    Optionally, will cast the 128 bit hexadecimal hash to a 64 bit integer.

    Warning: For floating point types, only bitwise identical data will
    produce matching MD5 fingerprints, so do not attempt to match sets
    of nearly identical data by rounding off beforehand.

    Example:

     >>> data = np.linspace(0,np.pi,1000000)
     >>> eps = 1e-6
     >>> a = md5_array(data.round(3))
     >>> b = md5_array((data+eps).round(3))
     >>> assert a==b, 'Mismatch between %s and %s' % (a,b)

    This is due to the inexact nature of the floating point representation.
    """

    if not isinstance(data, np.ndarray):
        data = np.asarray(data)

    # Only accepts float,complex,int,bool,...
    if (not np.issubdtype(data.dtype, np.number) and
        data.dtype not in [bool, np.bool, np.bool_]):
        raise TypeError('MD5 hex digest only accepts numeric/boolean arrays.')

    datahash = md5.md5(data.tostring())

    if numeric:
        xor = lambda a,b: chr(ord(a)^ord(b)) # bitwise xor on 2 bytes -> 1 byte
        sbuf128 = datahash.digest()
        sbuf64 = ''.join([xor(a,b) for a,b in zip(sbuf128[::2],sbuf128[1::2])])
        return np.fromstring(sbuf64, np.int64).item()
    else:
        return datahash.hexdigest()


def split_nodes(length, parrank, parsize):
    """Split length over nodes.

    Divide length into parsize even sized chunks, and return the start/end
    indices of the parrank'th chunk.
    """
    if parsize == 1:
        return 0, length
    pernode = int(round(length / float(parsize)))
    if parrank == parsize - 1:
        return parrank * pernode, length
    return parrank * pernode, (parrank + 1) * pernode


class Spline:
    def __init__(self, xi, yi, leftderiv=None, rightderiv=None):
        """Cubic spline approximation class.

        xi, yi specifies the known data points.

        leftderiv and rightderiv specifies the first derivative on the
        boundaries. If set to None, the second derivative is set to zero.

        Example usage::

          >>> xi = arange(.1, 5, .5)    # known data points
          >>> yi = cos(xi)              # known data points
          >>> sp = Spline(xi, yi)       # make spline
          >>> x = arange(-.5, 5.5, .05) # points to interpolate to
          >>> y = sp(x)  # get spline value on an entire list
          >>> y2 = sp(4) # get spline value at a single point

        Based on 'Numerical recipes in c'
        """
        self.xy = (xi, yi)
        N = len(xi)
        self.ypp = u = np.zeros(N) # The second derivatives y''
        tmp = np.zeros(N - 1)

        # Set left boundary condition
        if leftderiv is None: # natural spline - second derivative is zero
            tmp[0] = u[0] = 0.0
        else: # clamped spline - first derivative is fixed
            tmp[0] = 3 / (xi[1] - xi[0]) * (
                (yi[1] - yi[0]) / (xi[1] - xi[0]) - leftderiv)
            u[0] = -.5

        for i in range(1, N - 1):
            sig = (xi[i] - xi[i - 1]) / (xi[i + 1] - xi[i - 1])
            p = sig * u[i - 1] + 2
            u[i] = (sig - 1) / p
            tmp[i] = (yi[i + 1] - yi[i]) / (xi[i + 1] - xi[i]) - \
                     (yi[i] - yi[i - 1]) / (xi[i] - xi[i - 1])
            tmp[i] = (6 * tmp[i] / (xi[i +1] - xi[i-1]) - sig * tmp[i - 1]) / p

        # Set right boundary condition
        if rightderiv is None: # natural spline - second derivative is zero
            qn = tmpn = 0.0
        else: # clamped spline - first derivative is fixed
            qn = .5
            tmpn = 3 / (xi[N - 1] - xi[N - 2]) * (
                rightderiv - (yi[N - 1] - yi[N - 2]) / (xi[N - 1] - xi[N - 2]))

        u[N - 1] = (tmpn - qn * tmp[N - 2]) / (qn * u[N - 1] + 1)
        for k in range(N - 2, -1, -1): # backsubstitution step
            u[k] = u[k] * u[k + 1] + tmp[k]

    def __call__(self, x):
        """Evaluate spline for each point in input argument.

        The value in point x[i-1] < x <= x[i] is determined by::

                                    ''       ''
          y(x) = a y    + b y  + c y    + d y
                    i-1      i      i-1      i

        """
        x = np.array(x, float)
        if x.ndim == 0: x.shape = (1,)
        y = np.zeros_like(x)
        xi, yi = self.xy

        i = None
        for j, xval in enumerate(x):
            i = self.locate(xval, i)
            h = xi[i] - xi[i - 1]
            a = (xi[i] - xval) / h
            b = 1. - a
            c = (a**3 - a) * h**2 / 6.
            d = (b**3 - b) * h**2 / 6.
            y[j] = (a * yi[i - 1] + b * yi[i] +
                    c * self.ypp[i - 1] + d * self.ypp[i])
        return y

    def locate(self, x, guess=None):
        """return i such that x[i-1] < x <= xi[i]

        1 or len(xi) - 1 is returned if x is outside list range.
        """
        xi = self.xy[0]
        if x <= xi[0]: return 1
        elif x > xi[-1]: return len(xi) - 1
        elif guess and xi[guess - 1] < x <= xi[guess]: return guess
        else: return np.searchsorted(xi, x)
