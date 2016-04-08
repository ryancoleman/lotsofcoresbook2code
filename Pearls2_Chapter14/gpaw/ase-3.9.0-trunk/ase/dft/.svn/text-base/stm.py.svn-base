import pickle

import numpy as np


class STM:
    def __init__(self, atoms, symmetries=None):
        """Scanning tunneling microscope.

        atoms: Atoms object or filename
            Atoms to scan or name of file to read LDOS from.
        symmetries: list of int
            List of integers 0, 1, and/or 2 indicating which surface
            symmetries have been used to reduce the number of k-points
            for the DFT calculation.  The three integers correspond to
            the following three symmetry operations::

                 [-1  0]   [ 1  0]   [ 0  1]
                 [ 0  1]   [ 0 -1]   [ 1  0]
        """

        if isinstance(atoms, str):
            with open(atoms, 'rb') as f:
                self.ldos, self.bias, self.cell = pickle.load(f)
            self.atoms = None
        else:
            self.atoms = atoms
            self.cell = atoms.cell
            self.bias = None
            self.ldos = None
            assert not self.cell[2, :2].any() and not self.cell[:2, 2].any()

        self.symmetries = symmetries or []
                               
    def calculate_ldos(self, bias):
        """Calculate local density of states for given bias."""
        if self.ldos is not None and bias == self.bias:
            return

        if bias < 0:
            emin = bias
            emax = 0.0
        else:
            emin = 0
            emax = bias

        calc = self.atoms.calc

        nbands = calc.get_number_of_bands()
        weights = calc.get_k_point_weights()
        nkpts = len(weights)
        nspins = calc.get_number_of_spins()
        eigs = np.array([[calc.get_eigenvalues(k, s)
                          for k in range(nkpts)]
                         for s in range(nspins)])
        eigs -= calc.get_fermi_level()
        ldos = 0.0
        for s in range(nspins):
            for k in range(nkpts):
                for n in range(nbands):
                    e = eigs[s, k, n]
                    if emin < e < emax:
                        psi = calc.get_pseudo_wave_function(n, k, s)
                        ldos += weights[k] * (psi * np.conj(psi)).real

        if 0 in self.symmetries:
            # (x,y) -> (-x,y)
            ldos[1:] += ldos[:0:-1].copy()
            ldos[1:] *= 0.5

        if 1 in self.symmetries:
            # (x,y) -> (x,-y)
            ldos[:, 1:] += ldos[:, :0:-1].copy()
            ldos[:, 1:] *= 0.5
            
        if 2 in self.symmetries:
            # (x,y) -> (y,x)
            ldos += ldos.transpose((1, 0, 2)).copy()
            ldos *= 0.5
            
        self.ldos = ldos
        self.bias = bias

    def write(self, filename='stm.pckl'):
        """Write local density of states to pickle file."""
        with open(filename, 'wb') as f:
            pickle.dump((self.ldos, self.bias, self.cell), f,
                        protocol=pickle.HIGHEST_PROTOCOL)
        
    def get_averaged_current(self, bias, z):
        """Calculate avarage current at height z.

        Use this to get an idea of what current to use when scanning."""

        self.calculate_ldos(bias)
        nz = self.ldos.shape[2]

        # Find grid point:
        n = z / self.cell[2, 2] * nz
        dn = n - np.floor(n)
        n = int(n) % nz

        # Average and do linear interpolation:
        return ((1 - dn) * self.ldos[:, :, n].mean() +
                dn * self.ldos[:, :, (n + 1) % nz].mean())
    
    def scan(self, bias, current, z0=None, repeat=(1, 1)):
        """Constant current 2-d scan.
        
        Returns three 2-d arrays (x, y, z) containing x-coordinates,
        y-coordinates and heights.  These three arrays can be passed to
        matplotlibs contourf() function like this:
            
        >>> import matplotlib.pyplot as plt
        >>> plt.gca(aspect='equal')
        >>> plt.contourf(x, y, z)
        >>> plt.show()
        
        """
        
        self.calculate_ldos(bias)

        L = self.cell[2, 2]
        nz = self.ldos.shape[2]
        h = L / nz

        ldos = self.ldos.reshape((-1, nz))

        heights = np.empty(ldos.shape[0])
        for i, a in enumerate(ldos):
            heights[i] = find_height(a, current, h, z0)

        s0 = heights.shape = self.ldos.shape[:2]
        heights = np.tile(heights, repeat)
        s = heights.shape

        ij = np.indices(s, dtype=float).reshape((2, -1)).T
        x, y = np.dot(ij / s0, self.cell[:2, :2]).T.reshape((2,) + s)

        return x, y, heights
    
    def linescan(self, bias, current, p1, p2, npoints=50, z0=None):
        """Constant current line scan.

        Example::

            stm = STM(...)
            z = ...  # tip position
            c = stm.get_averaged_current(-1.0, z)
            stm.linescan(-1.0, c, (1.2, 0.0), (1.2, 3.0))
        """

        heights = self.scan(bias, current, z0)[2]

        p1 = np.asarray(p1, float)
        p2 = np.asarray(p2, float)
        d = p2 - p1
        s = np.dot(d, d)**0.5

        cell = self.cell[:2, :2]
        shape = np.array(heights.shape, float)
        M = np.linalg.inv(cell)
        line = np.empty(npoints)
        for i in range(npoints):
            p = p1 + i * d / (npoints - 1)
            q = np.dot(p, M) * shape
            line[i] = interpolate(q, heights)
        return np.linspace(0, s, npoints), line
        
        
def interpolate(q, heights):
    qi = q.astype(int)
    f = q - qi
    g = 1 - f
    qi %= heights.shape
    n0, m0 = qi
    n1, m1 = (qi + 1) % heights.shape
    z = (g[0] * g[1] * heights[n0, m0] +
         f[0] * g[1] * heights[n1, m0] +
         g[0] * f[1] * heights[n0, m1] +
         f[0] * f[1] * heights[n1, m1])
    return z
     
     
def find_height(ldos, current, h, z0=None):
    if z0 is None:
        n = len(ldos) - 1
    else:
        n = int(z0 / h)
    while n >= 0:
        if ldos[n] > current:
            break
        n -= 1
    else:
        return 0.0

    c2, c1 = ldos[n:n + 2]
    return (n + 1 - (current - c1) / (c2 - c1)) * h
