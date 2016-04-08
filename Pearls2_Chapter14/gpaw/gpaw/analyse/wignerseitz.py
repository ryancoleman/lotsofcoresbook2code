from math import sqrt, pi
import numpy as np

from ase.units import Bohr

from gpaw.utilities import pack, pack2
from gpaw.analyse.hirshfeld import HirshfeldDensity
from gpaw.utilities.tools import coordinates
from gpaw.mpi import MASTER


def wignerseitz(gd, atoms, scale=None):
    """Determine which atom is closest to each grid point.

    The atom distances might be scaled by the scale list."""
    if scale is None:
        scale = [1.] * len(atoms)
    else:
        assert(len(scale) == len(atoms))
    r_vG, R2min_G = coordinates(gd, atoms[0].position / Bohr)
    R2min_G *= scale[0] ** 2
    index_G = gd.zeros(dtype=int)
    for a, atom in enumerate(atoms[1:]):
        r_vG, r2_G = coordinates(gd, atom.position / Bohr)
        r2_G *= scale[a + 1] ** 2
        index_G = np.where(R2min_G > r2_G, a + 1, index_G)
        R2min_G = np.where(R2min_G > r2_G, r2_G, R2min_G)
    return index_G


class WignerSeitz:

    def __init__(self, gd, atoms,
                 calculator=None, scale=None):
        """Find the grid points nearest to the atoms"""

        self.atoms = atoms
        self.gd = gd
        self.calculator = calculator

        self.atom_index = wignerseitz(gd, atoms, scale)

    def expand(self, density):
        """Expand a smooth density in Wigner-Seitz cells around the atoms"""
        n = len(self.atoms)
        weights = np.empty((n,))
        for a in range(n):
            mask = np.where(self.atom_index == a, density, 0.0)
            # XXX Optimize! No need to integrate in zero-region
            weights[a] = self.gd.integrate(mask)

        return weights

    def expand_density(self, nt_G, s, nspins):
        """Get the weights of spin-density in Wigner-Seitz cells
        around the atoms. The spin index and number of spins are
        needed for the augmentation sphere corrections."""
        weights_a = self.expand(nt_G)
        for a, nucleus in enumerate(self.atoms):
            weights_a[a] += nucleus.get_density_correction(s, nspins)
        return weights_a

    def expand_wave_function(self, psit_G, u, n):
        """Get the weights of wave function in Wigner-Seitz cells
        around the atoms. The spin-k-point index u and band number n
        are needed for the augmentation sphere corrections."""

        assert psit_G.dtype == float
        # smooth part
        weigths = self.expand(psit_G ** 2)

        # add augmentation sphere corrections
        for a, nucleus in enumerate(self.atoms):
            P_i = nucleus.P_uni[u, n]
            P_p = pack(np.outer(P_i, P_i))
            Delta_p = sqrt(4 * pi) * nucleus.setup.Delta_pL[:, 0]
            weigths[a] += np.dot(Delta_p, P_p)

        return weigths

    def get_charges(self, den_g):
        """Charge on the atom according to the Wigner-Seitz partitioning

        Can be applied to any density den_g.
        """
        assert(den_g.shape == tuple(self.gd.n_c))
        charges = []
        for atom, q in zip(self.atoms, self.expand(den_g)):
            charges.append(atom.number - q)
        return charges

    def get_effective_volume_ratio(self, atom_index):
        """Effective volume to free volume ratio.

        After: Tkatchenko and Scheffler PRL 102 (2009) 073005
        """
        atoms = self.atoms
        finegd = self.gd

        den_g, gd = self.calculator.density.get_all_electron_density(atoms)
        assert(gd == finegd)
        denfree_g, gd = self.hdensity.get_density([atom_index])
        assert(gd == finegd)

        # the atoms r^3 grid
        position = self.atoms[atom_index].position / Bohr
        r_vg, r2_g = coordinates(finegd, origin=position)
        r3_g = r2_g * np.sqrt(r2_g)

        weight_g = np.where(self.atom_index == atom_index, 1.0, 0.0)

        nom = finegd.integrate(r3_g * den_g[0] * weight_g)
        denom = finegd.integrate(r3_g * denfree_g)

        return nom / denom

    def get_effective_volume_ratios(self):
        """Return the list of effective volume to free volume ratios."""
        ratios = []
        self.hdensity = HirshfeldDensity(self.calculator)
        for a, atom in enumerate(self.atoms):
            ratios.append(self.get_effective_volume_ratio(a))
        return np.array(ratios)


class LDOSbyBand:

    """Base class for a band by band LDOS"""

    def by_element(self):
        # get element indicees
        elemi = {}
        for i, nucleus in enumerate(self.paw.atoms):
            symbol = nucleus.setup.symbol
            if elemi.has_key(symbol):
                elemi[symbol].append(i)
            else:
                elemi[symbol] = [i]
        for key in elemi.keys():
            elemi[key] = self.get(elemi[key])
        return elemi


class WignerSeitzLDOS(LDOSbyBand):

    """Class to get the unfolded LDOS defined by Wigner-Seitz cells"""

    def __init__(self, paw):
        self.paw = paw
        self.ws = WignerSeitz(paw.gd, paw.atoms)

        nu = paw.nkpts * paw.nspins
        ldos = np.empty((nu, paw.nbands, len(paw.atoms)))
        for u, kpt in enumerate(paw.kpt_u):
            for n, psit_G in enumerate(kpt.psit_nG):
                ldos[u, n, :] = ws.expand_wave_function(psit_G, u, n)

    def write(self, filename, slavewrite=False):
        if self.world.rank == MASTER or slavewrite:
            paw = self.paw
            f = open(filename, 'w')

            nn = len(paw.atoms)

            for k in range(paw.nkpts):
                for s in range(paw.nspins):
                    u = s * paw.nkpts + k
                    for n in range(paw.nbands):
                        # avery: Added dummy loop body to make compiling work.
                        1
