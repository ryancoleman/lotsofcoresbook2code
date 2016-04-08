"""Part of the module for electrodynamic simulations

"""

from ase.units import Hartree, Bohr, _eps0, _c, _aut
from gpaw import PoissonConvergenceError
from gpaw.fd_operators import Gradient
from gpaw.io import open as gpaw_io_open
from gpaw.tddft.units import attosec_to_autime, autime_to_attosec, \
    autime_to_attosec
from gpaw.transformers import Transformer
from gpaw.utilities import mlsqr
from gpaw.utilities.blas import axpy
from math import pi
from numpy.fft import fftn, ifftn, fft2, ifft2
from string import split
import _gpaw
import gpaw.mpi as mpi
import numpy as np
import sys

# in atomic units, 1/(4*pi*e_0) = 1
_eps0_au = 1.0 / (4.0 * np.pi)

# Base class for the classical polarizable material:
#    -holds various functions at each point in the calculation grid:
#        1) the dielectric function (permittivity)
#        2) electric field
#        3) classical polarization charge density
#    -contains routines for calculating them from each other and/or external potential
class PolarizableMaterial():
    def __init__(self, components=None, sign = -1.0):
        self.gd          = None
        self.initialized = False
        self.sign        = sign
        self.messages    = []

        if components is None:
            self.components = []
        else:
            self.components  = components
    
    
    def add_component(self, component):
        self.components.append(component)

    def permittivity_value(self, omega=0.0):
        return self.eps_infty + _eps0_au * np.sum(self.beta / (self.bar_omega**2.0 - 1J * self.alpha * omega - omega**2.0), axis=0)

    def get_static_permittivity(self):
        return self.gd.collect(self.permittivity_value(0.0))

    def initialize(self, gd):
        if self.initialized: # double initialization leads to problems
            return
        self.initialized = True
        self.messages.append("Polarizable Material:")
        
        try:
            self.Nj = max(component.permittivity.Nj for component in self.components)
        except:
            self.Nj = 0
            
        self.gd = gd
        
        # 3-dimensional scalar array: rho, eps_infty
        self.charge_density = self.gd.zeros()
        self.eps_infty = np.ones(self.gd.empty().shape) * _eps0_au
        
        # 3-dimensional vector arrays:
        #        electric field, total polarization density
        dims = [3] + list(self.gd.empty().shape)
        self.electric_field = np.zeros(dims)
        self.polarization_total = np.zeros(dims)
        
        # 4-dimensional vector arrays:
        #        currents, polarizations
        dims = [3, self.Nj] + list(self.gd.empty().shape)
        self.currents      = np.zeros(dims)
        self.polarizations = np.zeros(dims)
        
        # 4-dimensional scalar arrays:
        #        oscillator parameters alpha, beta, bar_omega, eps_infty
        dims = [self.Nj] + list(self.gd.empty().shape)
        self.alpha     = np.zeros(dims)
        self.beta      = np.zeros(dims)
        self.bar_omega = np.ones(dims)
        
        # Set the permittivity for each grid point
        for component in self.components:
            self.apply_mask(mask = component.get_mask(self.gd),
                           permittivity = component.permittivity)
    
    # Restart by regenerating the structures
    # TODO: is everything always fully reproduced? Should there be a test for it (e.g. checksum[mask])
    def read(self, reader):
        r = reader
        num_components = r['classmat.num_components']
        for n in range(num_components):
            name = r['classmat.component_%i.name' % n]
            arguments = r['classmat.component_%i.arguments' % n]
            eps_infty = r['classmat.component_%i.eps_infty' % n]
            eps = r['classmat.component_%i.eps' % n] # Data as array
            eval('self.add_component(%s(permittivity = Permittivity(data=eps), %s))' % (name, arguments))
            
    # Save information on the structures for restarting
    def write(self, writer):
        w = writer
        ind = 0
        for component in self.components:
            w['classmat.component_%i.name' % ind] = component.name
            w['classmat.component_%i.arguments' % ind] = component.arguments
            w['classmat.component_%i.eps_infty' % ind] = component.permittivity.eps_infty
            w['classmat.component_%i.eps' % ind] = component.permittivity.data_eVA()
            ind += 1
        
            
    # Here the 3D-arrays are filled with material-specific information
    def apply_mask(self, mask, permittivity):
        for j in range(permittivity.Nj):
            self.bar_omega[j] = np.logical_not(mask) * self.bar_omega[j] + mask * permittivity.oscillators[j].bar_omega
            self.alpha    [j] = np.logical_not(mask) * self.alpha[j]    + mask * permittivity.oscillators[j].alpha
            self.beta     [j] = np.logical_not(mask) * self.beta[j]     + mask * permittivity.oscillators[j].beta
        
        # Add dummy oscillators if needed
        for j in range(permittivity.Nj, self.Nj):
            self.bar_omega  [j] = np.logical_not(mask) * self.bar_omega[j] + mask * 1.0
            self.alpha      [j] = np.logical_not(mask) * self.alpha[j]    + mask * 0.0
            self.beta       [j] = np.logical_not(mask) * self.beta[j]     + mask * 0.0

        # Print the permittivity information
        self.messages.append("  Permittivity data:")
        self.messages.append("    bar_omega         alpha          beta")
        self.messages.append("  ----------------------------------------")
        for j in range(permittivity.Nj):
            self.messages.append("%12.6f  %12.6f  %12.6f" % (permittivity.oscillators[j].bar_omega,
                                                             permittivity.oscillators[j].alpha,
                                                             permittivity.oscillators[j].beta))
        self.messages.append("  ----------------------------------------")
        self.messages.append("...done initializing Polarizable Material")
        masksum  = self.gd.comm.sum(int(np.sum(mask)))
        masksize = self.gd.comm.sum(int(np.size(mask)))
        self.messages.append("Fill ratio: %f percent" % (100.0 * float(masksum)/float(masksize)))
        
        
    # E(r) = -Grad V(r)
    # NB: Here -V(r) is used (if/when sign=-1), because in GPAW the
    #     electron unit charge is +1 so that the calculated V(r) is
    #     positive around negative charge. In order to get the correct
    #     direction of the electric field, the sign must be changed.
    def solve_electric_field(self, phi):
        for v in range(3):
            Gradient(self.gd, v, n=3).apply(-1.0 * self.sign * phi, self.electric_field[v])

    # n(r) = -Div P(r)
    def solve_rho(self):
        self.charge_density *= 0.0
        dmy         = self.gd.empty()
        for v in range(3):
            Gradient(self.gd, v, n=3).apply(self.polarization_total[v], dmy)
            self.charge_density -= dmy

    # P(r, omega) = [eps(r, omega) - eps0] E(r, omega)
    # P0(r) = [eps_inf(r) - eps0] E0(r) + sum_j P0_j(r) // Gao2012, Eq. 10
    def solve_polarizations(self):
        for v in range(3):
            self.polarizations[v] = _eps0_au * self.beta / (self.bar_omega**2.0) * self.electric_field[v]
        self.polarization_total = np.sum(self.polarizations, axis=1) + (self.eps_infty - _eps0_au ) * self.electric_field

    def propagate_polarizations(self, timestep):
        for v in range(3):
            self.polarizations[v] = self.polarizations[v] + timestep * self.currents[v]
        self.polarization_total = np.sum(self.polarizations, axis=1)
    
    def propagate_currents(self, timestep):
        c1 = (1.0 - 0.5 * self.alpha*timestep)/(1.0 + 0.5 * self.alpha*timestep)
        c2 = - timestep / (1.0 + 0.5 * self.alpha*timestep) * (self.bar_omega**2.0)
        c3 = - timestep / (1.0 + 0.5 * self.alpha*timestep) * (-1.0) * _eps0_au * self.beta
        for v in range(3):
            self.currents[v] = c1 * self.currents[v] + c2 * self.polarizations[v] + c3 * self.electric_field[v]

    def kick_electric_field(self, timestep, kick):
        for v in range(3):
            self.electric_field[v] = self.electric_field[v] + kick[v] / timestep

# Box-shaped classical material
class PolarizableBox():
    def __init__(self, corner1, corner2, permittivity):
        # sanity check
        assert(len(corner1)==3)
        assert(len(corner2)==3)

        self.corner1      = np.array(corner1)/Bohr # from Angstroms to atomic units
        self.corner2      = np.array(corner2)/Bohr # from Angstroms to atomic units
        self.permittivity = permittivity
        
        self.name = 'PolarizableBox'
        self.arguments = 'corner1=[%f, %f, %f], corner2=[%f, %f, %f]' % (corner1[0], corner1[1], corner1[2],
                                                                         corner2[0], corner2[1], corner2[2])        

    # Setup grid descriptor and the permittivity values inside the box
    def get_mask(self, gd):
        
        # 3D coordinates at each grid point
        r_gv = gd.get_grid_point_coordinates().transpose((1, 2, 3, 0))
        
        # inside or outside
        return np.logical_and(np.logical_and( # z
                np.logical_and(np.logical_and( # y
                 np.logical_and( #x
                                r_gv[:, :, :, 0] > self.corner1[0],
                                r_gv[:, :, :, 0] < self.corner2[0]),
                                 r_gv[:, :, :, 1] > self.corner1[1]),
                                 r_gv[:, :, :, 1] < self.corner2[1]),
                                  r_gv[:, :, :, 2] > self.corner1[2]),
                                  r_gv[:, :, :, 2] < self.corner2[2])


# Shape from atom positions (surrounding region)
class PolarizableAtomisticRegion():
    def __init__(self, atoms=None, atom_positions=None, distance=0.0, permittivity=None):

        if atoms is not None:
            self.atom_positions = np.array(atoms.get_positions())
        else:
            self.atom_positions = np.array(atom_positions)
        
        # sanity check
        assert(len(self.atom_positions)>1)
        
        self.permittivity = permittivity
        self.name = 'PolarizableAtomisticRegion'
        
        # use the minimum interatomic distance
        if distance<1e-10:
            self.distance = np.sqrt(np.min([np.sort(np.sum((self.atom_positions-ap)**2, axis=1))[1] for ap in self.atom_positions]))/Bohr
        else:
            self.distance   = distance/Bohr # from Angstroms to atomic units
       
        dbohr = self.distance*Bohr
        self.arguments = 'distance = %20.12e, atom_positions=[' % dbohr
        for ap in self.atom_positions:
            self.arguments += '[%20.12e, %20.12e, %20.12e],' % (ap[0], ap[1], ap[2])
        self.arguments = self.arguments[:-1] + ']'

    # Setup grid descriptor and the permittivity values inside the box
    def get_mask(self, gd):
        
        # 3D coordinates at each grid point
        r_gv = gd.get_grid_point_coordinates().transpose((1, 2, 3, 0))
        
        # inside or outside
        mask_tbl = False * np.ones(r_gv.shape[:-1])
        for ap in self.atom_positions:
            mask_tbl = np.logical_or(mask_tbl, np.array( (r_gv[:, :, :, 0] - ap[0]/Bohr)**2.0 +
                                                         (r_gv[:, :, :, 1] - ap[1]/Bohr)**2.0 +
                                                         (r_gv[:, :, :, 2] - ap[2]/Bohr)**2.0 < self.distance**2))
        return mask_tbl
    
# Sphere-shaped classical material
class PolarizableSphere():
    def __init__(self, center, radius, permittivity):
        self.permittivity = permittivity
        self.center      = np.array(center)/Bohr # from Angstroms to atomic units
        self.radius      = radius/Bohr # from Angstroms to atomic units
        
        # sanity check
        assert(len(self.center)==3)
        
        self.name = 'PolarizableSphere'
        self.arguments = 'center=[%20.12e, %20.12e, %20.12e], radius=%20.12e' % (center[0], center[1], center[2], radius)

    def get_mask(self, gd):
        
        # 3D coordinates at each grid point
        r_gv = gd.get_grid_point_coordinates().transpose((1, 2, 3, 0))
        
        # inside or outside
        return  np.array( (r_gv[:, :, :, 0] - self.center[0])**2.0 +
                          (r_gv[:, :, :, 1] - self.center[1])**2.0 +
                          (r_gv[:, :, :, 2] - self.center[2])**2.0 < self.radius**2)

# Sphere-shaped classical material
class PolarizableEllipsoid():
    def __init__(self, center, radii, permittivity):
        # sanity check
        assert(len(center)==3)
        assert(len(radii)==3)
        
        self.center       = np.array(center)/Bohr # from Angstroms to atomic units
        self.radii        = np.array(radii)/Bohr   # from Angstroms to atomic units
        self.permittivity = permittivity
        
        self.name = 'PolarizableEllipsoid'
        self.arguments = 'center=[%20.12e, %20.12e, %20.12e], radii=[%20.12e, %20.12e, %20.12e]' % (center[0], center[1], center[2], radii[0], radii[1], radii[2])
        

    def get_mask(self, gd):
        
        # 3D coordinates at each grid point
        r_gv = gd.get_grid_point_coordinates().transpose((1, 2, 3, 0))
        
        # inside or outside
        return  np.array( (r_gv[:, :, :, 0] - self.center[0])**2.0/self.radii[0]**2.0 +
                          (r_gv[:, :, :, 1] - self.center[1])**2.0/self.radii[1]**2.0 +
                          (r_gv[:, :, :, 2] - self.center[2])**2.0/self.radii[2]**2.0 < 1.0)

 # Rod-shaped classical material
class PolarizableRod():
    def __init__(self, corners, radius, permittivity, round_corners=True):
        # sanity check
        assert(np.array(corners).shape[0]>1)  # at least two points
        assert(np.array(corners).shape[1]==3) # 3D
        
        self.name = 'PolarizableRod'
        self.arguments = 'radius = %20.12e, corners=[' % radius
        for c in corners:
            self.arguments += '[%20.12e, %20.12e, %20.12e],' % (c[0], c[1], c[2])
        self.arguments = self.arguments[:-1] + ']'
        
        self.corners      = np.array(corners)/Bohr # from Angstroms to atomic units
        self.radius       = radius/Bohr  # from Angstroms to atomic units
        self.round_corners = round_corners
        self.permittivity = permittivity

    def get_mask(self, gd):
        
        # 3D coordinates at each grid point
        r_gv = gd.get_grid_point_coordinates().transpose((1, 2, 3, 0))
        ng = r_gv.shape[0:-1]
        ngv = r_gv.shape
        
        a = self.corners[0]

        mask = False * np.ones(ng)

        for p in self.corners[1:]:
            #http://en.wikipedia.org/wiki/Distance_from_a_point_to_a_line:
            # d = |(a-p)-((a-p).n)n|   point p, line a+tn  (|n|=1)
            n = (p-a)/np.sqrt((p-a).dot((p-a)))
            v1 = np.array([a[w]-r_gv[:, :, :, w] for w in range(3)]).transpose((1, 2, 3, 0)) # a-p

            v2 = np.sum(np.array([v1[:, :, :, w]*n[w] for w in range(3)]), axis=0)  # (a-p).n

            v3 = np.array([v2*n[w] for w in range(3)]).transpose(1, 2, 3, 0)       # ((a-p).n)n

            d = np.zeros(ng)
            for ind, idx in np.ndenumerate(v3[:, :, :, 0]):
                d[ind] = np.array(v1[ind]-v3[ind]).dot(v1[ind]-v3[ind])

            # angle between (p-a) and (r-a):
            pa = p-a # (3)
            ra = np.array([r_gv[:, :, :, w]-a[w] for w in range(3)]).transpose((1, 2, 3, 0)) # (ng1, ng2, ng3, 3)
            para = np.sum([pa[w]*ra[:, :, :, w] for w in range(3)], axis=0)
            ll2   = pa.dot(pa)*np.sum([ra[:, :, :, w]*ra[:, :, :, w] for w in range(3)], axis=0)
            angle1 = np.arccos(para/(1.0e-9+np.sqrt(ll2)))
            
            # angle between (a-p) and (r-p):
            ap = a-p # (3)
            rp = np.array([r_gv[:, :, :, w]-p[w] for w in range(3)]).transpose((1, 2, 3, 0)) # (ng1, ng2, ng3, 3)
            aprp = np.sum([ap[w]*rp[:, :, :, w] for w in range(3)], axis=0)
            ll2   = ap.dot(ap)*np.sum([rp[:, :, :, w]*rp[:, :, :, w] for w in range(3)], axis=0)
            angle2 = np.arccos(aprp/(1.0e-9+np.sqrt(ll2)))

            # Include in the mask
            this_mask = np.logical_and(np.logical_and(angle1 < 0.5*np.pi, angle2 < 0.5*np.pi),
                                      d < self.radius**2.0 )

            # Add spheres around current end points 
            if self.round_corners:
                # |r-a| and |r-p|
                raDist = np.sum([ra[:, :, :, w]*ra[:, :, :, w] for w in range(3)], axis=0)
                rpDist = np.sum([rp[:, :, :, w]*rp[:, :, :, w] for w in range(3)], axis=0)
                this_mask = np.logical_or(this_mask,
                                         np.logical_or(raDist < self.radius**2.0, rpDist < self.radius**2.0))

            mask =  np.logical_or(mask, this_mask)

            # move to next point
            a = p

        return mask

class PolarizableTetrahedron():
    #http://steve.hollasch.net/cgindex/geometry/ptintet.html
    #obrecht@imagen.com (Doug Obrecht) writes:
    #
    # Can someone point me to an algorithm that determines if a point is within a tetrahedron?
    #
    # Let the tetrahedron have vertices
    #     V1 = (x1, y1, z1)
    #    V2 = (x2, y2, z2)
    #    V3 = (x3, y3, z3)
    #    V4 = (x4, y4, z4)
    #
    #and your test point be
    #
    #        P = (x, y, z).
    #Then the point P is in the tetrahedron if following five determinants all have the same sign.
    #
    #             |x1 y1 z1 1|
    #        D0 = |x2 y2 z2 1|
    #             |x3 y3 z3 1|
    #             |x4 y4 z4 1|
    #
    #             |x  y  z  1|
    #        D1 = |x2 y2 z2 1|
    #             |x3 y3 z3 1|
    #             |x4 y4 z4 1|
    #
    #             |x1 y1 z1 1|
    #        D2 = |x  y  z  1|
    #             |x3 y3 z3 1|
    #             |x4 y4 z4 1|
    # 
    #             |x1 y1 z1 1|
    #        D3 = |x2 y2 z2 1|
    #             |x  y  z  1|
    #             |x4 y4 z4 1|
    #
    #             |x1 y1 z1 1|
    #        D4 = |x2 y2 z2 1|
    #             |x3 y3 z3 1|
    #             |x  y  z  1|
    #
    # Some additional notes:
    #
    # If by chance the D0=0, then your tetrahedron is degenerate (the points are coplanar).
    # If any other Di=0, then P lies on boundary i (boundary i being that boundary formed by the three points other than Vi).
    # If the sign of any Di differs from that of D0 then P is outside boundary i.
    # If the sign of any Di equals that of D0 then P is inside boundary i.
    # If P is inside all 4 boundaries, then it is inside the tetrahedron.
    # As a check, it must be that D0 = D1+D2+D3+D4.
    # The pattern here should be clear; the computations can be extended to simplicies of any dimension. (The 2D and 3D case are the triangle and the tetrahedron).
    # If it is meaningful to you, the quantities bi = Di/D0 are the usual barycentric coordinates.
    # Comparing signs of Di and D0 is only a check that P and Vi are on the same side of boundary i.

    def __init__(self, corners, permittivity):
        # sanity check
        assert(len(corners)==4)     # exactly 4 points
        assert(len(corners[0])==3)  # 3D
        
        self.name = 'PolarizableTetrahedron'
        self.arguments = 'corners=['
        for c in corners:
            self.arguments += '[%20.12e, %20.12e, %20.12e],' % (c[0], c[1], c[2])
        self.arguments += ']'
        
        self.corners      = np.array(corners)/Bohr # from Angstroms to atomic units
        self.permittivity = permittivity

    def determinant_value(self, x, y, z,
                       x1, y1, z1,
                       x2, y2, z2,
                       x3, y3, z3,
                       x4, y4, z4, ind):
        mat = np.array([[x1, y1, z1, 1], [x2, y2, z2, 1], [x3, y3, z3, 1], [x4, y4, z4, 1]])
        mat[ind][:] = np.array([x, y, z, 1])
        return np.linalg.det(mat)

    def get_mask(self, gd):
        r_gv = gd.get_grid_point_coordinates().transpose((1, 2, 3, 0))
        ng   = r_gv.shape[0:-1]
        ngv  = r_gv.shape
        x1, y1, z1 = self.corners[0]
        x2, y2, z2 = self.corners[1]
        x3, y3, z3 = self.corners[2]
        x4, y4, z4 = self.corners[3]
        
        mask = np.ones(ng)==0
        
        # TODO: associate a determinant for each point, and use numpy tools to determine
        #       the mask without the *very slow* loop over grid points
        for ind, pt in np.ndenumerate(r_gv[:, :, :, 0]):
            x, y, z = r_gv[ind][:]
            d0 = np.array([[x1, y1, z1, 1],
                           [x2, y2, z2, 1],
                           [x3, y3, z3, 1],
                           [x4, y4, z4, 1]])
            d1 = np.array([[x,  y,  z, 1],
                           [x2, y2, z2, 1],
                           [x3, y3, z3, 1],
                           [x4, y4, z4, 1]])
            d2 = np.array([[x1, y1, z1, 1],
                           [x,  y,  z,  1],
                           [x3, y3, z3, 1],
                           [x4, y4, z4, 1]])
            d3 = np.array([[x1, y1, z1, 1],
                           [x2, y2, z2, 1],
                           [x,  y,  z,  1],
                           [x4, y4, z4, 1]])
            d4 = np.array([[x1, y1, z1, 1],
                           [x2, y2, z2, 1],
                           [x3, y3, z3, 1],
                           [x,  y,  z,  1]])
            s0 = np.linalg.det(d0)
            s1 = np.linalg.det(d1)
            s2 = np.linalg.det(d2)
            s3 = np.linalg.det(d3)
            s4 = np.linalg.det(d4)
            
            if (np.sign(s0) == np.sign(s1) or abs(s1) < 1e-12) and \
               (np.sign(s0) == np.sign(s2) or abs(s2) < 1e-12) and \
               (np.sign(s0) == np.sign(s3) or abs(s3) < 1e-12) and \
               (np.sign(s0) == np.sign(s4) or abs(s4) < 1e-12):
                mask[ind] = True
        return mask
       
# Lorentzian oscillator function: L(omega) = eps0 * beta / (w**2 - i*alpha*omega - omega**2)    // Coomar2011, Eq. 2
class LorentzOscillator:
    def __init__(self, bar_omega, alpha, beta):
        self.bar_omega = bar_omega
        self.alpha    = alpha
        self.beta     = beta

    def value(self, omega):
        return _eps0_au * self.beta / (self.bar_omega**2 - 1J * self.alpha * omega - omega**2)

# Dieletric function: e(omega) = eps_inf + sum_j L_j(omega) // Coomar2011, Eq. 2
class Permittivity:
    def __init__(self, fname=None, data=None, eps_infty = _eps0_au ):

        # Initialize to vacuum permittivity
        self.eps_infty = eps_infty
        self.Nj = 0
        self.oscillators = []

        # Input as data array
        if data is not None:
            assert fname is None
            self.Nj = len(data)
            for v in range(self.Nj):
                self.oscillators.append(LorentzOscillator(data[v][0]/Hartree, data[v][1]/Hartree, data[v][2]/Hartree/Hartree))
            return

        # Input as filename
        if fname != None: # read permittivity from a 3-column file
            fp = open(fname, 'r')
            lines = fp.readlines()
            fp.close()

            self.Nj = len(lines)

            for line in lines:
                bar_omega = float(split(line)[0]) / Hartree
                alpha     = float(split(line)[1]) / Hartree
                beta      = float(split(line)[2]) / Hartree / Hartree
                self.oscillators.append(LorentzOscillator(bar_omega, alpha, beta))
            return

            

    def value(self, omega = 0):
        return self.eps_infty + sum([osc.value(omega) for osc in self.oscillators])
    
    def data(self):
        return [[osc.bar_omega, osc.alpha, osc.beta] for osc in self.oscillators]
    
    def data_eVA(self):
        return [[osc.bar_omega*Hartree, osc.alpha*Hartree, osc.beta*Hartree*Hartree] for osc in self.oscillators]

# Dieletric function that renormalizes the static permittivity to the requested value (usually epsZero) 
class PermittivityPlus(Permittivity):
    def __init__(self, fname=None, data=None, eps_infty = _eps0_au, epsZero = _eps0_au, newbar_omega = 0.01, new_alpha = 0.10, **kwargs):
        Permittivity.__init__(self, fname=fname, data=data, eps_infty=eps_infty)
                
        # Convert given values from eVs to Hartrees
        _newbar_omega = newbar_omega / Hartree
        _new_alpha    = new_alpha / Hartree
        
        # Evaluate the new value    
        _new_beta = ((epsZero - self.value(0.0))*_newbar_omega**2.0/_eps0_au).real
        self.oscillators.append(LorentzOscillator(_newbar_omega, _new_alpha, _new_beta))
        self.Nj = len(self.oscillators)
        
