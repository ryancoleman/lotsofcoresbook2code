
import numpy as np

from gpaw.utilities import erf
from gpaw.gaunt import gaunt as G_LLL
from gpaw.spherical_harmonics import Y, nablarlYL as nablaY #XXX use gpaw.sphere.rsh
from gpaw.sphere.rsh import intYY, intYY_ex, intYY_ey, intYY_ez

def bilinear_concentric_potential(r_g, dr_g, f_g, ft_g, l1, l2, alpha, rfilt=None):
    """Calculate corrections for concentric functions and potentials::

                 /     _   _a    _   _a    ~   _   _a  ~ _   _a   _
        v      = | f  (r - R ) V(r - R ) - f  (r - R ) V(r - R ) dr
         m1,m2   /  L1,L2                    L1,L2

    where f(r) and ft(r) are bilinear product of two localized functions which
    are radial splines times real spherical harmonics (l1,m1) or (l2,m2) and::

          _       1       _ -1              ~ _    erf(alpha*r)  _ -1
        V(r) = --------- |r|        ^       V(r) = ------------ |r|
               4 pi eps0                            4 pi eps0

    Note that alpha (and rfilt) should conform with the cutoff radius.
    """
    work_g = erf(alpha*r_g)

    if rfilt is None:
        M = np.vdot(f_g - ft_g * work_g, r_g * dr_g)
    else:
        M = np.vdot((f_g - ft_g * work_g)[r_g>=rfilt], \
            (r_g * dr_g)[r_g>=rfilt])

        # Replace 1/r -> (3-r^2/rfilt^2)/(2*rfilt) for r < rfilt
        M += np.vdot((f_g - ft_g * work_g)[r_g<rfilt], \
            (r_g**2/(2*rfilt) * (3-(r_g/rfilt)**2) * dr_g)[r_g<rfilt])

    v_mm = np.empty((2*l1+1, 2*l2+1), dtype=float)
    for m1 in range(2*l1+1):
        for m2 in range(2*l2+1):
            v_mm[m1,m2] = M * intYY(l1, m1-l1, l2, m2-l2)
    return v_mm

def bilinear_concentric_force(r_g, dr_g, f_g, ft_g, l1, l2, alpha, rfilt=None):
    """Calculate corrections for concentric functions and potentials::

        _        /     _   _a  __    _   _a    ~   _   _a  __  ~ _   _a   _
        F      = | f  (r - R ) \/_ V(r - R ) - f  (r - R ) \/_ V(r - R ) dr
         m1,m2   /  L1,L2        r               L1,L2       r

    where f(r) and ft(r) are bilinear product of two localized functions which
    are radial splines times real spherical harmonics (l1,m1) or (l2,m2) and::

          _       1       _ -1              ~ _    erf(alpha*r)  _ -1
        V(r) = --------- |r|        ^       V(r) = ------------ |r|
               4 pi eps0                            4 pi eps0

    Note that alpha (and rfilt) should conform with the cutoff radius.
    """
    work_g = erf(alpha*r_g) - 2*alpha/np.pi**0.5 \
        * r_g * np.exp(-alpha**2 * r_g**2)

    if rfilt is None:
        M = - np.vdot(f_g - ft_g * work_g, dr_g)
    else:
        M = - np.vdot((f_g - ft_g * work_g)[r_g>=rfilt], dr_g[r_g>=rfilt])

        # Replace 1/r -> (3-r^2/rfilt^2)/(2*rfilt) for r < rfilt
        work_g = (r_g/rfilt) * erf(alpha*r_g) - alpha*rfilt/np.pi**0.5 \
            * (3-(r_g/rfilt)**2) * np.exp(-alpha**2 * r_g**2)
        M += - np.vdot((f_g * (r_g/rfilt) - ft_g * work_g)[r_g<rfilt], \
            ((r_g/rfilt)**2 * dr_g)[r_g<rfilt])

    F_mmv = np.empty((2*l1+1, 2*l2+1, 3), dtype=float)
    for m1 in range(2*l1+1):
        for m2 in range(2*l2+1):
            F_mmv[m1,m2,0] = M * intYY_ex(l1, m1-l1, l2, m2-l2)
            F_mmv[m1,m2,1] = M * intYY_ey(l1, m1-l1, l2, m2-l2)
            F_mmv[m1,m2,2] = M * intYY_ez(l1, m1-l1, l2, m2-l2)
    return F_mmv

# -------------------------------------------------------------------

def bilinear_expansion_potential(r_g, dr_g, f_g, l1, l2, dR_v, lmax=4):
    """Calculate multipole expansion of off-center function and potential::

                 /     _   _a'   _   _a   _         _       1       _ -1
        v      = | f  (r - R ) V(r - R ) dr     , V(r) = --------- |r|
         m1,m2   /  L1,L2                                4 pi eps0

    where f(r) is a bilinear product of two localized functions which are
    radial splines times real spherical harmonics (l1,m1) or (l2,m2).

    The seperation between atoms a and a' is specified in atomic units as::

         _   _a  _a'
        dR = R - R

    """
    r = np.sum(dR_v**2)**0.5
    Lmax = (lmax + 1)**2
    v_mm = np.zeros((2*l1+1, 2*l2+1), dtype=float)
    for m1 in range(2*l1+1):
        L1 = l1**2 + m1
        for m2 in range(2*l2+1):
            L2 = l2**2 + m2
            Q_L = np.empty(Lmax, dtype=float)
            for l in range(lmax+1):
                for m in range(2*l + 1):
                    L = l**2 + m
                    Q_L[L] = np.vdot(f_g, r_g**(2+l) * dr_g) * G_LLL[L1,L2,L]
                    v_mm[m1,m2] += 4 * np.pi / ((2*l+1) * r**(2*l+1)) \
                        * Q_L[L] * Y(L, *tuple(dR_v))
    return v_mm


def bilinear_expansion_force(r_g, dr_g, f_g, l1, l2, dR_v, lmax=4):
    """Calculate multipole expansion of off-center function and potential::

        _        /     _   _a' __    _   _a   _        _       1       _ -1
        F      = | f  (r - R ) \/_ V(r - R ) dr    , V(r) = --------- |r|
         m1,m2   /  L1,L2        r                          4 pi eps0

    where f(r) is a bilinear product of two localized functions which are
    radial splines times real spherical harmonics (l1,m1) or (l2,m2).

    The seperation between atoms a and a' is specified in atomic units as::

         _   _a  _a'
        dR = R - R

    """
    r = np.sum(dR_v**2)**0.5
    Lmax = (lmax + 1)**2
    F_mmv = np.zeros((2*l1+1, 2*l2+1, 3), dtype=float)
    for m1 in range(2*l1+1):
        L1 = l1**2 + m1
        for m2 in range(2*l2+1):
            L2 = l2**2 + m2
            Q_L = np.empty(Lmax, dtype=float)
            for l in range(lmax+1):
                for m in range(2*l + 1):
                    L = l**2 + m
                    Q_L[L] = np.vdot(f_g, r_g**(2+l) * dr_g) * G_LLL[L1,L2,L]
                    F_mmv[m1,m2,:] += 4 * np.pi / r**(2*l+1) * Q_L[L] \
                        * ( dR_v / r**2 * Y(L, *tuple(dR_v)) \
                           - np.array(nablaY(L, dR_v)) / (2*l+1) )
    return F_mmv

