
import numpy as np

from gpaw.utilities import fact
from gpaw.sphere import lmfact
from gpaw.sphere.legendre import ilegendre, legendre, dlegendre

# Define the Heaviside function
heaviside = lambda x: (1.0+np.sign(x))/2.0

# Define spherical harmoncics and normalization coefficient
C = lambda l,m: ((2.*l+1.)/(4*np.pi*lmfact(l,m)))**0.5

def Y(l,m,theta,phi):
    if m == 0:
        return C(l,m)*legendre(l,abs(m),np.cos(theta))
    elif m > 0:
        return C(l,m)*legendre(l,abs(m),np.cos(theta))*np.cos(m*phi)*2**0.5
    else:
        return C(l,m)*legendre(l,abs(m),np.cos(theta))*np.sin(abs(m)*phi)*2**0.5

# Define theta-derivative of spherical harmoncics
def dYdtheta(l,m,theta,phi):
    if m == 0:
        return C(l,m)*dlegendre(l,abs(m),np.cos(theta))
    elif m > 0:
        return C(l,m)*dlegendre(l,abs(m),np.cos(theta))*np.cos(m*phi)*2**0.5
    else:
        return C(l,m)*dlegendre(l,abs(m),np.cos(theta))*np.sin(abs(m)*phi)*2**0.5

# Define phi-derivative of spherical harmoncics
def dYdphi(l,m,theta,phi):
    if m == 0:
        return np.zeros_like(theta)
    elif m > 0:
        return -m*C(l,m)*legendre(l,abs(m),np.cos(theta))*np.sin(m*phi)*2**0.5
    else:
        return -m*C(l,m)*legendre(l,abs(m),np.cos(theta))*np.cos(m*phi)*2**0.5

# -------------------------------------------------------------------

def intYY(l1, m1, l2, m2):
    """Calculates::

                  pi 2pi
                  /  /
           A    = |  |    y (u,v) y (u,v) sin(u) dv du
            LL'   /  /     lm      l'm'
                  0  0

    where u = theta and v = phi in the usual notation. Note that the result
    is only non-zero if l1 = l2 and m1 = m2.
    """
    if (l1,m1) == (l2,m2):
        return C(l1,m1)**2*2*np.pi*ilegendre(l1,m1) # == 1 always, XXX right?
    else:
        return 0.0

# -------------------------------------------------------------------

from gpaw.sphere.csh import intYY_ex as _intYY_ex, intYY_ey as _intYY_ey, intYY_ez as _intYY_ez #TODO make independent

def mix(l1, m1, l2, m2, func):
    # m1 == 0
    if m1 == 0 and m2 == 0:
        return func(l1,abs(m1),l2,abs(m2))
    elif m1 == 0 and m2 > 0:
        return 1/2**0.5*((-1)**m2*func(l1,abs(m1),l2,abs(m2)) \
                        +func(l1,abs(m1),l2,-abs(m2)))
    elif m1 == 0 and m2 < 0:
        return 1/(2**0.5*1j)*((-1)**m2*func(l1,abs(m1),l2,abs(m2)) \
                             -func(l1,abs(m1),l2,-abs(m2)))
    # m1 > 0
    elif m1 > 0 and m2 == 0:
        return 1/2**0.5*((-1)**m1*func(l1,abs(m1),l2,abs(m2)) \
                        +func(l1,-abs(m1),l2,abs(m2)))
    elif m1 > 0 and m2 > 0:
        return 0.5*((-1)**(m1+m2)*func(l1,abs(m1),l2,abs(m2)) \
                   +(-1)**m1*func(l1,abs(m1),l2,-abs(m2)) \
                   +(-1)**m2*func(l1,-abs(m1),l2,abs(m2)) \
                   +func(l1,-abs(m1),l2,-abs(m2)))
    elif m1 > 0 and m2 < 0:
        return 1/2j*((-1)**(m1+m2)*func(l1,abs(m1),l2,abs(m2)) \
                    -(-1)**m1*func(l1,abs(m1),l2,-abs(m2)) \
                    +(-1)**m2*func(l1,-abs(m1),l2,abs(m2)) \
                    -func(l1,-abs(m1),l2,-abs(m2)))
    # m1 < 0
    elif m1 < 0 and m2 == 0:
        return -1/(2**0.5*1j)*((-1)**m1*func(l1,abs(m1),l2,abs(m2)) \
                             -func(l1,-abs(m1),l2,abs(m2)))
    elif m1 < 0 and m2 > 0:
        return -1/2j*((-1)**(m1+m2)*func(l1,abs(m1),l2,abs(m2)) \
                     +(-1)**m1*func(l1,abs(m1),l2,-abs(m2)) \
                     -(-1)**m2*func(l1,-abs(m1),l2,abs(m2)) \
                     -func(l1,-abs(m1),l2,-abs(m2)))
    elif m1 < 0 and m2 < 0:
        return 0.5*((-1)**(m1+m2)*func(l1,abs(m1),l2,abs(m2)) \
                   -(-1)**m1*func(l1,abs(m1),l2,-abs(m2)) \
                   -(-1)**m2*func(l1,-abs(m1),l2,abs(m2)) \
                   +func(l1,-abs(m1),l2,-abs(m2)))
    else:
        raise ValueError('Invalid arguments (l1=%s, m1=%s, l2=%s, m2=%s)' \
                         % (l1,m1,l2,m2))

def intYY_ex(l1, m1, l2, m2):
    """Calculates::

                  pi 2pi
                  /  /
           A    = |  |    y (u,v) sin(u)*cos(v) y (u,v) sin(u) dv du
            LL'   /  /     lm                    l'm'
                  0  0

    where u = theta and v = phi in the usual notation. Note that the result
    is only non-zero if `|l1-l2|` = 1 and `|m1-m2|` = 1.
    """
    #if abs(l1-l2) != 1 or abs(m1-m2) != 1:
    #    return 0.0
    #scale = C(l1,m1)*C(l2,m2)*np.pi
    #if abs(m1) == abs(m2)+1:
    #    return scale*np.sign(l1-l2)/(2.*l2+1.)*ilegendre(l1,m1)
    #else:
    #    assert abs(m1) == abs(m2)-1
    #    return scale*np.sign(l2-l1)/(2.*l1+1.)*ilegendre(l2,m2)
    return mix(l1,m1,l2,m2,_intYY_ex)

def intYY_ey(l1, m1, l2, m2):
    """Calculates::

                  pi 2pi
                  /  /
           A    = |  |    y (u,v) sin(u)*sin(v) y (u,v) sin(u) dv du
            LL'   /  /     lm                    l'm'
                  0  0

    where u = theta and v = phi in the usual notation. Note that the result
    is only non-zero if `|l1-l2|` = 1 and `|m1-m2|` = 1.
    """
    #if abs(l1-l2) != 1 or abs(m1-m2) != 1:
    #    return 0.0
    #scale = C(l1,m1)*C(l2,m2)*np.sign(m1-m2)*np.pi/1j
    #if abs(m1) == abs(m2)+1:
    #    return scale*np.sign(l1-l2)/(2.*l2+1.)*ilegendre(l1,m1)
    #else:
    #    assert abs(m1) == abs(m2)-1
    #    return scale*np.sign(l2-l1)/(2.*l1+1.)*ilegendre(l2,m2)
    return mix(l1,m1,l2,m2,_intYY_ey)

def intYY_ez(l1, m1, l2, m2):
    """Calculates::

                  pi 2pi
                  /  /
           A    = |  |    y (u,v) cos(u) y (u,v) sin(u) dv du
            LL'   /  /     lm             l'm'
                  0  0

    where u = theta and v = phi in the usual notation. Note that the result
    is only non-zero if `|l1-l2|` = 1 and m1 = m2.
    """
    return mix(l1,m1,l2,m2,_intYY_ez)

# -------------------------------------------------------------------

from gpaw.sphere.csh import intYdYdtheta_ex as _intYdYdtheta_ex, intYdYdtheta_ey as _intYdYdtheta_ey, intYdYdtheta_ez as _intYdYdtheta_ez #TODO make independent

def intYdYdtheta_ex(l1, m1, l2, m2):
    """Calculates::

                  pi 2pi
                  /  /                           d y(u,v)
           A    = |  |    y (u,v) cos(u)*cos(v) --- l'm'  sin(u) dv du
            LL'   /  /     lm                    du  
                  0  0

    where u = theta and v = phi in the usual notation. Note that the result
    is only non-zero if `|l1-l2|` is odd and `|m1-m2|` = 1 (stricter rule applies).
    """
    #if abs(l1-l2) % 2 != 1 or abs(m1-m2) != 1:
    #    return 0.0
    #scale = -C(l1,m1)*C(l2,m2)*np.pi
    #if abs(m1) == abs(m2)+1:
    #    if l1+1 == l2:
    #        return scale*2/(2.0*l2+1)*(l2+1)/(2.0*l2-1.0)*fact(l2+abs(m2))/fact(l2-abs(m2)-1)*(l2-abs(m2)-1)
    #    elif l1-1 == l2:
    #        return scale*2/(2.0*l2+1)*((l2+1)*fact(l2+abs(m2))/fact(l2-abs(m2))*(l2-abs(m2))-l2/(2.0*l2+3.0)*fact(l2+abs(m2)+1)/fact(l2-abs(m2))*(l2-abs(m2)+1))
    #    elif l1-l2 > 2: # and (l1-l2)%2 == 1 which is always true
    #        return -scale*2*abs(m2)*fact(l2+abs(m2))/fact(l2-abs(m2))
    #    else:
    #        return 0.0
    #else:
    #    assert abs(m1) == abs(m2)-1
    #    if l1 == l2+1:
    #        return -scale*2/(2.0*l1+1.0)*(l1-1)/(2.0*l1-1.0)*fact(l1+abs(m1))/fact(l1-abs(m1)-1)*(l1-abs(m1)-1)
    #    elif l1 == l2-1:
    #        return -scale*2/(2.0*l1+1.0)*((l1+1)/(2.0*l1+3.0)*fact(l1+abs(m1)+1)/fact(l1-abs(m1))*(l1-abs(m1)+1)-(abs(m1)+1)*fact(l1+abs(m1))/fact(l1-abs(m1))*(l1-abs(m1)))
    #    elif l2-l1 > 2: # and (l2-l1)%2 == 1 which is always true
    #        return scale*2*(abs(m1)+1)*fact(l1+abs(m1))/fact(l1-abs(m1))
    #    else:
    #        return 0.0
    return mix(l1,m1,l2,m2,_intYdYdtheta_ex)

def intYdYdtheta_ey(l1, m1, l2, m2):
    """Calculates::

                  pi 2pi
                  /  /                           d y(u,v)
           A    = |  |    y (u,v) cos(u)*sin(v) --- l'm'  sin(u) dv du
            LL'   /  /     lm                    du  
                  0  0

    where u = theta and v = phi in the usual notation. Note that the result
    is only non-zero if `|l1-l2|` is odd and `|m1-m2|` = 1 (stricter rule applies).
    """
    #if abs(l1-l2) % 2 != 1 or abs(m1-m2) != 1:
    #    return 0.0
    #scale = -C(l1,m1)*C(l2,m2)*np.sign(m1-m2)*np.pi/1j
    #if abs(m1) == abs(m2)+1:
    #    if l1+1 == l2:
    #        return scale*2/(2.0*l2+1)*(l2+1)/(2.0*l2-1.0)*fact(l2+abs(m2))/fact(l2-abs(m2)-1)*(l2-abs(m2)-1)
    #    elif l1-1 == l2:
    #        return scale*2/(2.0*l2+1)*((l2+1)*fact(l2+abs(m2))/fact(l2-abs(m2))*(l2-abs(m2))-l2/(2.0*l2+3.0)*fact(l2+abs(m2)+1)/fact(l2-abs(m2))*(l2-abs(m2)+1))
    #    elif l1-l2 > 2: # and (l1-l2)%2 == 1 which is always true
    #        return -scale*2*abs(m2)*fact(l2+abs(m2))/fact(l2-abs(m2))
    #    else:
    #        return 0.0
    #else:
    #    assert abs(m1) == abs(m2)-1
    #    if l1 == l2+1:
    #        return -scale*2/(2.0*l1+1.0)*(l1-1)/(2.0*l1-1.0)*fact(l1+abs(m1))/fact(l1-abs(m1)-1)*(l1-abs(m1)-1)
    #    elif l1 == l2-1:
    #        return -scale*2/(2.0*l1+1.0)*((l1+1)/(2.0*l1+3.0)*fact(l1+abs(m1)+1)/fact(l1-abs(m1))*(l1-abs(m1)+1)-(abs(m1)+1)*fact(l1+abs(m1))/fact(l1-abs(m1))*(l1-abs(m1)))
    #    elif l2-l1 > 2: # and (l2-l1)%2 == 1 which is always true
    #        return scale*2*(abs(m1)+1)*fact(l1+abs(m1))/fact(l1-abs(m1))
    #    else:
    #        return 0.0
    return mix(l1,m1,l2,m2,_intYdYdtheta_ey)

def intYdYdtheta_ez(l1, m1, l2, m2):
    """Calculates::

                     pi 2pi
                     /  /                    d y(u,v)
           A    =  - |  |    y (u,v) sin(u) --- l'm'  sin(u) dv du
            LL'      /  /     lm             du  
                     0  0

    where u = theta and v = phi in the usual notation. Note that the result
    is only non-zero if `|l1-l2|` = 1 and m1 = m2.
    """
    #if abs(l1-l2) != 1 or m1 != m2:
    #    return 0.0
    #scale = -C(l1,m1)*C(l2,m2)*2*np.pi
    #return scale*np.sign(l1-l2)*(l2+1-heaviside(l1-l2))*(l2-np.sign(l1-l2)*abs(m1)+heaviside(l1-l2))/(2*l2+1)*ilegendre(l1,m1)
    return mix(l1,m1,l2,m2,_intYdYdtheta_ez)

# -------------------------------------------------------------------

from gpaw.sphere.csh import intYdYdphi_ex as _intYdYdphi_ex, intYdYdphi_ey as _intYdYdphi_ey #TODO make independent

def intYdYdphi_ex(l1, m1, l2, m2):
    """Calculates::

                     pi 2pi
                     /  /              -1           d y(u,v)
           A    =  - |  |    y (u,v) sin(u)*sin(v) --- l'm'  sin(u) dv du
            LL'      /  /     lm                    dv  
                     0  0

    where u = theta and v = phi in the usual notation. Note that the result
    is only non-zero if `|l1-l2|` is odd and `|m1-m2|` = 1 (stricter rule applies).
    """
    #if abs(l1-l2) % 2 != 1 or abs(m1-m2) != 1:
    #    return 0.0
    #scale = -C(l1,m1)*C(l2,m2)*np.sign(m1-m2)*np.pi/1j*1j*m2
    #if abs(m1) == abs(m2)+1 and l1 > l2: # and (l1-l2)%2 == 1 which is always true
    #    return scale*2*lmfact(l2,m2)
    #elif abs(m1) == abs(m2)-1 and l1 < l2: # and (l2-l1)%2 == 1 which is always true
    #    return scale*2*lmfact(l1,m1)
    #else:
    #    return 0.0
    return mix(l1,m1,l2,m2,_intYdYdphi_ex)

def intYdYdphi_ey(l1, m1, l2, m2):
    """Calculates::

                  pi 2pi
                  /  /              -1           d y(u,v)
           A    = |  |    y (u,v) sin(u)*cos(v) --- l'm'  sin(u) dv du
            LL'   /  /     lm                    dv  
                  0  0

    where u = theta and v = phi in the usual notation. Note that the result
    is only non-zero if `|l1-l2|` is odd and `|m1-m2|` = 1 (stricter rule applies).
    """
    #if abs(l1-l2) % 2 != 1 or abs(m1-m2) != 1:
    #    return 0.0
    #scale = C(l1,m1)*C(l2,m2)*np.pi*1j*m2
    #if abs(m1) == abs(m2)+1 and l1 > l2: # and (l1-l2)%2 == 1 which is always true
    #    return scale*2*lmfact(l2,m2)
    #elif abs(m1) == abs(m2)-1 and l1 < l2: # and (l2-l1)%2 == 1 which is always true
    #    return scale*2*lmfact(l1,m1)
    #else:
    #    return 0.0
    return mix(l1,m1,l2,m2,_intYdYdphi_ey)

def intYdYdphi_ez(l1, m1, l2, m2): #XXX this is silly
    return 0.0

# -------------------------------------------------------------------

def intYgradY(l1, m1, l2, m2, r_g, dr_g, A_g, B_g, dBdr_g, v=None):
    """Calculates::

                  pi 2pi
                  /  /    *           / __               \  2
          A     = |  |   A(r) y (u,v) | \/  B(r) y (u,v) | r sin(u) dv du dr
           vLL'   /  /         lm     \   v       l'm'   /
                  0  0

    where u = theta and v = phi in the usual notation. Note that the result
    is only non-zero if `|l1-l2|` is odd and `|m1-m2|` <= 1 (stricter rule applies).
    """
    if v is None:
        D = [intYY_ex(l1,m1,l2,m2), \
             intYY_ey(l1,m1,l2,m2), \
             intYY_ez(l1,m1,l2,m2)]
        G = [intYdYdtheta_ex(l1,m1,l2,m2) + intYdYdphi_ex(l1,m1,l2,m2), \
             intYdYdtheta_ey(l1,m1,l2,m2) + intYdYdphi_ey(l1,m1,l2,m2), \
             intYdYdtheta_ez(l1,m1,l2,m2) + intYdYdphi_ez(l1,m1,l2,m2)]
        D, G = np.array(D), np.array(G)
    elif v == 0:
        D = intYY_ex(l1,m1,l2,m2)
        G = intYdYdtheta_ex(l1,m1,l2,m2) + intYdYdphi_ex(l1,m1,l2,m2)
    elif v == 1:
        D = intYY_ey(l1,m1,l2,m2)
        G = intYdYdtheta_ey(l1,m1,l2,m2) + intYdYdphi_ey(l1,m1,l2,m2)
    elif v == 2:
        D = intYY_ez(l1,m1,l2,m2)
        G = intYdYdtheta_ez(l1,m1,l2,m2) + intYdYdphi_ez(l1,m1,l2,m2)
    else:
        raise ValueError
    return D * np.vdot(A_g, dBdr_g * r_g**2 * dr_g) \
        + G * np.vdot(A_g, B_g * r_g * dr_g)

