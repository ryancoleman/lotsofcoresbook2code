
import numpy as np

from gpaw.utilities import fact
from gpaw.sphere import lmfact
from gpaw.sphere.legendre import ilegendre, legendre, dlegendre

# Define the Heaviside function
heaviside = lambda x: (1.0+np.sign(x))/2.0

# Define spherical harmoncics and normalization coefficient
C = lambda l,m: (-1)**((m+abs(m))//2)*((2.*l+1.)/(4*np.pi*lmfact(l,m)))**0.5
Y = lambda l,m,theta,phi: C(l,m)*legendre(l,abs(m),np.cos(theta))*np.exp(1j*m*phi)

# Define theta-derivative of spherical harmoncics
dYdtheta = lambda l,m,theta,phi: C(l,m)*dlegendre(l,abs(m),np.cos(theta))*np.exp(1j*m*phi)

# Define phi-derivative of spherical harmoncics
dYdphi = lambda l,m,theta,phi: 1j*m*C(l,m)*legendre(l,abs(m),np.cos(theta))*np.exp(1j*m*phi)

# -------------------------------------------------------------------

def intYY(l1, m1, l2, m2):
    """Calculates::

                  pi 2pi
                  /  /      *
           A    = |  |    Y (u,v) Y (u,v) sin(u) dv du
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

def intYY_ex(l1, m1, l2, m2):
    """Calculates::

                  pi 2pi
                  /  /      *
           A    = |  |    Y (u,v) sin(u)*cos(v) Y (u,v) sin(u) dv du
            LL'   /  /     lm                    l'm'
                  0  0

    where u = theta and v = phi in the usual notation. Note that the result
    is only non-zero if `|l1-l2|` = 1 and `|m1-m2|` = 1.
    """
    if abs(l1-l2) != 1 or abs(m1-m2) != 1:
        return 0.0
    scale = C(l1,m1)*C(l2,m2)*np.pi
    if abs(m1) == abs(m2)+1:
        return scale*np.sign(l1-l2)/(2.*l2+1.)*ilegendre(l1,m1)
    else:
        assert abs(m1) == abs(m2)-1
        return scale*np.sign(l2-l1)/(2.*l1+1.)*ilegendre(l2,m2)

def intYY_ey(l1, m1, l2, m2):
    """Calculates::

                  pi 2pi
                  /  /      *
           A    = |  |    Y (u,v) sin(u)*sin(v) Y (u,v) sin(u) dv du
            LL'   /  /     lm                    l'm'
                  0  0

    where u = theta and v = phi in the usual notation. Note that the result
    is only non-zero if `|l1-l2|` = 1 and `|m1-m2|` = 1.
    """
    if abs(l1-l2) != 1 or abs(m1-m2) != 1:
        return 0.0
    scale = C(l1,m1)*C(l2,m2)*np.sign(m1-m2)*np.pi/1j
    if abs(m1) == abs(m2)+1:
        return scale*np.sign(l1-l2)/(2.*l2+1.)*ilegendre(l1,m1)
    else:
        assert abs(m1) == abs(m2)-1
        return scale*np.sign(l2-l1)/(2.*l1+1.)*ilegendre(l2,m2)

def intYY_ez(l1, m1, l2, m2):
    """Calculates::

                  pi 2pi
                  /  /      *
           A    = |  |    Y (u,v) cos(u) Y (u,v) sin(u) dv du
            LL'   /  /     lm             l'm'
                  0  0

    where u = theta and v = phi in the usual notation. Note that the result
    is only non-zero if `|l1-l2|` = 1 and m1 = m2.
    """
    if abs(l1-l2) != 1 or m1 != m2:
        return 0.0
    scale = C(l1,m1)*C(l2,m2)*2*np.pi
    return scale*(l2-np.sign(l1-l2)*abs(m1)+heaviside(l1-l2))/(2.*l2+1.)*ilegendre(l1,m1)

# -------------------------------------------------------------------

def intYdYdtheta_ex(l1, m1, l2, m2):
    """Calculates::

                  pi 2pi
                  /  /      *                    d Y(u,v)
           A    = |  |    Y (u,v) cos(u)*cos(v) --- l'm'  sin(u) dv du
            LL'   /  /     lm                    du  
                  0  0

    where u = theta and v = phi in the usual notation. Note that the result
    is only non-zero if `|l1-l2|` is odd and `|m1-m2|` = 1 (stricter rule applies).
    """
    if abs(l1-l2) % 2 != 1 or abs(m1-m2) != 1:
        return 0.0
    scale = -C(l1,m1)*C(l2,m2)*np.pi
    if abs(m1) == abs(m2)+1:
        if l1+1 == l2:
            return scale*2/(2.0*l2+1)*(l2+1)/(2.0*l2-1.0)*fact(l2+abs(m2))/fact(l2-abs(m2)-1)*(l2-abs(m2)-1)
        elif l1-1 == l2:
            return scale*2/(2.0*l2+1)*((l2+1)*fact(l2+abs(m2))/fact(l2-abs(m2))*(l2-abs(m2))-l2/(2.0*l2+3.0)*fact(l2+abs(m2)+1)/fact(l2-abs(m2))*(l2-abs(m2)+1))
        elif l1-l2 > 2: # and (l1-l2)%2 == 1 which is always true
            return -scale*2*abs(m2)*fact(l2+abs(m2))/fact(l2-abs(m2))
        else:
            return 0.0
    else:
        assert abs(m1) == abs(m2)-1
        if l1 == l2+1:
            return -scale*2/(2.0*l1+1.0)*(l1-1)/(2.0*l1-1.0)*fact(l1+abs(m1))/fact(l1-abs(m1)-1)*(l1-abs(m1)-1)
        elif l1 == l2-1:
            return -scale*2/(2.0*l1+1.0)*((l1+1)/(2.0*l1+3.0)*fact(l1+abs(m1)+1)/fact(l1-abs(m1))*(l1-abs(m1)+1)-(abs(m1)+1)*fact(l1+abs(m1))/fact(l1-abs(m1))*(l1-abs(m1)))
        elif l2-l1 > 2: # and (l2-l1)%2 == 1 which is always true
            return scale*2*(abs(m1)+1)*fact(l1+abs(m1))/fact(l1-abs(m1))
        else:
            return 0.0

def intYdYdtheta_ey(l1, m1, l2, m2):
    """Calculates::

                  pi 2pi
                  /  /      *                    d Y(u,v)
           A    = |  |    Y (u,v) cos(u)*sin(v) --- l'm'  sin(u) dv du
            LL'   /  /     lm                    du  
                  0  0

    where u = theta and v = phi in the usual notation. Note that the result
    is only non-zero if `|l1-l2|` is odd and `|m1-m2|` = 1 (stricter rule applies).
    """
    if abs(l1-l2) % 2 != 1 or abs(m1-m2) != 1:
        return 0.0
    scale = -C(l1,m1)*C(l2,m2)*np.sign(m1-m2)*np.pi/1j
    if abs(m1) == abs(m2)+1:
        if l1+1 == l2:
            return scale*2/(2.0*l2+1)*(l2+1)/(2.0*l2-1.0)*fact(l2+abs(m2))/fact(l2-abs(m2)-1)*(l2-abs(m2)-1)
        elif l1-1 == l2:
            return scale*2/(2.0*l2+1)*((l2+1)*fact(l2+abs(m2))/fact(l2-abs(m2))*(l2-abs(m2))-l2/(2.0*l2+3.0)*fact(l2+abs(m2)+1)/fact(l2-abs(m2))*(l2-abs(m2)+1))
        elif l1-l2 > 2: # and (l1-l2)%2 == 1 which is always true
            return -scale*2*abs(m2)*fact(l2+abs(m2))/fact(l2-abs(m2))
        else:
            return 0.0
    else:
        assert abs(m1) == abs(m2)-1
        if l1 == l2+1:
            return -scale*2/(2.0*l1+1.0)*(l1-1)/(2.0*l1-1.0)*fact(l1+abs(m1))/fact(l1-abs(m1)-1)*(l1-abs(m1)-1)
        elif l1 == l2-1:
            return -scale*2/(2.0*l1+1.0)*((l1+1)/(2.0*l1+3.0)*fact(l1+abs(m1)+1)/fact(l1-abs(m1))*(l1-abs(m1)+1)-(abs(m1)+1)*fact(l1+abs(m1))/fact(l1-abs(m1))*(l1-abs(m1)))
        elif l2-l1 > 2: # and (l2-l1)%2 == 1 which is always true
            return scale*2*(abs(m1)+1)*fact(l1+abs(m1))/fact(l1-abs(m1))
        else:
            return 0.0

def intYdYdtheta_ez(l1, m1, l2, m2):
    """Calculates::

                     pi 2pi
                     /  /      *             d Y(u,v)
           A    =  - |  |    Y (u,v) sin(u) --- l'm'  sin(u) dv du
            LL'      /  /     lm             du  
                     0  0

    where u = theta and v = phi in the usual notation. Note that the result
    is only non-zero if `|l1-l2|` = 1 and m1 = m2.
    """
    if abs(l1-l2) != 1 or m1 != m2:
        return 0.0
    scale = -C(l1,m1)*C(l2,m2)*2*np.pi
    return scale*np.sign(l1-l2)*(l2+1-heaviside(l1-l2))*(l2-np.sign(l1-l2)*abs(m1)+heaviside(l1-l2))/(2*l2+1)*ilegendre(l1,m1)

# -------------------------------------------------------------------

def intYdYdphi_ex(l1, m1, l2, m2):
    """Calculates::

                     pi 2pi
                     /  /      *       -1           d Y(u,v)
           A    =  - |  |    Y (u,v) sin(u)*sin(v) --- l'm'  sin(u) dv du
            LL'      /  /     lm                    dv  
                     0  0

    where u = theta and v = phi in the usual notation. Note that the result
    is only non-zero if `|l1-l2|` is odd and `|m1-m2|` = 1 (stricter rule applies).
    """
    if abs(l1-l2) % 2 != 1 or abs(m1-m2) != 1:
        return 0.0
    scale = -C(l1,m1)*C(l2,m2)*np.sign(m1-m2)*np.pi/1j*1j*m2
    if abs(m1) == abs(m2)+1 and l1 > l2: # and (l1-l2)%2 == 1 which is always true
        return scale*2*lmfact(l2,m2)
    elif abs(m1) == abs(m2)-1 and l1 < l2: # and (l2-l1)%2 == 1 which is always true
        return scale*2*lmfact(l1,m1)
    else:
        return 0.0

def intYdYdphi_ey(l1, m1, l2, m2):
    """Calculates::

                  pi 2pi
                  /  /      *       -1           d Y(u,v)
           A    = |  |    Y (u,v) sin(u)*cos(v) --- l'm'  sin(u) dv du
            LL'   /  /     lm                    dv  
                  0  0

    where u = theta and v = phi in the usual notation. Note that the result
    is only non-zero if `|l1-l2|` is odd and `|m1-m2|` = 1 (stricter rule applies).
    """
    if abs(l1-l2) % 2 != 1 or abs(m1-m2) != 1:
        return 0.0
    scale = C(l1,m1)*C(l2,m2)*np.pi*1j*m2
    if abs(m1) == abs(m2)+1 and l1 > l2: # and (l1-l2)%2 == 1 which is always true
        return scale*2*lmfact(l2,m2)
    elif abs(m1) == abs(m2)-1 and l1 < l2: # and (l2-l1)%2 == 1 which is always true
        return scale*2*lmfact(l1,m1)
    else:
        return 0.0

def intYdYdphi_ez(l1, m1, l2, m2): #XXX this is silly
    return 0.0

