import numpy as np
from Numeric import pi, sqrt
from tools import factorial
from tools import Rational as Q

"""
This is a script designed for construction of the real solid spherical
harmonics (RSSH) in cartesian form. These can be written as::

         m    |m|  l  |m|
  Y  =  Y  = C    r  P  (cos theta) Phi (phi)
   L     l    l       l                m

where C_l^|m| is a normalization constant
P_l^|m| is the associatied legendre polynomial
and:

              / cos(m phi) , m > 0
  Phi (phi) = |   1        , m = 0
     m        \ sin(-m phi), m < 0

The first few harmonics are listed below::
 +----+---------------------+-__---------------------------------------+
 |  L | l |  m | r^l * Y    | \/ (r^l * Y)                             |
 +----+---s----+------------+------------------------------------------+
 |  0 | 0 |  0 | 1          | (0, 0, 0)                                |
 +----+---p----+------------+------------------------------------------+
 |  1 | 1 | -1 | y          | (0, 1, 0)                                |
 |  2 | 1 |  0 | z          | (0, 0, 1)                                |
 |  3 | 1 |  1 | x          | (1, 0, 0)                                |
 +----+---d----+------------+------------------------------------------+
 |  4 | 2 | -2 | xy         | ( y,  x,  0)                             |
 |  5 | 2 | -1 | yz         | ( 0,  z,  y)                             |
 |  6 | 2 |  0 | 3z^2-r^2   | (-x, -y, 2z)                             |
 |  7 | 2 |  1 | xz         | ( z,  0,  x)                             |
 |  8 | 2 |  2 | x^2-y^2    | ( x, -y,  0)                             |
 +----+---f----+------------+------------------------------------------+
 |  9 | 3 | -3 | 3x^2y-y^3  | (          2xy,       x^2-y^2,        0) |
 | 10 | 3 | -2 | xyz        | (           yz,            xz,       xy) |
 | 11 | 3 | -1 | 5yz^2-yr^2 | (         -2xy, 4z^2-x^2-3y^2,      8yz) |
 | 12 | 3 |  0 | 5z^3-3zr^2 | (         -2xz,          -2yz, 3z^2-r^2) |
 | 13 | 3 |  1 | 5xz^2-xr^2 | (4z^2-3x^2-y^2,          -2xy,      8xz) |
 | 14 | 3 |  2 | x^2z-y^2z  | (          2xz,          -2yz,  x^2-y^2) |
 | 15 | 3 |  3 | x^3-3xy^2  | (      x^2-y^2,          -2xy,        0) |
 +----+--------+----------+--------------------------------------------+

Y_lm is represented as a polynomial in x, y, and z

The function consists of three parts: a normalization constant accessed by
class 'Normalization(l, m)', a polynomial in z accessed with method
'legendre(l, m)', and a polynomial in x and y accessed with method 'Phi(l, m)'

The normalization and the z-polynomial are both invariant of the sign of m
The z-polynomial has powers l-|m|, l-|m|-2, l-|m|-4, l-..., i.e. it is strictly odd (even) if l-|m| is odd (even)
The combined power of x and y is |m| in all terms of Phi
"""

Y_lp = [{}, {}] # Global list of dictionaries for storing calculated
                # Legendre polynomials, and Phi functions

#--------------------------- RELEVANT USER METHODS ---------------------------
def L_to_lm(L):
    """convert L index to (l, m) index"""
    l = int(sqrt(L))
    m = L - l**2 - l
    return l, m

def lm_to_L(l,m):
    """convert (l, m) index to L index"""
    return l**2 + l + m

def Y_to_string(l, m, deriv=None, multiply=None, numeric=False):
    # for L in range(40): print L, Y_to_string(*L_to_lm(L))
    """                                                   l    m
       If deriv is None, return string representation of r  * Y (x, y, z)
                                                               l
       
       If deriv == q, return string is the derivative of above with respect
       to x, y or z if q is 0, 1 or 2 respectively.

       multiply=q indicates that the entire expression should be multiplied by
       x, y or z if q is 0, 1 or 2 respectively.

       numeric=True/False indicates whether the normalization constant should
       be written as a numeric or an algebraic expression.
    """
    assert deriv is None or deriv in range(3)
    assert multiply is None or multiply in range(3)
    
    if deriv is None:
        norm, xyzs = Y_collect(l, m)
    else:
        norm, xyzs = dYdq(l, m, deriv)

    if multiply is not None:
        xyzs = q_times_xyzs(xyzs, multiply)

    string = to_string(l, xyzs, deriv is not None, multiply is not None)
    if string == '0': return '0'
    else: return norm.tostring(numeric) + (' * ' + string) * (string != '1')

def gauss_to_string(l, m, numeric=False):
    """Return string representation of the generalized gaussian::
    
                       _____                           2  
         m            /  1       l!         l+3/2  -a r   l  m
        g (x,y,z) =  / ----- --------- (4 a)      e      r  Y (x,y,z)
         l         \/  4 pi  (2l + 1)!                       l
         
       numeric=True/False indicates whether the normalization constant should 
       be written as a number or an algebraic expression.
    """
    norm, xyzs = Y_collect(l, m)

    ng = Q(2**(2*l+3) * factorial(l), 2 * factorial(2 * l + 1))
    norm.multiply(ng)

    string = to_string(l, xyzs)
    string = (' * ' + string) * (string != '1')
    if numeric:
        snorm = repr(eval(repr(norm.norm)))
    else:
        snorm = repr(norm.norm)
    string = 'sqrt(a**%s*%s)/pi'%(2*l+3, snorm) + string
    string += ' * exp(-a*r2)'

    return string

def gauss_potential_to_string(l, m, numeric=False):
    """Return string representation of the potential of  a generalized
       gaussian.

       The potential is determined by::

          m        m ^    _           m ^
         v [g (r) Y (r) ](r) = v (r) Y (r)
          l  l     l         l     l  l

       where::
               4 pi /  -l-1 /r    l+2         l /oo   1-l      \ 
       v (r) = ---- | r     | dx x   g (r) + r  | dx x   g (r) |
        l      2l+1 \       /0        l         /r        l    /
    """            
    v_l = [[Q(4,1), 1],
           [Q(4,3), 1, 2],
           [Q(4,15), 3, 6, 4],
           [Q(4,105), 15, 30, 20, 8],
           [Q(4,945), 105, 210, 140, 56, 16],
           [Q(4,10395), 945, 1890, 1260, 504, 144, 32],
           ]

    norm, xyzs = Y_collect(l, m)
    norm.multiply(v_l[l][0])

    string = txt_sqrt(norm.norm, numeric) + '*' + (l!=0)*'('
    if numeric:
        string += repr(v_l[l][1] * sqrt(pi))
    else:
        string += str(v_l[l][1]) + '*sqrt(pi)'
    string += '*erf(sqrt(a)*r)'

    if len(v_l[l]) > 2:
        string += '-('
        for n, coeff in enumerate(v_l[l][2:]):
            if n == 0:
                string += str(coeff)
            else:
                string += '+' + str(coeff) + '*(sqrt(a)*r)**%d'%(2*n)
        string += ')*sqrt(a)*r*exp(-a*r2)'

    if l == 0:
        string += '/r'
    elif l == 1:
        string += ')/r/r2*' + to_string(l, xyzs)
    else:
        string += ')/r/r2**%d*'%l + to_string(l, xyzs)

    return string

#----------------------------- TECHNICAL METHODS -----------------------------
def to_string(l, xyzs, deriv=False, multiply=False):
    """Return string representation of an xyz dictionary"""
    if xyzs == {}: return '0'
    out = ''

    for xyz, coef in xyzs.items():
        x, y, z = xyz
        r = l - x - y - z - deriv + multiply
        one = abs(coef) != 1 or (x == 0 and y == 0 and z == 0 and r == 0)
        out += sign(coef) + str(abs(coef)) * one
        out += ('*x'*x + '*y'*y + '*z'*z + '*r2'*(r/2))[1 - one:]

    if out[0] == '+': out = out[1:]
    if len(xyzs) > 1: out = '(' + out + ')'
    return out

def sign(x):
    """Return string representation of the sign of x"""
    if x >= 0: return '+'
    else: return '-'

def txt_sqrt(norm, numeric=False):
    if numeric:
        return repr(sqrt(norm))
    else:
        if sqrt(norm) % 1 == 0:
            return str(sqrt(norm))
        else:
            return 'sqrt(' + str(norm.nom) + \
                   ('./' + str(norm.denom)) * (norm.denom != 1) + ')'

class Normalization:
    """Determine normalization factor of spherical harmonic
                   ______________
             /    / 2l+1   (l-m)!
             |   /  ---- * ------  , m != 0
             | \/   2 pi   (l+m)!
       C  = <      _____
        L    |    / 2l+1
             |   /  ----           , m = 0
             \ \/   4 pi
    """
    def __init__(self, l, m):
        m = abs(m)
        if m == 0:
            self.norm = Q(2 * l + 1, 4)
        else:
            self.norm = Q((2 * l + 1) * factorial(l - m), 2 * factorial(l + m))

    def __str__(self):
        n = self.norm
        sn = sqrt(n)
        if int(sn) == sn:
            string = repr(sn) + '/sqrt(pi)'
        else:
            string = 'sqrt(' + repr(n.nom) + \
                     ('./' + repr(n.denom)) * (n.denom != 1) + '/pi)'
        return string

    def __repr__(self):
        return repr(self.__float__())

    def __float__(self):
        return sqrt(self.norm / pi)

    def multiply(self, x):
        self.norm *= x**2

    def tostring(self, numeric=False):
        if numeric:
            return self.__repr__()
        else:
            return self.__str__()

def legendre(l, m):
    """Determine z dependence of spherical harmonic.
       Returns vector, where the p'th element is the coefficient of
       z^p r^(l-|m|-p).
    """
    # Check if requested has already been calculated
    if (l, m) in Y_lp[0]:
        return Y_lp[0][(l, m)]
    
    m = abs(m)
    assert l >= 0 and 0 <= m <=l
    result = np.zeros(l - m + 1, 'O')
    if l == m == 0:
        """Use that
             0
            P (z) = 1
             0
        """
        result[0] = Q(1)
    elif l == m:
        """Use the recursion relation
            m              m-1
           P (z) = (2m-1) P   (z)
            m              m-1
        """
        result[:] += (2 * m - 1) * legendre(l - 1, m - 1)
    elif l == m + 1:
        """Use the recursion relation
            l-1              l-1
           P  (z) = (2l-1)z P   (z)
            l                l-1
            
        """
        result[1:] += (2 * l - 1) * legendre(l-1, l-1)
    else:
        """Use the recursion relation
            m     2l-1    m       l+m-1  2  m
           P (z)= ---- z P  (z) - ----- r  P  (z)
            l      l-m    l-1      l-m      l-2
        """
        result[1:] += np.multiply(legendre(l - 1, m), Q(2 * l - 1, l - m))
        result[:(l - 2) - m + 1] -= np.multiply(legendre(l - 2, m),
                                                 Q(l + m - 1, l - m))
    # Store result in global dictionary
    Y_lp[0][(l, m)] = result
    return result

def Phi(m):
    """Determine the x and y dependence of the spherical harmonics from
                      |m|   |m|
                   / r   sin  (theta) cos(|m| phi), m >= 0
       Phi (phi) = |
          m        |  |m|   |m|
                   \ r   sin  (theta) sin(|m| phi), m < 0
       Returns dictionary of format {(i, j): c} where c is the coefficient
       of x^i y^j
    """
    # Check if requested has already been calculated
    if m in Y_lp[1]:
        return Y_lp[1][m]
    
    if   m ==  0:
        xys = {(0, 0): 1} # use that Phi_0  = 1
    elif m ==  1:
        xys = {(1, 0): 1} # use that Phi_1  = x
    elif m == -1:
        xys = {(0, 1): 1} # use that Phi_-1 = y
    else:
        """Use the recurrence formula
        
           m > 0:  Phi (x,y) = x Phi   (x,y) - y Phi   (x,y)
                     |m|           |m|-1           1-|m|

           m < 0:  Phi (x,y) = y Phi   (x,y) + x Phi   (x,y)
                     |m|           |m|-1           1-|m|           
        """
        xys  = {}
        phi1 = Phi(abs(m) - 1)
        phi2 = Phi(1 - abs(m))
        for x, y in phi1:
            new = (x + (m > 0), y + (m < 0))
            xys[new] = xys.get(new, 0) +  phi1[(x, y)]
        for x,y in phi2:
            new = (x + (m < 0), y + (m > 0))
            sign = 2 * (m < 0) - 1
            xys[new] = xys.get(new, 0) + sign * phi2[(x, y)]

    # Store result in global dictionary
    Y_lp[1][m] = xys
    return xys

def Y_collect(l, m):
    """Collect all necessary parts of spherical harmonic and return in
       simplified format.
       Return dictionary xyzs has format {(i, j, k): c} where c is the
       coefficient of x^i y^j z^k r^(l-|m|-k), or (since i+j = |m|) the
       coefficient of x^i y^j z^k r^(l-i-j-k), from which it is clear that all
       terms are of power l in x, y and z collectively.
    """
    zs = legendre(l, m)
    xys = Phi(m)

    xyzs = {}
    for xy in xys:
        if xys[xy] != 0:
            for p in range(len(zs)):
                if zs[p] != 0:
                    xyzs[xy + (p,)] = xys[xy] * zs[p]

    # get normalization constant and simplify
    norm = Normalization(l, m)
    norm.multiply(simplify(xyzs))
    
    return norm, xyzs

def Y_collect2(l, m):
    """Same as Y_collect, but collective power of x, y, and z are
    adjusted, such the it is always equal to l (thus avoiding
    multiplication by r)
    """
    norm, p = Y_collect(l, m)
    done = False
    while not done:
        p2 = {}
        done = True
        for (nx, ny, nz), c in p.items():
            n = nx + ny + nz
            if n < l:
                p2[(nx + 2, ny, nz)] = p2.get((nx + 2, ny, nz), 0) + c
                p2[(nx, ny + 2, nz)] = p2.get((nx, ny + 2, nz), 0) + c
                p2[(nx, ny, nz + 2)] = p2.get((nx, ny, nz + 2), 0) + c
                if n + 2 < l:
                    done = False
            else:
                assert n == l
                p2[(nx, ny, nz)] = p2.get((nx, ny, nz), 0) + c
        p = p2
    p2 = p.copy()
    for n, c in p.items():
        if c == 0:
            del p2[n]
    return norm, p2

def dYdq(l, m, q):
    """Returns a normalization constant, and a dictionary discribing
       the functional form of the derivative of r^l Y_l^m(x,y,z) with
       respect to x, y or z if q is either 0, 1 or 2 respectively. The
       format of the output dictionary is {(i, j, k): c}, where c is the
       coefficient of x^i y^j z^k r^(l-i-j-k-1).
    """
    norm, xyzs = Y_collect(l, m)
    dxyzs = {}
    
    for xyz, coef in xyzs.items():
        x, y, z = xyz
        r = l - x - y - z

        # chain rule: diff coordinate q only
        if xyz[q] != 0:
            dxyz = list(xyz)
            dxyz[q] -= 1
            dxyz = tuple(dxyz)
            dxyzs[dxyz] = dxyzs.get(dxyz, 0) + xyz[q] * coef
        
        # chain rule: diff coordinate r only
        if r != 0:
            dxyz = list(xyz)
            dxyz[q] += 1
            dxyz = tuple(dxyz)
            dxyzs[dxyz] = dxyzs.get(dxyz, 0) + r * coef

    # remove zeros from list
    for dxyz in dxyzs.keys():
        if dxyzs[dxyz] == 0: dxyzs.pop(dxyz)

    # simplify
    if dxyzs != {}: norm.multiply(simplify(dxyzs))
    return norm, dxyzs

def simplify(xyzs):
    """Rescale coefficients to smallest integer value"""
    norm = Q(1)
    numxyz = np.array(xyzs.values())

    # up-scale all 'xyz' coefficients to integers
    for xyz in numxyz:
        numxyz *= xyz.denom
        norm /= xyz.denom

    # determine least common divisor for 'xyz' coefficients
    dmax = 1
    num_max = max(abs(np.floor(numxyz)))
    for d in range(2, num_max + 1):
        test = numxyz / d
        if np.alltrue(test == np.floor(test)): dmax = d

    # Update simplified dictionary
    norm *= dmax
    for i, xyz in enumerate(xyzs):
        xyzs[xyz] = numxyz[i] / dmax

    return norm

def q_times_xyzs(xyzs, q):
    """multiply xyz dictionary by x, y, or z according to q = 0, 1, or 2"""
    qxyzs = {}
    for xyz, c in xyzs.items():
        qxyz = list(xyz)
        qxyz[q] += 1
        qxyz = tuple(qxyz)
        
        qxyzs[qxyz] = c
    return qxyzs

#--------------------- TEST AND CODE CONSTRUCTING METHODS ---------------------
def orthogonal(L1, L2):
    """Perform the integral
          2pi pi
           /  /
       I = |  |sin(theta) d(theta) d(phi) Y (theta, phi) * Y (theta, phi)
           /  /                            L1               L2
           0  0
       which should be a kronecker delta in L1 and L2
    """
    I = 0.0
    N = 40

    for theta in np.arange(0, pi, pi / N):
        for phi in np.arange(0, 2 * pi, 2 * pi / N):
            x = np.cos(phi) * np.sin(theta)
            y = np.sin(phi) * np.sin(theta)
            z = np.cos(theta)
            r2 = x*x + y*y + z*z
            
            Y1 = eval(Y_to_string(*L_to_lm(L1)))
            Y2 = eval(Y_to_string(*L_to_lm(L2)))

            I += np.sin(theta) * Y1 * Y2
    I *= 2 * (pi / N)**2

    return I

def check_orthogonality(Lmax=10):
    """Check orthogonality for all combinations of the first few harmonics"""
    all_passed = True
    for L1 in range(Lmax+1):
        for L2 in range(L1, Lmax+1):
            I = orthogonal(L1, L2)
            passed =  abs(I - (L1 == L2)) < 3e-3
            all_passed *= passed
            print('L1 = %s,  L2 = %s, passed = %s, I = %s' %(L1, L2, passed, I))
    if all_passed: print('All tests passed')
    else: print('Some tests failed')

def symmetry1(lmax, display=True):
    """Make dictionary of format
       diff = {(l1, m1, q1): (nrel, l2, m2, q2)}
       indicating that
            m1              m2
         d Y             d Y
            l1              l2
         ------ = nrel * ------
          d q1            d q2       
    """
    diff = {} # diff[(l1, m1, q1)] = (nrel, l2, m2, q2)
    unique_L = [] # unique_L[L] = (l, m, q, norm, dxyzs)
    for L in range((lmax + 1)**2):
        l, m = L_to_lm(L)
        for q in range(3):
            identical = False
            name = (l, m, 'xyz'[q])

            norm, dxyzs = dYdq(l, m, q)

            for unique in unique_L:
                if dxyzs == unique[4]:
                    diff[name] = (norm.eval() / unique[3],) + unique[0:3]
                    identical = True
                    break
            if identical == False:
                unique_L.append(name + (norm.eval(), dxyzs))
    if display:
        for key, value in diff.items():
            print(str(key) + ' = ' + str(value[0]) + ' * ' + str(value[1:]))
    else: return diff

def symmetry2(l, display=True):
    """Make dictionary of format
       diff = {(l1, m1, q1): (nrel, l2, m2, q2)}
       indicating that
              m1              m2
           d Y             d Y
              l1              l2
           ------ = nrel * ------
            d q1            d q2
       and
                m1                m2
          q1 * Y   = nrel * q2 * Y
                l1                l2
    """
    diff = {} # diff[(l1, m1, q1)] = (nrel, l2, m2, q2)
    unique_L = [] # unique_L[L] = (l, m, q, dnorm, dxyzs, qnorm, qxyzs)
    for m in range(-l, l+1):
        for q in range(3):
            identical = False
            name = (l, m, q)

            qnorm, xyzs = Y_collect(l, m)
            qxyzs = q_times_xyzs(xyzs, q)
            dnorm, dxyzs = dYdq(l, m, q)
            
            for unique in unique_L:
                if dxyzs == unique[4] and qxyzs == unique[6]:
                    dnrel = dnorm.eval() / unique[3]
                    qnrel = qnorm.eval() / unique[5]
                    print(dnrel == qnrel)
                    if dnrel == qnrel:
                        diff[name] = (dnrel,) + unique[0:3]
                        identical = True
                        break
            if identical == False:
                unique_L.append(name + (dnorm.eval(), dxyzs,
                                        qnorm.eval(), qxyzs))
    if display:
        for key, value in diff.items():
            print(str(key) + ' = ' + str(value[0]) + ' * ' + str(value[1:]))
    else: return diff

def construct_spherical_harmonic_c_function(file, lmax, funcname,
                                            multiply=None, deriv=None):
    """Construct a macro for evaluating values of spherical harmonics,
    or the derivative of any spherical harmonic with respect to some axis.
    The deriv keyword corresponds to that of the Y_to_string function."""
    w = file.write
    indent = 0
    def wn(string=''):
        w(2 * indent * ' ')
        w(string)
        w('\\\n')
    wn('#define %s(l, f, x, y, z, r2, p) (' % funcname)
    indent = 2
    wn('{')
    wn('  switch(l)')
    wn('    {')
    switchindent = 3
    indent += switchindent
    for l in range(lmax + 1):
        wn('case %d:' % l)
        indent += 1
        for M, m in enumerate(range(-l, l + 1)):
            Ystr = Y_to_string(l, m, numeric=True, deriv=deriv)
            wn('p[%d] = f * %s;' % (M, Ystr))
        wn('break;')
        indent -= 1
    wn('default:')
    wn('  assert(0 == 1);')
    indent -= switchindent
    wn('    }')
    wn('}')
    indent = 0
    wn(')')
    w('\n')


def construct_spherical_harmonic_c_code(filename='spherical_harmonics.h',
                                        lmax=4):
    """Construct macros for evaluating spherical harmonics as well as their
    derivatives."""
    file = open(filename, 'w')
    construct = construct_spherical_harmonic_c_function
    construct(file, lmax, 'spherical_harmonics')
    for c in range(3):
        construct(file, lmax, 'spherical_harmonics_derivative_%s' % 'xyz'[c],
                  multiply=c, deriv=c)
    file.close()

    
def construct_c_code(file='temp.c', lmax=3):
    """Method for generating the code in c/spline.c"""
    txt = '//Computer generated code! Hands off!'
    start_func = """
    
// inserts values of f(r) r^l Y_lm(theta, phi) in elements of input array 'a'
void bmgs_radial3(const bmgsspline* spline, int m, 
                  const int n[3], 
                  const double C[3],
                  const double h[3],
                  const double* f, double* a)
{
  int l = spline->l;
  if (l == 0)
    for (int q = 0; q < n[0] * n[1] * n[2]; q++)
      a[q] = 0.28209479177387814 * f[q];
"""
    start_deriv = """

// insert values of
// d( f(r) * r^l Y_l^m )                           d( r^l Y_l^m )
// --------------------- = g(r) q r^l Y_l^m + f(r) --------------
//        dq                                             dq
// where q={x, y, z} and g(r) = 1/r*(df/dr)
void bmgs_radiald3(const bmgsspline* spline, int m, int c, 
                  const int n[3], 
                  const double C[3],
                  const double h[3],
                  const double* f, const double* g, double* a)
{
  int l = spline->l;
"""
    start_case = """
    {
      int q = 0;
      double x = C[0];
      for (int i0 = 0; i0 < n[0]; i0++)
        {
          double y = C[1];
          for (int i1 = 0; i1 < n[1]; i1++)
            {
              double z = C[2];
              for (int i2 = 0; i2 < n[2]; i2++, q++)
                {
"""
    end_case = """
                  z += h[2];
                }
              y += h[1];
            }
          x += h[0];
        }
    }
"""
    
    # insert code for evaluating the function
    txt += start_func
    for l in range(1, lmax + 1):
        txt += '  else if (l == %s)' %l
        txt += start_case
        case = ''
        for m in range(-l, l+1):
            if m == -l: case += ' ' * 18 + 'if (m == %s)\n' %m
            elif m == l: case += '\n' + ' ' * 18 +'else\n'
            else: case += '\n' + ' ' * 18 + 'else if (m == %s)\n' %m
            case += ' ' * 20 + 'a[q] = f[q] * '
            case += Y_to_string(l,m, numeric=True) + ';'
        if 'r2' in case: txt += ' ' * 18 + 'double r2 = x*x+y*y+z*z;\n'
        txt += case
        txt += end_case
    txt += """  else
    assert(0 == 1);
}
"""
    
    # insert code for evaluating the derivative
    txt += start_deriv
    for q in range(3):
        txt += '  // ' + 'xyz'[q] + '\n'
        for l in range(0, lmax + 1):
            if l == 0 and q == 0:
                txt += '  if (c == 0 && l == 0)'
            else: txt += '  else if (c == %s && l == %s)' %(q, l)
        
            txt += start_case
            case = ''
            for m in range(-l, l+1):
                if m == -l: case += ' ' * 18 + 'if (m == %s)\n' %m
                elif m == l: case += '\n' + ' ' * 18 + 'else\n'
                else: case += '\n' + ' ' * 18 + 'else if (m == %s)\n' %m
                case += ' ' * 20 + 'a[q] = g[q] * '
                case += Y_to_string(l, m, multiply=q, numeric=True)
                diff = Y_to_string(l, m, deriv=q, numeric=True)
                if diff != '0':
                    case += ' + f[q] * ' + diff
                case += ';'
            if 'r2' in case: txt += ' ' * 18 + 'double r2 = x*x+y*y+z*z;\n'
            txt += case
            txt += end_case
    txt += """  else
      assert(0 == 1);
}
"""
    f = open(file, 'w')
    print(txt, file=f)
    f.close()
    
def construct_gauss_code(lmax=2):
    """Method for generating the code in gpaw/utilities/gauss.py"""
    Lmax = (lmax + 1)**2
    out= 'Y_L = [\n'
    for L in range(Lmax):
        l, m = L_to_lm(L)
        out+= '  \'' + Y_to_string(l, m, numeric=True) + '\',\n'
    out += ']'

    out += '\ngauss_L = [\n'
    for L in range(Lmax):
        l, m = L_to_lm(L)
        out += '  \'' + gauss_to_string(l, m, numeric=True) + '\',\n'
    out += ']'
    
    out += '\ngausspot_L = [\n'
    for L in range(Lmax):
        l, m = L_to_lm(L)
        out += '  \'' + gauss_potential_to_string(l, m, numeric=True) + '\',\n'
    out += ']'
    
    print(out)

def construct_spherical_code(lmax=3):
    """Method for generating the code in gpaw/spherical_harmonics.py"""
    YL = []
    norms = []
    for L in range((lmax+1)**2):
        #norm, xyzs = Y_collect(*L_to_lm(L))
        norm, xyzs = Y_collect2(*L_to_lm(L))
        norms.append(str(norm))
        YL.append(zip(xyzs.values(), xyzs.keys()))

    print('Y_L = [')
    for L, Y in enumerate(YL):
        l = sqrt(L)
        if l % 1 == 0:
            print('  #' + 'spdfghijklmn'[int(l)] + ':')
        print('  %s,' % Y)
    print(']')
    print('norms =', norms)


if __name__ == '__main__':
    construct_spherical_harmonic_c_code()
