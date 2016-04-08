
import numpy as np

from gpaw.sphere import lmfact

# Highest l for which associated Legendre polynomials are implemented
lmax = 9

# Integral norm of the associated Legendre polynomials
ilegendre = lambda l,m: 2./(2.*l+1.)*lmfact(l,m)

def legendre(l, m, w):
    if m < 0 or l < m:
        return np.zeros_like(w)
    if (l,m) == (0,0):
        return np.ones_like(w)
    elif (l,m) == (1,0):
        return w.copy()
    elif (l,m) == (1,1):
        return (1-w**2)**0.5
    elif (l,m) == (2,0):
        return 1.5*w**2-0.5
    elif (l,m) == (2,1):
        return 3*(1-w**2)**0.5*w
    elif (l,m) == (2,2):
        return 3*(1-w**2)
    elif (l,m) == (3,0):
        return 2.5*w**3-1.5*w
    elif (l,m) == (3,1):
        return (1-w**2)**0.5*(7.5*w**2-1.5)
    elif (l,m) == (3,2):
        return 15*w-15*w**3
    elif (l,m) == (3,3):
        return 15*(1-w**2)**1.5
    elif (l,m) == (4,0):
        return 4.375*w**4-3.75*w**2+0.375
    elif (l,m) == (4,1):
        return (1-w**2)**0.5*(17.5*w**3-7.5*w)
    elif (l,m) == (4,2):
        return (1-w**2)*(52.5*w**2-7.5)
    elif (l,m) == (4,3):
        return 105*(1-w**2)**1.5*w
    elif (l,m) == (4,4):
        return 105*(1-w**2)**2
    elif (l,m) == (5,0):
        return 7.875*w**5-8.75*w**3+1.875*w
    elif (l,m) == (5,1):
        return (1-w**2)**0.5*(39.375*w**4-26.25*w**2+1.875)
    elif (l,m) == (5,2):
        return (1-w**2)*(157.5*w**3-52.5*w)
    elif (l,m) == (5,3):
        return (1-w**2)**1.5*(472.5*w**2-52.5)
    elif (l,m) == (5,4):
        return 945*(1-w**2)**2*w
    elif (l,m) == (5,5):
        return 945*(1-w**2)**2.5
    elif (l,m) == (6,0):
        return 14.4375*w**6-19.6875*w**4+6.5625*w**2-0.3125
    elif (l,m) == (6,1):
        return (1-w**2)**0.5*(86.625*w**5-78.75*w**3+13.125*w)
    elif (l,m) == (6,2):
        return (1-w**2)*(433.125*w**4-236.25*w**2+13.125)
    elif (l,m) == (6,3):
        return (1-w**2)**1.5*(1732.5*w**3-472.5*w)
    elif (l,m) == (6,4):
        return (1-w**2)**2*(5197.5*w**2-472.5)
    elif (l,m) == (6,5):
        return 10395*(1-w**2)**2.5*w
    elif (l,m) == (6,6):
        return 10395*(1-w**2)**3
    elif (l,m) == (7,0):
        return 26.8125*w**7-43.3125*w**5+19.6875*w**3-2.1875*w
    elif (l,m) == (7,1):
        return (1-w**2)**0.5*(187.6875*w**6-216.5625*w**4+59.0625*w**2-2.1875)
    elif (l,m) == (7,2):
        return (1-w**2)*(1126.125*w**5-866.25*w**3+118.125*w)
    elif (l,m) == (7,3):
        return (1-w**2)**1.5*(5630.625*w**4-2598.75*w**2+118.125)
    elif (l,m) == (7,4):
        return (1-w**2)**2*(22522.5*w**3-5197.5*w)
    elif (l,m) == (7,5):
        return (1-w**2)**2.5*(67567.5*w**2-5197.5)
    elif (l,m) == (7,6):
        return 135135*(1-w**2)**3*w
    elif (l,m) == (7,7):
        return 135135.*(1-w**2)**3.5
    elif (l,m) == (8,0):
        return 50.2734375*w**8-93.84375*w**6+54.140625*w**4-9.84375*w**2+0.2734375
    elif (l,m) == (8,1):
        return (1-w**2)**0.5*(402.1875*w**7-563.0625*w**5+216.5625*w**3-19.6875*w)
    elif (l,m) == (8,2):
        return (1-w**2)*(2815.3125*w**6-2815.3125*w**4+649.6875*w**2-19.6875)
    elif (l,m) == (8,3):
        return (1-w**2)**1.5*(16891.875*w**5-11261.25*w**3+1299.375*w)
    elif (l,m) == (8,4):
        return (1-w**2)**2*(84459.375*w**4-33783.75*w**2+1299.375)
    elif (l,m) == (8,5):
        return (1-w**2)**2.5*(337837.5*w**3-67567.5*w)
    elif (l,m) == (8,6):
        return (1-w**2)**3*(1013512.5*w**2-67567.5)
    elif (l,m) == (8,7):
        return 2027025*(1-w**2)**3.5*w
    elif (l,m) == (8,8):
        return 2027025*(1-w**2)**4
    elif (l,m) == (9,0):
        return 94.9609375*w**9-201.09375*w**7+140.765625*w**5-36.09375*w**3+2.4609375*w
    elif (l,m) == (9,1):
        return (1-w**2)**0.5*(854.6484375*w**8-1407.65625*w**6+703.828125*w**4-108.28125*w**2+2.4609375)
    elif (l,m) == (9,2):
        return (1-w**2)*(6837.1875*w**7-8445.9375*w**5+2815.3125*w**3-216.5625*w)
    elif (l,m) == (9,3):
        return (1-w**2)**1.5*(47860.3125*w**6-42229.6875*w**4+8445.9375*w**2-216.5625)
    elif (l,m) == (9,4):
        return (1-w**2)**2*(287161.875*w**5-168918.75*w**3+16891.875*w)
    elif (l,m) == (9,5):
        return (1-w**2)**2.5*(1435809.375*w**4-506756.25*w**2+16891.875)
    elif (l,m) == (9,6):
        return (1-w**2)**3*(5743237.5*w**3-1013512.5*w)
    elif (l,m) == (9,7):
        return (1-w**2)**3.5*(17229712.5*w**2-1013512.5)
    elif (l,m) == (9,8):
        return 34459425*(1-w**2)**4*w
    elif (l,m) == (9,9):
        return 34459425*(1-w**2)**4.5
    else:
        raise ValueError('Unsupported arguments (l,m)=(%d,%d)' % (l,m))


def dlegendre(l, m, w): # XXX not d/dw but rather -sin(theta)*d/dw
    if m < 0 or l < m:
        return np.zeros_like(w)
    if (l,m) == (0,0):
        return np.zeros_like(w)
    elif (l,m) == (1,0):
        return -(1-w**2)**0.5
    elif (l,m) == (1,1):
        return w.copy()
    elif (l,m) == (2,0):
        return -3*w*(1-w**2)**0.5
    elif (l,m) == (2,1):
        return 6*w**2-3
    elif (l,m) == (2,2):
        return 6*(1-w**2)**0.5*w
    elif (l,m) == (3,0):
        return (1-w**2)**0.5*(-7.5*w**2+1.5)
    elif (l,m) == (3,1):
        return 22.5*w**3-16.5*w
    elif (l,m) == (3,2):
        return (1-w**2)**0.5*(45*w**2-15)
    elif (l,m) == (3,3):
        return 45*(1-w**2)*w
    elif (l,m) == (4,0):
        return (1-w**2)**0.5*(-17.5*w**3+7.5*w)
    elif (l,m) == (4,1):
        return 70*w**4-67.5*w**2+7.5
    elif (l,m) == (4,2):
        return (1-w**2)**0.5*(210*w**3-120*w)
    elif (l,m) == (4,3):
        return (1-w**2)*(420*w**2-105)
    elif (l,m) == (4,4):
        return 420*(1-w**2)**1.5*w
    elif (l,m) == (5,0):
        return (1-w**2)**0.5*(-39.375*w**4+26.25*w**2-1.875)
    elif (l,m) == (5,1):
        return (196.875*w**5-236.25*w**3+54.375*w)
    elif (l,m) == (5,2):
        return (1-w**2)**0.5*(787.5*w**4-630*w**2+52.5)
    elif (l,m) == (5,3):
        return (1-w**2)*(2362.5*w**3-1102.5*w)
    elif (l,m) == (5,4):
        return (1-w**2)**1.5*(4725*w**2-945)
    elif (l,m) == (5,5):
        return 4725*(1-w**2)**2*w
    elif (l,m) == (6,0):
        return (1-w**2)**0.5*(-86.625*w**5+78.75*w**3-13.125*w)
    elif (l,m) == (6,1):
        return 519.75*w**6-748.125*w**4+262.5*w**2-13.125
    elif (l,m) == (6,2):
        return (1-w**2)**0.5*(2598.75*w**5-2677.5*w**3+498.75*w)
    elif (l,m) == (6,3):
        return (1-w**2)*(10395.0*w**4-7087.5*w**2+472.5)
    elif (l,m) == (6,4):
        return (1-w**2)**1.5*(31185*w**3-12285*w)
    elif (l,m) == (6,5):
        return (1-w**2)**2*(62370*w**2-10395)
    elif (l,m) == (6,6):
        return 62370*(1-w**2)**2.5*w
    elif (l,m) == (7,0):
        return (1-w**2)**0.5*(-187.6875*w**6+216.5625*w**4-59.0625*w**2+2.1875)
    elif (l,m) == (7,1):
        return 1313.8125*w**7-2208.9375*w**5+1043.4375*w**3-120.3125*w
    elif (l,m) == (7,2):
        return (1-w**2)**0.5*(7882.875*w**6-9961.875*w**4+2953.125*w**2-118.125)
    elif (l,m) == (7,3):
        return (1-w**2)*(39414.375*w**5-35516.25*w**3+5551.875*w)
    elif (l,m) == (7,4):
        return (1-w**2)**1.5*(157657.5*w**4-93555*w**2+5197.5)
    elif (l,m) == (7,5):
        return (1-w**2)**2*(472972.5*w**3-161122.5*w)
    elif (l,m) == (7,6):
        return (1-w**2)**2.5*(945945*w**2-135135)
    elif (l,m) == (7,7):
        return 945945*(1-w**2)**3*w
    elif (l,m) == (8,0):
        return (1-w**2)**0.5*(-402.1875*w**7+563.0625*w**5-216.5625*w**3+19.6875*w)
    elif (l,m) == (8,1):
        return 3217.5*w**8-6193.6875*w**6+3681.5625*w**4-689.0625*w**2+19.6875
    elif (l,m) == (8,2):
        return (1-w**2)**0.5*(22522.5*w**7-33783.75*w**5+13860.0*w**3-1338.75*w)
    elif (l,m) == (8,3):
        return (1-w**2)*(135135.0*w**6-152026.875*w**4+38981.25*w**2-1299.375)
    elif (l,m) == (8,4):
        return (1-w**2)**1.5*(675675.0*w**5-540540.0*w**3+72765.0*w)
    elif (l,m) == (8,5):
        return (1-w**2)**2*(2702700.0*w**4-1418917.5*w**2+67567.5)
    elif (l,m) == (8,6):
        return (1-w**2)**2.5*(8108100*w**3-2432430*w)
    elif (l,m) == (8,7):
        return (1-w**2)**3*(16216200*w**2-2027025)
    elif (l,m) == (8,8):
        return 16216200*(1-w**2)**3.5*w
    elif (l,m) == (9,0):
        return (1-w**2)**0.5*(-854.6484375*w**8+1407.65625*w**6-703.828125*w**4+108.28125*w**2-2.4609375)
    elif (l,m) == (9,1):
        return 7691.8359375*w**9-16690.78125*w**7+11965.078125*w**5-3140.15625*w**3+219.0234375*w
    elif (l,m) == (9,2):
        return (1-w**2)**0.5*(61534.6875*w**8-106981.875*w**6+56306.25*w**4-9095.625*w**2+216.5625)
    elif (l,m) == (9,3):
        return (1-w**2)*(430742.8125*w**7-582769.6875*w**5+211148.4375*w**3-17541.5625*w)
    elif (l,m) == (9,4):
        return (1-w**2)**1.5*(2584456.875*w**6-2618240.625*w**4+591215.625*w**2-16891.875)
    elif (l,m) == (9,5):
        return (1-w**2)**2*(12922284.375*w**5-9290531.25*w**3+1097971.875*w)
    elif (l,m) == (9,6):
        return (1-w**2)**2.5*(51689137.5*w**4-24324300*w**2+1013512.5)
    elif (l,m) == (9,7):
        return (1-w**2)**3*(155067412.5*w**3-41554012.5*w)
    elif (l,m) == (9,8):
        return (1-w**2)**3.5*(310134825*w**2-34459425)
    elif (l,m) == (9,9):
        return 310134825*(1-w**2)**4*w
    else:
        raise ValueError('Unsupported arguments (l,m)=(%d,%d)' % (l,m))

# -------------------------------------------------------------------

if __name__ == '__main__':
    from gpaw.mpi import world
    from gpaw.sphere import lmiter

    def equal(x, y, tol=0, msg=''):
        if abs(x-y) > tol:
            msg = (msg + '(%.9g,%.9g) != (%.9g,%.9g) (error: |%.9g| > %.9g)' %
                   (np.real(x),np.imag(x),np.real(y),np.imag(y),abs(x-y),tol))
            raise AssertionError(msg)

    nw = 10000
    dw = 2.0/nw
    w = np.linspace(-1+dw/2,1-dw/2,nw)
    equal(w[1]-w[0], dw, 1e-12)
    tol = 1e-5

    # Test orthogonality of the associated Legendre polynomials
    if world.rank == 0:
        print('\n%s\nAssociated Legendre orthogonality\n%s' % ('-'*40,'-'*40))
    for l1,m1 in lmiter(lmax, False, comm=world):
        p1 = legendre(l1, m1, w)
        for l2,m2 in lmiter(lmax, False):
            p2 = legendre(l2, m2, w)
            scale = (ilegendre(l1,m1)*ilegendre(l2,m2))**0.5
            v = np.dot(p1,p2)*dw / scale
            if (l1,m1) == (l2,m2):
                v0 = ilegendre(l1,m1) / scale
                print('l1=%2d, m1=%2d, l2=%2d, m2=%2d, v=%12.9f, err=%12.9f' % (l1,m1,l2,m2,v,abs(v-v0)))
                equal(v, v0, tol, 'l1=%2d, m1=%2d, l2=%2d, m=%2d: ' % (l1,m1,l2,m2))
            elif m1 == m2:
                print('l1=%2d, m1=%2d, l2=%2d, m2=%2d, v=%12.9f, err=%12.9f' % (l1,m1,l2,m2,v,abs(v)))
                equal(v, 0, tol, 'l1=%2d, m1=%2d, l2=%2d, m=%2d: ' % (l1,m1,l2,m2))
    del nw, dw, w, tol, p1, p2, v, v0
    world.barrier()

    # =======================

    nw = 50000
    dw = 1.9/nw
    w = np.linspace(-0.95+dw/2,0.95-dw/2,nw)
    tol = 1e-5

    # Test theta-derivative of the associated Legendre polynomials
    if world.rank == 0:
        print('\n%s\nAssociated Legendre theta-derivative\n%s' % ('-'*40,'-'*40))
    for l,m in lmiter(lmax, False, comm=world):
        scale = l**(-abs(m))*lmfact(l,abs(m)) # scale to fit ~ [-1; 1]
        p = legendre(l, m, w) / scale
        dpdtheta = -(1-w[1:-1]**2)**0.5*(p[2:]-p[:-2])/(2.0*dw)
        dpdtheta0 = dlegendre(l, m, w[1:-1]) / scale
        e = np.sum((dpdtheta-dpdtheta0)**2)**0.5/(nw-2.0)**0.5
        print('l=%2d, m=%2d, err=%12.9f, max=%12.9f' % (l,m,e,np.abs(dpdtheta-dpdtheta0).max()))
        #if e > tol:
        #    import pylab as pl
        #    w, p = w[1:-1], p[1:-1]
        #    fig = pl.figure(lm_to_L(l,m))
        #    ax = pl.axes()
        #    ax.plot(w, p,'-r', w, dpdtheta, '-b', w, dpdtheta0, '-g')
        #    fig = pl.figure((lmax+1)**2+lm_to_L(l,m))
        #    ax = pl.axes()
        #    ax.plot(w, dpdtheta-dpdtheta0, '-k')
        #    pl.show()
        equal(e, 0, tol, 'l=%2d, m=%2d: ' % (l,m))
    del nw, dw, w, tol, p, dpdtheta, dpdtheta0, e
    world.barrier()

    # =======================

    R = 1.0
    npts = 1000
    tol = 1e-9
    # ((R-dR)/(R+dR))**(lmax+1) = tol
    # (lmax+1)*np.log((R-dR)/(R+dR)) = np.log(tol)
    # (R-dR)/(R+dR) = np.exp(np.log(tol)/(lmax+1))
    # R-dR = (R+dR) * tol**(1/(lmax+1))
    # R * (1-tol**(1/(lmax+1))) = dR * (1+tol**(1/(lmax+1)))
    dR = R * (1-tol**(1./(lmax+1))) / (1+tol**(1./(lmax+1)))
    assert abs(((R-dR)/(R+dR))**(lmax+1) - tol) < 1e-12

    r_g = np.random.uniform(R+dR, 10*R, size=npts)
    world.broadcast(r_g, 0)
    theta_g = np.random.uniform(0, np.pi, size=npts)
    world.broadcast(theta_g, 0)
    phi_g = np.random.uniform(0, np.pi, size=npts)
    world.broadcast(phi_g, 0)

    r_vg = np.empty((3, npts), dtype=float)
    r_vg[0] = r_g*np.cos(phi_g)*np.sin(theta_g)
    r_vg[1] = r_g*np.sin(phi_g)*np.sin(theta_g)
    r_vg[2] = r_g*np.cos(theta_g)

    rm_g = np.random.uniform(0, R-dR, size=npts)
    world.broadcast(rm_g, 0)
    thetam_g = np.random.uniform(0, np.pi, size=npts)
    world.broadcast(thetam_g, 0)
    phim_g = np.random.uniform(0, np.pi, size=npts)
    world.broadcast(phim_g, 0)

    rm_vg = np.empty((3, npts), dtype=float)
    rm_vg[0] = rm_g*np.cos(phim_g)*np.sin(thetam_g)
    rm_vg[1] = rm_g*np.sin(phim_g)*np.sin(thetam_g)
    rm_vg[2] = rm_g*np.cos(thetam_g)

    # Cosine of the angle between r and r' using 1st spherical law of cosines
    w_g = np.cos(theta_g)*np.cos(thetam_g) \
        + np.sin(theta_g)*np.sin(thetam_g)*np.cos(phi_g-phim_g)

    f0_g = np.sum((r_vg-rm_vg)**2, axis=0)**(-0.5)
    f_g = np.zeros_like(f0_g)

    # Test truncated multi-pole expansion using Legendre polynomials
    if world.rank == 0:
        print('\n%s\nLegendre multi-pole expansion\n%s' % ('-'*40,'-'*40))
    for l in range(lmax+1):
        f_g += 1/r_g * (rm_g/r_g)**l * legendre(l, 0, w_g)

    e = np.abs(f_g-f0_g).max()
    equal(e, 0, tol, 'lmax=%02d, dR=%g: ' % (lmax, dR))
    del R, dR, npts, tol, r_g, theta_g, phi_g, r_vg, rm_g, thetam_g, phim_g, rm_vg, w_g, f0_g, f_g, e

