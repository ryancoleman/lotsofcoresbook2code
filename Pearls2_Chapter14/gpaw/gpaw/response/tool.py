import numpy as np
from scipy.optimize import leastsq
import pylab as pl

def check_degenerate_bands(filename, etol):

    from gpaw import GPAW
    calc = GPAW(filename,txt=None)
    print('Number of Electrons   :', calc.get_number_of_electrons())
    nibzkpt = calc.get_ibz_k_points().shape[0]
    nbands = calc.get_number_of_bands()
    print('Number of Bands       :', nbands)
    print('Number of ibz-kpoints :', nibzkpt)
    e_kn = np.array([calc.get_eigenvalues(k) for k in range(nibzkpt)])
    f_kn = np.array([calc.get_occupation_numbers(k) for k in range(nibzkpt)])
    for k in range(nibzkpt):
        for n in range(1,nbands):
            if (f_kn[k,n-1] - f_kn[k,n] > 1e-5) and (np.abs(e_kn[k,n] - e_kn[k,n-1]) < etol):
                print(k, n, e_kn[k,n], e_kn[k, n-1])
    return
    

def get_orbitals(calc):
    """Get LCAO orbitals on 3D grid by lcao_to_grid method."""

    bfs_a = [setup.phit_j for setup in calc.wfs.setups]
    
    from gpaw.lfc import BasisFunctions
    bfs = BasisFunctions(calc.wfs.gd, bfs_a, calc.wfs.kd.comm, cut=True)
    spos_ac = calc.atoms.get_scaled_positions()
    bfs.set_positions(spos_ac)

    nLCAO = calc.get_number_of_bands()
    orb_MG = calc.wfs.gd.zeros(nLCAO)
    C_M = np.identity(nLCAO)
    bfs.lcao_to_grid(C_M, orb_MG,q=-1)
    
    return orb_MG

def find_peaks(x,y,threshold = None):
    """ Find peaks for a certain curve.

    Usage:
    threshold = (xmin, xmax, ymin, ymax)

    """

    assert type(x) == np.ndarray and type(y) == np.ndarray
    assert x.ndim == 1 and y.ndim == 1
    assert x.shape[0] == y.shape[0]

    if threshold is None:
        threshold = (x.min(), x.max(), y.min(), y.max())
        
    if type(threshold) is not tuple:
        threshold = (threshold, )

    if len(threshold) == 1:
        threshold += (x.max(), y.min(), y.max())
    elif len(threshold) == 2:
        threshold += (y.min(), y.max())
    elif len(threshold) == 3:
        threshold += (y.max(),)
    else:
        pass

    xmin = threshold[0]
    xmax = threshold[1]
    ymin = threshold[2]
    ymax = threshold[3]

    peak = {}
    npeak = 0
    for i in range(1, x.shape[0]-1):
        if (y[i] >= ymin and y[i] <= ymax and
            x[i] >= xmin and x[i] <= xmax ):
            if y[i] > y[i-1] and y[i] > y[i+1]:
                peak[npeak] = np.array([x[i], y[i]])
                npeak += 1

    peakarray = np.zeros([npeak, 2])
    for i in range(npeak):
        peakarray[i] = peak[i]
        
    return peakarray

def lorz_fit(x,y, npeak=1, initpara = None):
    """ Fit curve using Lorentzian function

    Note: currently only valid for one and two lorentizian

    The lorentzian function is defined as::

                      A w
        lorz = --------------------- + y0
                (x-x0)**2 + w**2

    where A is the peak amplitude, w is the width, (x0,y0) the peak position

    Parameters:

    x, y: ndarray
        Input data for analyze
    p: ndarray
        Parameters for curving fitting function. [A, x0, y0, w]
    p0: ndarray
        Parameters for initial guessing. similar to p
    
    """

    
    def residual(p, x, y):
        
        err = y - lorz(x, p, npeak)
        return err

    def lorz(x, p, npeak):

        if npeak == 1:        
            return p[0] * p[3] / ( (x-p[1])**2 + p[3]**2 ) + p[2]
        if npeak == 2:
            return ( p[0] * p[3] / ( (x-p[1])**2 + p[3]**2 ) + p[2]
                   + p[4] * p[7] / ( (x-p[5])**2 + p[7]**2 ) + p[6] )
        else:
            raise ValueError('Larger than 2 peaks not supported yet!')

    if initpara is None:
        if npeak == 1:
            initpara[i] = np.array([1., 0., 0., 0.1])
        if npeak == 2:
            initpara[i] = np.array([1., 0., 0., 0.1,
                                    3., 0., 0., 0.1])
    p0 = initpara
    
    result = leastsq(residual, p0, args=(x, y), maxfev=2000)

    yfit = lorz(x, result[0],npeak)

    return yfit, result[0]

def linear_fit(x,y, initpara = None):
    
    def residual(p, x, y):
        err = y - linear(x, p)
        return err

    def linear(x, p):
        return p[0]*x + p[1]

    if initpara is None:
        initpara = np.array([1.0,1.0])

    p0 = initpara
    result = leastsq(residual, p0, args=(x, y), maxfev=2000)
    yfit = linear(x, result[0])

    return yfit, result[0]


def plot_setfont():

    params = {
              'axes.labelsize': 18,
              'text.fontsize': 18,
              'legend.fontsize': 18,
              'xtick.labelsize': 18,
              'ytick.labelsize': 18,
              'text.usetex': True}
    #          'figure.figsize': fig_size}
    pl.rcParams.update(params)


def plot_setticks(x=True, y=True):

    pl.minorticks_on()
    ax = pl.gca()
    if x:
        ax.xaxis.set_major_locator(pl.AutoLocator())
        x_major = ax.xaxis.get_majorticklocs()
        dx_minor =  (x_major[-1]-x_major[0])/(len(x_major)-1) /5.
        ax.xaxis.set_minor_locator(pl.MultipleLocator(dx_minor))
    else:
        pl.minorticks_off()
        
    if y:
        ax.yaxis.set_major_locator(pl.AutoLocator())
        y_major = ax.yaxis.get_majorticklocs()
        dy_minor =  (y_major[-1]-y_major[0])/(len(y_major)-1) /5.
        ax.yaxis.set_minor_locator(pl.MultipleLocator(dy_minor)) 
    else:
        pl.minorticks_off()
