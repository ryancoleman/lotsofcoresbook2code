# -*- coding: utf-8 -*-

import numpy as np

from ase.test import NotAvailable

try:
    import scipy
except ImportError:
    raise NotAvailable('This needs scipy module.')

try:
    from scipy.optimize import curve_fit
except ImportError:
    from scipy.optimize import leastsq

    # this part comes from
    # http://projects.scipy.org/scipy/browser/trunk/scipy/optimize/minpack.py
    def _general_function(params, xdata, ydata, function):
        return function(xdata, *params) - ydata
    # end of this part

    def curve_fit(f, x, y, p0):
        func = _general_function
        args = (x, y, f)
        # this part comes from
        # http://projects.scipy.org/scipy/browser/trunk/scipy/optimize/minpack.py
        popt, pcov, infodict, mesg, ier = leastsq(func, p0, args=args, full_output=1)

        if ier not in [1,2,3,4]:
            raise RuntimeError, "Optimal parameters not found: " + mesg
        # end of this part
        return popt, pcov

def taylor(V, E0, beta, alpha, V0):
    'Taylor Expansion up to 3rd order about V0'

    E = E0 + beta/2.*(V-V0)**2/V0 + alpha/6.*(V-V0)**3/V0
    return E

def murnaghan(V, E0, B0, BP, V0):
    'From PRB 28,5480 (1983'

    E = E0 + B0*V/BP*(((V0/V)**BP)/(BP-1)+1) - V0*B0/(BP-1)
    return E

def birch(V, E0, B0, BP, V0):
    '''
    From Intermetallic compounds: Principles and Practice, Vol. I: Principles
    Chapter 9 pages 195-210 by M. Mehl. B. Klein, D. Papaconstantopoulos
    paper downloaded from Web

    case where n=0
    '''

    E = (E0
         + 9.0/8.0*B0*V0*((V0/V)**(2.0/3.0) - 1.0)**2
         + 9.0/16.0*B0*V0*(BP-4.)*((V0/V)**(2.0/3.0) - 1.0)**3)
    return E

def birchmurnaghan(V, E0, B0, BP, V0):
    'BirchMurnaghan equation from PRB 70, 224107'

    eta = (V/V0)**(1./3.)
    E = E0 + 9.*B0*V0/16.*(eta**2-1)**2*(6 + BP*(eta**2-1.) - 4.*eta**2)
    return E

def pouriertarantola(V, E0, B0, BP, V0):
    'Pourier-Tarantola equation from PRB 70, 224107'

    eta = (V/V0)**(1./3.)
    squiggle = -3.*np.log(eta)

    E = E0 + B0*V0*squiggle**2/6.*(3. + squiggle*(BP - 2))
    return E

def vinet(V, E0, B0, BP, V0):
    'Vinet equation from PRB 70, 224107'

    eta = (V/V0)**(1./3.)

    E = (E0 + 2.*B0*V0/(BP-1.)**2
         * (2. - (5. +3.*BP*(eta-1.)-3.*eta)*np.exp(-3.*(BP-1.)*(eta-1.)/2.)))
    return E

def antonschmidt(V, Einf, B, n, V0):
    '''From Intermetallics 11, 23-32 (2003)

    Einf should be E_infinity, i.e. infinite separation, but
    according to the paper it does not provide a good estimate
    of the cohesive energy. They derive this equation from an
    empirical formula for the volume dependence of pressure,

    E(vol) = E_inf + int(P dV) from V=vol to V=infinity

    but the equation breaks down at large volumes, so E_inf
    is not that meaningful

    n should be about -2 according to the paper.

    I find this equation does not fit volumetric data as well
    as the other equtions do.
    '''

    E = B*V0/(n+1.) * (V/V0)**(n+1.)*(np.log(V/V0)-(1./(n+1.))) + Einf
    return E

def p3(V, c0, c1, c2, c3):
    'polynomial fit'

    E = c0 + c1*V + c2*V**2 + c3*V**3
    return E

def parabola(x, a, b, c):
    '''
    parabola polynomial function

    this function is used to fit the data to get good guesses for
    the equation of state fits

    a 4th order polynomial fit to get good guesses for
    was not a good idea because for noisy data the fit is too wiggly
    2nd order seems to be sufficient, and guarentees a single minimum'''

    return a + b*x + c*x**2

class EquationOfStateASE2:
    """Fit equation of state for bulk systems.

    The following equation is used::

       taylor
           A third order Taylor series expansion about the minimum volume

       murnaghan
           PRB 28, 5480 (1983)

       birch
           Intermetallic compounds: Principles and Practice,
           Vol I: Principles. pages 195-210

       birchmurnaghan
           PRB 70, 224107

       pouriertarantola
           PRB 70, 224107

       vinet
           PRB 70, 224107

       antonschmidt
           Intermetallics 11, 23-32 (2003)

       p3
           A third order polynomial fit

    Use::

       eos = EquationOfState(volumes, energies, eos='murnaghan')
       v0, e0, B = eos.fit()
       eos.plot()

    """
    def __init__(self, volumes, energies, eos='murnaghan'):
        self.v = np.array(volumes)
        self.e = np.array(energies)
        self.eos_string = eos

        self.v0 = None

    def fit(self):
        """Calculate volume, energy, and bulk modulus.

        Returns the optimal volume, the minumum energy, and the bulk
        modulus.  Notice that the ASE units for the bulk modulus is
        eV/Angstrom^3 - to get the value in GPa, do this::

          v0, e0, B = eos.fit()
          print B / kJ * 1.0e24, 'GPa'

        """

        # old ASE2 implementation

        p0 = [min(self.e), 1, 1]
        popt, pcov = curve_fit(parabola, self.v, self.e, p0)

        parabola_parameters = popt
        ## Here I just make sure the minimum is bracketed by the volumes
        ## this if for the solver
        minvol = min(self.v)
        maxvol = max(self.v)

        # the minimum of the parabola is at dE/dV = 0, or 2*c V +b =0
        c = parabola_parameters[2]
        b = parabola_parameters[1]
        a = parabola_parameters[0]
        parabola_vmin = -b/2/c

        if not (minvol < parabola_vmin and parabola_vmin < maxvol):
            print 'Warning the minimum volume of a fitted parabola is not in your volumes. You may not have a minimum in your dataset'

        # evaluate the parabola at the minimum to estimate the groundstate energy
        E0 = parabola(parabola_vmin, a, b, c)
        # estimate the bulk modulus from Vo*E''.  E'' = 2*c
        B0 = 2*c*parabola_vmin

        if self.eos_string == 'antonschmidt':
            BP = -2
        else:
            BP = 4

        initial_guess = [E0, B0, BP, parabola_vmin]

        # now fit the equation of state
        p0 = initial_guess
        popt, pcov = curve_fit(eval(self.eos_string), self.v, self.e, p0)

        self.eos_parameters = popt

        if self.eos_string == 'p3':
            c0, c1, c2, c3 = self.eos_parameters
            # find minimum E in E = c0 + c1*V + c2*V**2 + c3*V**3
            # dE/dV = c1+ 2*c2*V + 3*c3*V**2 = 0
            # solve by quadratic formula with the positive root

            a = 3 * c3
            b = 2 * c2
            c = c1

            self.v0 = (-b + np.sqrt(b**2 - 4*a*c))/(2*a)
            self.e0 = p3(self.v0, c0, c1, c2, c3)
            self.B = (2*c2 + 6*c3*self.v0)*self.v0
        else:
            self.v0 = self.eos_parameters[3]
            self.e0 = self.eos_parameters[0]
            self.B = self.eos_parameters[1]

        return self.v0, self.e0, self.B

    def plot(self, filename=None, show=None):
        """Plot fitted energy curve.

        Uses Matplotlib to plot the energy curve.  Use *show=True* to
        show the figure and *filename='abc.png'* or
        *filename='abc.eps'* to save the figure to a file."""

        #import matplotlib.pyplot as plt
        import pylab as plt

        if self.v0 is None:
            self.fit()

        if filename is None and show is None:
            show = True

        x = 3.95
        f = plt.figure(figsize=(x * 2.5**0.5, x))
        f.subplots_adjust(left=0.12, right=0.9, top=0.9, bottom=0.15)
        plt.plot(self.v, self.e, 'o')
        x = np.linspace(min(self.v), max(self.v), 100)
        y = eval(self.eos_string)(x,
                                  self.eos_parameters[0],
                                  self.eos_parameters[1],
                                  self.eos_parameters[2],
                                  self.eos_parameters[3],
                                  )
        plt.plot(x, y, '-r')
        try:
            from ase.units import kJ
            plt.xlabel(u'volume [Å^3]')
            plt.ylabel(u'energy [eV]')
            plt.title(u'%s: E: %.3f eV, V: %.3f Å^3, B: %.3f GPa' %
                      (self.eos_string, self.e0, self.v0, self.B / kJ * 1.e24))
        except ImportError:
            plt.xlabel(u'volume [L(length)^3]')
            plt.ylabel(u'energy [E(energy)]')
            plt.title(u'%s: E: %.3f E, V: %.3f L^3, B: %.3e E/L^3' %
                      (self.eos_string, self.e0, self.v0, self.B))

        if show:
            plt.show()
        if filename is not None:
            f.savefig(filename)

        return f
