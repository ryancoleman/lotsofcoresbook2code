# -*- coding: utf-8 -*-

try:
    import numpy as np
except ImportError:
    # required due to ase/test/eoswoase.py
    pass

class EquationOfStateSJEOS:
    """Fit equation of state for bulk systems.

    10.1103/PhysRevB.67.026103

    The following equation is used::

       A third order inverse polynomial fit

                           2      3        -1/3
       E(V) = c + c t + c t  + c t ,  t = V
               0   1     2      3

    More methods available in ase.utils.eosase2

    Use::

       eos = EquationOfState(volumes, energies)
       v0, e0, B = eos.fit()
       eos.plot()

    """
    def __init__(self, volumes, energies, eos='sjeos'):
        assert eos == 'sjeos', eos + ' eos not available. Probably scipy missing.'
        self.v = np.array(volumes)
        self.e = np.array(energies)
        self.eos_string = 'sjeos'

        self.v0 = None

    def fit(self):
        """Calculate volume, energy, and bulk modulus.

        Returns the optimal volume, the minumum energy, and the bulk
        modulus.  Notice that the ASE units for the bulk modulus is
        eV/Angstrom^3 - to get the value in GPa, do this::

          v0, e0, B = eos.fit()
          print B / kJ * 1.0e24, 'GPa'

        """

        fit0 = np.poly1d(np.polyfit(self.v**-(1.0 / 3), self.e, 3))
        fit1 = np.polyder(fit0, 1)
        fit2 = np.polyder(fit1, 1)

        self.v0 = None
        for t in np.roots(fit1):
            if isinstance(t, float) and t > 0 and fit2(t) > 0:
                self.v0 = t**-3
                break

        if self.v0 is None:
            raise ValueError('No minimum!')

        self.e0 = fit0(t)
        self.B = t**5 * fit2(t) / 9
        self.fit0 = fit0

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
        y = self.fit0(x**-(1.0 / 3))
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

if __name__ == '__main__':
    try:
        import numpy as np
        # from ase/test/eos.py
        volumes = [29.205536, 30.581492, 32.000000, 33.461708, 34.967264]
        energies = [0.0190898, -0.0031172, -0.0096925, -0.0004014, 0.0235753]
        sjeos = (31.867118229937798, -0.0096410046694188622, 0.23984474782755572)
        eos = EquationOfStateSJEOS(volumes, energies)
        v0, e0, B = eos.fit()
        assert abs(v0 - sjeos[0]) < 5.e-6
        assert abs(B - sjeos[2]) < 5.e-6
    except ImportError:
        pass
