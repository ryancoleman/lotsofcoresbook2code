import numpy as np
from scipy.special import dawsn

from gpaw.gauss import Gauss, Lorentz


class ComplexLorentz:
    def __init__(self, width=0.08):
        self.dtype = complex
        self.set_width(width)
        
    def get(self, x, x0):
        return 0.5 / x0 * (((x + x0) / ((x + x0)**2 + self.width2)
                            - (x - x0) / ((x - x0)**2 + self.width2))
                           - 1.0j * self.width *
                           (1 / ((x + x0)**2 + self.width2)
                            - 1 / ((x - x0)**2 + self.width2))
                           )
        
    def set_width(self, width=0.08):
        self.width = width
        self.width2 = width**2


class ComplexGauss:
    def __init__(self, width=0.08):
        self.dtype = complex
        self.set_width(width)
        
    def get(self, x, x0):
        return 0.5 / x0 * (2 * self.wm1 *
                           (dawsn((x + x0) * self.wm1)
                            - dawsn((x - x0) * self.wm1))
                            - 1.0j * self.norm *
                            (np.exp(-((x + x0) * self.wm1)**2)
                             - np.exp(-((x - x0) * self.wm1)**2))
                           )
    
    def set_width(self, width=0.08):
        self.norm = 1. / width * np.sqrt(np.pi / 2)
        self.wm1 = np.sqrt(.5) / width


class Folder:
    """Fold a function with normalised Gaussians or Lorentzians"""
    def __init__(self, width,
                 folding='Gauss'):
        self.width = width
        if folding == 'Gauss':
            self.func = Gauss(width)
        elif folding == 'Lorentz':
            self.func = Lorentz(width)
        elif folding == 'ComplexGauss':
            self.func = ComplexGauss(width)
        elif folding == 'ComplexLorentz':
            self.func = ComplexLorentz(width)
        elif folding is None:
            self.func = None
        else:
            raise RuntimeError('unknown folding "' + folding + '"')

    def fold(self, x, y, dx=None, xmin=None, xmax=None):
        X = np.array(x)
        assert len(X.shape) == 1
        Y = np.array(y)
        assert X.shape[0] == Y.shape[0]

        if self.func is None:
            xl = X
            yl = Y
        else:
            if xmin is None:
                xmin = np.min(X) - 4 * self.width
            if xmax is None:
                xmax = np.max(X) + 4 * self.width
            if dx is None:
                dx = self.width / 4.

            xl = np.arange(xmin, xmax + 0.5 * dx, dx)
            
            xl, yl = self.fold_values(x, y, xl)
            
        return xl, yl

    def fold_values(self, x, y, xl=None):
        X = np.array(x)
        assert len(X.shape) == 1
        Y = np.array(y)
        assert X.shape[0] == Y.shape[0]

        if self.func is None:
            xl = X
            yl = Y
        else:
            if xl is None:
                Xl = np.unique(X)
            else:
                Xl = np.array(xl)
                assert len(Xl.shape) == 1
            
            # weight matrix
            weightm = np.empty((Xl.shape[0], X.shape[0]),
                               dtype=self.func.dtype)
            for i, x in enumerate(X):
                weightm[:, i] = self.func.get(Xl, x)

#            yl = np.dot(weightm, Y)
            yl = np.tensordot(weightm, Y, axes=(1, 0))

        return Xl, yl
