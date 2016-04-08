def factorial(x):
    """Return x!, where x is a non-negative integer"""
    if x < 2: return 1
    else: return x * factorial(x - 1)

def gcd(a, b):
    """Return greatest common divisor of a and b.
       Uses the euclidian algorithm
    """
    if b == 0: return a
    else: return gcd(b, a % b)

class Rational:
    """Class used to represent rational numbers as fractions, such that
       no precision is lost during calculation operations.
       Example usage with Numeric:
       import numpy as np
       from tools import Rational as Q

       n = np.zeros(4, 'O')                   array([0 , 0 , 0 , 0 ],'O')
       n[2:4] = [Q(35,12), Q(36,12)]           array([0 , 0 , 35./12 , 3 ],'O')
       24 * n                                  array([0 , 0 , 70 , 72 ],'O')
       np.multiply(n, Q(3,9))                 array([0 , 0 , 35./36 , 1 ],'O')
       
    """
    def __init__(self, nom=0, denom=1):
##        assert type(nom) == type(denom) == int

        # ensure that sign is in the nominator
        nom = cmp(denom, 0) * nom
        denom = abs(denom)

        # reduce fraction
        q = gcd(nom, denom)
        self.nom = nom / q
        self.denom = denom / q

    def __add__(self, x):
        if type(x) == float:
            return float(self) + x
        elif type(x) == int:
            x = Rational(x)
        nom = self.nom * x.denom + x.nom * self.denom
        denom = self.denom * x.denom
        return Rational(nom, denom)

    def __radd__(self, x):
        return self.__add__(x)

    def __mul__(self, x):
        if type(x) == float:
            return float(self) * x
        elif type(x) == int:
            x = Rational(x)
        return Rational(self.nom * x.nom, self.denom * x.denom)

    def __rmul__(self, x):
        return self.__mul__(x)

    def __neg__(self):
        return Rational(-self.nom, self.denom)

    def __pos__(self):
        return self.copy()

    def __sub__(self, x):
        return self.__add__(-x)

    def __rsub__(self, x):
        return -self.__sub__(x)
    
    def __div__(self, x):
        if type(x) == float:
            return float(self) / x
        elif type(x) == int:
            x = Rational(x)
        return self.__mul__(Rational(x.denom, x.nom))

    def __rdiv__(self, x):
        if type(x) == float:
            return x / float(self)
        elif type(x) == int:
            x = Rational(x)
        return x.__mul__(Rational(self.denom, self.nom))

    def __pow__(self, p):
        if p == 0: return Rational(1)
        if p >= 0 and type(p) == int:
            return Rational(self.nom**p, self.denom**p)
        else:
            return float(self)**p

    def __mod__(self, x):
        if type(x) == float:
            return float(self) % x
        return Rational(self.nom % (x * self.denom), self.denom)

    def __rmod__(self, x):
        if type(x) == int:
            x = Rational(x)
        i = self.__int__()
        return x.__mod__(i)

    def __abs__(self):
        return Rational(abs(self.nom), self.denom)
    
    def __nonzero__(self):
        return self.nom.__nonzero__()

    def __cmp__(self, x):
        return cmp(float(self), float(x))    
        
    def __str__(self):
        out = str(self.nom)
        if self.denom != 1:
            out += './' + str(self.denom)
        return out

    def __int__(self):
        assert self.denom == 1
        return self.nom

    def __float__(self):
        return float(self.nom) / self.denom

    def __repr__(self):
        out = repr(self.nom)
        if self.denom != 1:
            out += './' + repr(self.denom)
        return out

    def __copy__(self):
        return Rational(self.nom, self.denom)

    def floor(self):
        return int(float(self))

    def sqrt(self):
        return self**.5

    def abs(self):
        return Rational(abs(self.nom), self.denom)
    
    def copy(self):
        return Rational(self.nom, self.denom)
