from math import acos, cos, sin, sqrt
import numpy as np
from ase.atoms import string2vector

class Vector3d(list):
    def __init__(self,vector=None):
        if vector is None or vector == []:
            vector = [0,0,0]
        vector = string2vector(vector)
        list.__init__(self)
        for c in range(3):
            self.append(float(vector[c]))
        self.l = False

    def __add__(self, other):
        result = self.copy()
        for c in range(3):
            result[c] += other[c]
        return result

    def __div__(self,other):
        return Vector3d(np.array(self) / other)

    def __mul__(self, x):
        if type(x) == type(self):
            return np.dot( self, x )
        else:
            return Vector3d(x * np.array(self))
        
    def __rmul__(self, x):
        return self.__mul__(x)
        
    def __lmul__(self, x):
        return self.__mul__(x)

    def __neg__(self):
        return -1 * self
        
    def __str__(self):
        return "(%g,%g,%g)" % tuple(self)

    def __sub__(self, other):
        result = self.copy()
        for c in range(3):
            result[c] -= other[c]
        return result

    def angle(self, other):
        """Return the angle between the directions of yourself and the
        other vector in radians."""
        other = Vector3d(other)
        ll = self.length() * other.length()
        if not ll > 0:
            return None
        return acos((self * other) / ll)
        
    def copy(self):
        return Vector3d(self)

    def distance(self,vector):
        if type(vector) is not type(self):
            vector=Vector3d(vector)
        dv = self - vector
        return (self - vector).length()

    def length(self,value=None):
        if value:
            fac = value / self.length()
            for c in range(3):
                self[c] *= fac
            self.l = False
        if not self.l:
            self.l = sqrt(self.norm())
        return self.l

    def norm(self):
        #return np.sum( self*self )
        return self*self  #  XXX drop this class and use numpy arrays ...

    def rotation_axis(self, other):
        """Return the rotation axis to rotate yourself to the direction
        of the other vector. The length of the axis gives the rotation angle
        (in radians)."""
        other = Vector3d(other)
        angle = self.angle(other)
        if angle is None:
            return None
        axis = self.vprod(other)
        axis.length(angle)
        return axis

    def rotate(self, axis, angle=None):
        """Rotate the vector about the given axis with the given angle.
        Note, that the right hand rule applies: If your right thumb points
        into the direction of the axis, the other fingers show the rotation
        direction."""
        axis=Vector3d(axis)
        if angle is None:
            angle=axis.length()
        axis.length(1.) 
        res = (cos(angle) * self - sin(angle) * self.vprod(axis) +
               ((self * axis) * (1. - cos(angle))) * axis          )
        for c in range(3):
            self[c] = res[c]
            
    def vprod(self, a, b=None):
        """vector product"""
        if b is None:
            # [self x a]
            return Vector3d([self[1]*a[2]-self[2]*a[1],
                             self[2]*a[0]-self[0]*a[2],
                             self[0]*a[1]-self[1]*a[0]])
        else:
            # [a x b]
            return Vector3d([a[1]*b[2]-a[2]*b[1],
                             a[2]*b[0]-a[0]*b[2],
                             a[0]*b[1]-a[1]*b[0]])
                         
    def x(self):
        return self[0]

    def y(self):
        return self[1]

    def z(self):
        return self[2]
