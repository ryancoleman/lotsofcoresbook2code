"""Different kinds of laser pulses.
"""

import sys
import time
from math import log

import numpy as np
from ase.units import Bohr, Hartree

import gpaw.io
from gpaw.tddft.units import attosec_to_autime, autime_to_attosec, \
                             eV_to_aufrequency, aufrequency_to_eV

###########################
# Main class
###########################
class LaserField:
  """
    
  """    
  def __init__(self):
      pass
  
  def strength(self, t):
      return np.array([0.0, 0.0, 0.0])



class CWField(LaserField):
  """
  Continuously oscillating laser field which is switch on linearly.

  Parameters:
    e0   field strength  (in atomic units)
    w    field frequency (in atomic units)
    ts   switch on time  (in atomic units)    
  """
  def __init__(self, e0, w, ts):
    # FIXME: use eV, ang and attosec
    #           attosec_to_autime, autime_to_attosec, \
    #           eV_to_aufrequency, aufrequency_to_eV
    self.e0 = e0
    self.w  = w
    self.ts = ts
    
  def strength(self, t):
    if t < self.ts:
      c = self.e0 * (t / self.ts) * np.sin(self.w * t)
    else:
      c = self.e0 * np.sin(self.w * t)
    return np.array([0.0, 0.0, c])
