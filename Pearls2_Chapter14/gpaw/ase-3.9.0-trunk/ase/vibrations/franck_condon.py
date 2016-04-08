from operator import mul
from itertools import combinations, product, chain
from math import factorial

import numpy as np
from array import *

from ase.units import kg, C, _hbar, kB
from ase.vibrations import Vibrations

class FranckCondon:
    def __init__(self, atoms, vibname, minfreq=None, maxfreq=None):
        """Input is a atoms object and the corresponding vibrations. 
        With minfreq and maxfreq frequencies can    
        be excluded from the calculation"""
     
        self.atoms = atoms
        # V = a * v is the combined atom and xyz-index
        self.mm05_V = np.repeat(1. / np.sqrt(atoms.get_masses()), 3)
        self.minfreq = minfreq
        self.maxfreq = maxfreq
        self.shape = (len(self.atoms), 3) 

        vib = Vibrations(atoms, name=vibname)
        self.energies = np.real(vib.get_energies(method='frederiksen')) # [eV] 
        self.frequencies = np.real(
            vib.get_frequencies(method='frederiksen')) # [cm^-1]
        self.modes = vib.modes
        self.H = vib.H
  
    def get_Huang_Rhys_factors(self, forces):
        """Evaluate Huang-Rhys factors and corresponding frequencies 
        from forces on atoms in the exited electronic state.        
        The double harmonic approximation is used. HR factors are 
        the first approximation of FC factors, 
        no combinations or higher quanta (>1) exitations are considered"""

        assert(forces.shape == self.shape) 

        # Hesse matrix
        H_VV = self.H
        # sqrt of inverse mass matrix
        mm05_V = self.mm05_V
        # mass weighted Hesse matrix 
        Hm_VV = mm05_V[:, None] * H_VV * mm05_V
        # mass weighted displacements
        Fm_V = forces.flat * mm05_V
        X_V = np.linalg.solve(Hm_VV, Fm_V)
        # projection onto the modes
        modes_VV = self.modes
        #d_V=np.dot([mode for mode in modes_VV] ,X_V)
        d_V = np.dot(modes_VV, X_V)
        # Huang-Rhys factors S
        s = 1.e-20 / kg / C / _hbar**2 # SI units
        S = s * d_V**2 * self.energies / 2

        # reshape for minfreq
        indices = np.where(self.frequencies <= self.minfreq)
        np.append(indices, np.where(self.frequencies >= self.maxfreq)) 
        S=np.delete(S, indices)
        frequencies = np.delete(self.frequencies, indices)
        
        return S, frequencies

    def get_Franck_Condon_factors(self, order, temp, forces):
        """Return FC factors and corresponding frequencies up to given order.
        
        order= number of quanta taken into account
        T= temperature in K. Vibronic levels are occupied by a 
        Boltzman distribution.
        forces= forces on atoms in the exited electronic state"""
        
        S, f = self.get_Huang_Rhys_factors(forces)
        n = order + 1
        T = temp
        freq = np.array(f)
        
        # frequencies
        freq_n = [[]*i for i in range(n-1)]
        freq_neg = [[]*i for i in range(n-1)]

        for i in range(1,n):
            freq_n[i-1]=freq*i
            freq_neg[i-1]=freq*(-i)
            
        # combinations
        freq_nn=[x for x in combinations(chain(*freq_n), 2)] 
        for i in range(len(freq_nn)):
            freq_nn[i] = freq_nn[i][0] + freq_nn[i][1]
   
        indices2 = []
        for i, y in enumerate(freq):
            ind = [j for j,x in enumerate(freq_nn) if x%y==0.]# or x/y==3.]
            indices2.append(ind)
        indices2 = [x for x in chain(*indices2)]    
        freq_nn = np.delete(freq_nn, indices2)

        frequencies= [[] * x for x in range(3)]
        frequencies[0].append(freq_neg[0])
        frequencies[0].append([0])
        frequencies[0].append(freq_n[0])
        frequencies[0]=[x for x in chain(*frequencies[0])]

        for i in range(1,n-1):
            frequencies[1].append(freq_neg[i])
            frequencies[1].append(freq_n[i])
        frequencies[1]=[x for x in chain(*frequencies[1])]

        frequencies[2]=freq_nn    


        ##Franck-Condon factors
        E = freq / 8065.5
        f_n = [[] * i for i in range(n)]
        
        for j in range(0,n):
            f_n[j] = np.exp(-E * j /( kB * T))
        
        #partition function
        Z=np.empty(len(S))
        Z=np.sum(f_n,0)
        
        #occupation probability
        w_n=[[]*k for k in range(n)]
        for l in range(n):
            w_n[l]=f_n[l]/Z
        
        # overlap wavefunctions
        O_n = [[] * m for m in range(n)]
        O_neg = [[] * m for m in range(n)]
        for o in range(n):
            O_n[o] = [[] * p for p in range(n)]
            O_neg[o] = [[] * p for p in range(n - 1)]
            for q in range(o, n + o):
                a = np.minimum(o, q)
                summe=[]
                for k in range(a+1):
                    s = ((-1)**(q - k) * np.sqrt(S)**(o + q - 2 * k) *
                         factorial(o) * factorial(q) / 
                         (factorial(k) * factorial(o - k) * factorial(q - k)))
                    summe.append(s)
                summe = np.sum(summe, 0)
                O_n[o][q-o] = (np.exp(-S / 2) / 
                               (factorial(o) * factorial(q))**(0.5) * 
                               summe)**2 * w_n[o]
            for q in range(n - 1):
                O_neg[o][q] = [0 * b for b in range(len(S))]
            for q in range(o-1, -1, -1):
                a = np.minimum(o,q)
                summe = []
                for k in range(a+1):
                    s=((-1)**(q - k) * np.sqrt(S)**(o + q - 2 * k) * 
                       factorial(o) * factorial(q) / 
                       (factorial(k) * factorial(o - k) * factorial(q - k)))
                    summe.append(s)
                summe = np.sum(summe, 0)
                O_neg[o][q] = (np.exp(-S / 2) / 
                               (factorial(o) * factorial(q))**(0.5) * 
                               summe)**2 * w_n[o]
        O_neg = np.delete(O_neg, 0, 0)

        #Franck-Condon factors
        FC_n=[[]*i for i in range(n)]
        FC_n=np.sum(O_n,0)
        zero=reduce(mul,FC_n[0])
        FC_neg=[[]*i for i in range(n-2)]
        FC_neg=np.sum(O_neg,0)    
        FC_n=np.delete(FC_n,0,0)

        #combination FC factors
        FC_nn=[x for x in combinations(chain(*FC_n),2)]
        for i in range(len(FC_nn)):
            FC_nn[i]=FC_nn[i][0]*FC_nn[i][1]

        FC_nn=np.delete(FC_nn, indices2)

        FC=[[]*x for x in range(3)]
        FC[0].append(FC_neg[0])
        FC[0].append([zero])
        FC[0].append(FC_n[0])
        FC[0]=[x for x in chain(*FC[0])]

        for i in range(1,n-1):
            FC[1].append(FC_neg[i])
            FC[1].append(FC_n[i])
        FC[1]=[x for x in chain(*FC[1])]

        FC[2]=FC_nn    
 
        """Returned are two 3-dimensional lists. First inner list contains 
frequencies and FC-factors of vibrations exited with |1| quanta and 
the 0-0 transition.
        Second list contains frequencies and FC-factors from higher 
quanta exitations. Third list are combinations of two normal modes 
(including combinations of higher quanta exitations). """ 
        return FC, frequencies
