
import numpy as np
from pylab import *
from gpaw.test import equal

output = True

############################# H2O ############################

filename = 'H2O'

data = loadtxt(filename+'-td'+'.dat')
be_t = data[:,0]
f_t = data[:,1]

data = loadtxt(filename+'-dos'+'.dat')
be_d = data[:,0]
f_d = data[:,1]

if output:
    fig = figure(1, figsize=(14,8))
    fig.subplots_adjust(left=0.08, right=0.98)

    subplot(3, 1, 1)
    bar(be_d-0.1, f_d/f_d.max(), color='k', width=0.2, label='DFT eigenvalue')
    bar(be_t-0.1, f_t/f_t.max(), color='r', width=0.2, label='lr-TDDFT method')

    data = loadtxt(filename+'.dat')
    be_e = data[:,0]
    f_e = 2.1*data[:,1]

    plot(be_e, f_e, '-b', label='Exp')
    text(6.5, 0.8, 'H$_2$O', fontsize=22)
    #xlabel('Binding energy [eV]', fontsize=20)
    ylabel('Spectroscopic factor', fontsize=20)
    axis((5,45,0,2))
    xticks(fontsize=18)
    yticks(fontsize=18)

etol = 1.0
peaks = np.array([12.8, 14.9, 18.5])
ube_t = np.unique(np.round(be_t,1))[:len(peaks)]
ube_d = np.unique(np.round(be_d,1))[:len(peaks)]

for t, d, ref in zip(ube_t, ube_d, peaks):
    equal(t, ref, etol)
    equal(t, d, etol)

############################# CO ############################

filename='CO'

data = loadtxt(filename+'-td'+'.dat')
be_t = data[:,0]
f_t = data[:,1]

data = loadtxt(filename+'-dos'+'.dat')
be_d = data[:,0]
f_d = data[:,1]

if output:
    subplot(3, 1, 2)
    bar(be_d-0.1, f_d/f_d.max(), color='k', width=0.2, label='DFT eigenvalue')
    bar(be_t-0.1, f_t/f_t.max(), color='r', width=0.2, label='lr-TDDFT method')

    data = loadtxt(filename+'.dat')
    be_e = data[:,0]
    f_e = 2.1*data[:,1]

    plot(be_e, f_e, '-b', label='Exp')
    text(6.5,0.8,'CO',fontsize=22)
    #xlabel('Binding energy [eV]', fontsize=20)
    #ylabel('Spectroscopic factor', fontsize=20)
    axis((5,45,0,2))
    xticks(fontsize=18)
    yticks(fontsize=18)

etol = 1.0
peaks = np.array([13.7, 17.7, 19.0])
ube_t = np.unique(np.round(be_t,1))[:len(peaks)]
ube_d = np.unique(np.round(be_d,1))[:len(peaks)]

for t, d, ref in zip(ube_t, ube_d, peaks):
    equal(t, ref, etol)
    equal(t, d, etol)

############################# NH3 ############################

filename = 'NH3'

data = loadtxt(filename+'-td'+'.dat')
be_t = data[:,0]
f_t = data[:,1]

data = loadtxt(filename+'-dos'+'.dat')
be_d = data[:,0]
f_d = data[:,1]

if output:
    subplot(3, 1, 3)
    bar(be_d-0.1, f_d/f_d.max(), color='k', width=0.2, label='DFT eigenvalue')
    bar(be_t-0.1, f_t/f_t.max(), color='r', width=0.2, label='lr-TDDFT method')

    data = loadtxt(filename+'.dat')
    be_e = data[:,0]
    f_e = 2.1*data[:,1]

    plot(be_e, f_e, '-b', label='Exp')
    leg = legend()
    text(6.0, 0.8, 'NH$_3$', fontsize=22)
    xlabel('Binding energy [eV]', fontsize=20)
    #ylabel('Spectroscopic factor', fontsize=20)
    axis((5,45,0,2))
    xticks(fontsize=18)
    yticks(fontsize=18)

    savefig('PES_fig.png')
    #show()

etol = 1.0
peaks = np.array([10.8, 16.2])
ube_t = np.unique(np.round(be_t,1))[:len(peaks)]
ube_d = np.unique(np.round(be_d,1))[:len(peaks)]

for t, d, ref in zip(ube_t, ube_d, peaks):
    equal(t, ref, etol)
    equal(t, d, etol)

