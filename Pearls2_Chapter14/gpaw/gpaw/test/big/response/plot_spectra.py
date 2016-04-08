import numpy as np
from gpaw.test import wrap_pylab
wrap_pylab()
import pylab as pl

def plot_EELS(head):
# plot EELS spectra
    pl.figure(figsize=(4,7))

    d = np.loadtxt(head + '_q_list')
    q = d[:,0]
    ndata = q.shape[0] + 1
    w = np.zeros(ndata)
    w2 = np.zeros(ndata)
    
    for i in range(1,ndata):
        filename = head + '_EELS_' + str(i) 

        d = np.loadtxt(filename)
        pl.plot(d[:,0], d[:,2] + 0.4*(ndata-i-1),'-k', label=str(q[i-1])[:4])
        
    fsize = 14
    pl.xlabel('Energy (eV)', fontsize=fsize)
    pl.ylabel('Loss Function', fontsize=fsize)
    pl.title('EELS spectra of ' + head +': $\Gamma-\mathrm{M}$', fontsize=fsize)
    pl.ylim(0,)
    pl.legend(loc=2)
    
def plot_ABS(head):    
    # plot optical absorption specctrum
    pl.figure()
    d = np.loadtxt(head+'_abs.z')
    pl.plot(d[:,0], d[:,1], '-k', label='$\mathrm{Re}\epsilon(\omega)$')
    pl.plot(d[:,0], d[:,2], '-r', label='$\mathrm{Im}\epsilon(\omega)$')

    if head == 'si':
        # data from G.Kresse, PRB 73, 045112 (2006) 
        x = np.array([2.53, 2.71, 3.08, 3.72, 4.50])
        # mine
        x2 = np.array([2.53,2.75, 3.05, 3.67, 4.48])
        y = np.array([1.5,  18,  24,   47,   28  ])
        pl.plot(x,y,'or')

    fsize = 14
    pl.title('Dielectric function of ' + head)
    pl.legend()
    pl.xlabel('Energy (eV)', fontsize=fsize)
    pl.ylabel('$\epsilon$', fontsize=fsize)
    pl.xlim(0,10)
    
plot_EELS('graphite')
plot_ABS('si')


pl.show()
