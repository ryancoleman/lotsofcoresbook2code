from __future__ import print_function
import numpy as np
from pylab import *
from gpaw.response.tool import linear_fit
from gpaw.test import equal

tag = 'Nabulk'
d = np.loadtxt('%s_RPA.dat' %(tag))
x = d[:,0]
y = d[:,1]

xscale = x**(-1.5)

plot(xscale, y, 'o')

yfit, para = linear_fit(xscale,y)

xx = np.linspace(0, xscale.max(), 1000)
yy = para[0] * xx + para[1]
plot(xx, yy, '-k')

t = [int(x[i]) for i in range(len(d))]

xticks(xscale, t)
axis([0, xscale.max(), None, None])

xlabel('Response function cutoff (eV)')
ylabel('RPA correlation energy (eV)')
print('RPA correlation at infinite cutoff:', para[1])

equal(para[1], -1.27, 0.01)
#show()

