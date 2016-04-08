import StringIO
from itertools import chain

import numpy as np

import matplotlib.pyplot as plt
import matplotlib 

from ase.vibrations import Vibrations
from ase.vibrations.franck_condon import FranckCondon
from ase.io.xyz import read_xyz
from ase.optimize import FIRE

from gpaw import GPAW
from gpaw.cluster import Cluster
from gpaw.utilities.folder import Folder

butadiene = """10

C       3.649801161546418      5.442281389577507      3.863313703750026
C       5.051651240044169      5.368220758269772      4.162165876906096
C       5.750174626862403      4.162261915959347      4.240449977068684
C       7.150130182125531      4.155384186721486      4.537328602062397
H       3.218154657585170      4.565210696328925      3.522601038049320
H       3.077656122062729      6.375092902842770      3.826039498180272
H       5.478464901706067      6.370680001794822      4.422235395756437
H       5.320549047980879      3.220584852467720      3.974551561510350
H       7.723359150977955      3.224855971783890      4.574146712279462
H       7.580803493981530      5.034479218283977      4.877211530909463
"""

h = 0.3
atoms = Cluster(read_xyz(StringIO.StringIO(butadiene)))
atoms.minimal_box(3., h)
atoms.set_calculator(GPAW(h=h))
if 0:
    dyn = FIRE(atoms)
    dyn.run(fmax=0.05)
    atoms.write('butadiene.xyz')

vibname = 'fcvib'
vib = Vibrations(atoms, name=vibname)
vib.run()

# Modul
a = FranckCondon(atoms, vibname, minfreq=250)

# excited state forces
F = np.array([
        [-2.11413,    0.07317,   -0.91682],
        [3.23569,   -0.74520,    0.76758],
        [-3.44847,    0.63846,  -0.81080],
        [2.77345,    0.01272,    0.74811],
        [-0.06544,   -0.01078,   -0.03209],
        [-0.01245,   -0.01123,   -0.00040],
        [0.00186,   -0.05864,   -0.00371],
        [-0.00151,    0.05815,    0.00141],
        [0.01625,    0.00781,   -0.00202],
        [0.06253,   0.00902,    0.03381]
        ])

# Huang-Rhys factors
S, fq = a.get_Huang_Rhys_factors(F)

#Temperature and #quanta taken into account
T = 300
n = 10

colors=['red','blue','green']
labels=['normal modes with #quanta=1,-1 and the 0-0 transition', 'normal modes with #quanta>1,<-1', 'combination of two normal modes']

S = np.append(1, S)
fq = np.append(0, fq)
plt.vlines(fq,0,S,color='darkblue',label='Huang-Rhys factors', linewidth=2.0)

plt.legend(loc='upper right')
plt.axis([-1000,6000,0,1.05])
plt.ylabel('HR factors [a.u.]')
plt.xlabel('frequency [cm-1]')
plt.show()

#Plot Franck-Condon factors
FC, f = a.get_Franck_Condon_factors(n,T,F)
for l in range(len(FC)):
    plt.vlines(f[l],0,FC[l],color=colors[l],label=labels[l],linewidth=2.0)

#Fold the spectrum with a gaussian function
width = 300
folding = 'Gauss'
x=[j for j in chain(*f)]
y=[j for j in chain(*FC)]

X, Y = Folder(width, folding).fold(x, y)
plt.plot(X, Y*750,color='black',linewidth=2.0,label='Temperature='+str(T)+'K, #quanta='+str(n)+', width of gaussian='+str(width)+'cm^-1',linestyle='dashed')

plt.legend(loc='upper right')
plt.axis([-1000,6000,0,0.6])
plt.ylabel('FC factors [a.u.]')
plt.xlabel('frequency [cm-1]')
plt.show()


