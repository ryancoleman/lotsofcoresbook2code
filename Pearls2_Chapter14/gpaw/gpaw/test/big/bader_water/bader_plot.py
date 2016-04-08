#*******************PYLAB STUFF************************************************
from matplotlib import rc # This must be before "imort pylab"
rc('font', **{'family':'serif','sans-serif':['cm'],'serif':['cm']})
rc('text', usetex=True)
import pylab as pl

# this test required texlive-latex

width_pt = 360.0               # Get this from LaTeX using \showthe\columnwidth
inch = 1.0 / 72.27             # Convert pt to inch
golden_mean = (pl.sqrt(5)-1.0) / 2.0     # Aesthetic ratio
width_inch = width_pt * inch             # width in inches
height_inch = width_inch * golden_mean   # height in inches
figsize = [width_inch, height_inch]
params = {'lines.linewidth': 1.2,
          'font.size': 10,
          'figure.figsize': figsize,
          'figure.dpi': 200,
          'savefig.dpi': 200,
          'figure.subplot.right': 1.0,
          'figure.subplot.top': 0.88,
          'figure.subplot.left': 0.0,
          'figure.subplot.bottom': 0.0,
          'figure.subplot.wspace': 0.0,
          'figure.subplot.hspace': 0.0,
          'font.size': 10.0,
          'legend.fontsize': 'medium',
          'legend.loc': 'upper right',
          }
pl.rcParams.update(params)
#*******************PYLAB STUFF************************************************

from ase.io.cube import read_cube
import numpy as np

# Contour values and colors
vals = np.linspace(0.1, 1.8, 16)**2
vals = list(vals) + [99,]
colors = pl.cm.Reds(np.linspace(0.15, 1.0, len(vals)))

#Pseudo density
nt, atoms = read_cube('water_pseudo_density.cube', read_data=True)
x = len(nt) // 2
nt = nt[x]

# All electron density and bader volumes
n, atoms = read_cube('water_density.cube', read_data=True)
#bader, atoms2 = read_cube('AtIndex.cube', read_data=True)
x = len(n) // 2
n = n[x]
#bader = bader[x]

# plot
fig = pl.figure(figsize=(6.2, 3))
pl.subplot(121)
pl.contourf(nt.T, vals, origin='lower', extend='neither', colors=colors)
pl.axis('equal')
pl.axis([52-25, 52+25, 52-25, 52+25])
pl.axis('off')
pl.text(52.5, 55, '$7.07e$', size=20, ha='center', va='center')
pl.title('Pseudo density', size=20)

pl.subplot(122)
pl.contourf(n.T, vals, colors=colors,
           origin='lower', extend='neither')
#pl.contour(bader.T, [1.5], origin='lower', extend='neither', colors='k')
pl.axis('equal')
pl.axis([104-50, 104+50, 104-50, 104+50])
pl.axis('off')
pl.text(104.0, 112.0, '$9.12e$', size=20, ha='center', va='center')
pl.text( 86.5,  97.5, '$0.44e$', size=20, ha='right', va='center')
pl.text(122.0,  97.5, '$0.44e$', size=20, ha='left', va='center')
pl.title('All-electron density', size=20)

pl.savefig('water_divide_surf.eps')
pl.show()
