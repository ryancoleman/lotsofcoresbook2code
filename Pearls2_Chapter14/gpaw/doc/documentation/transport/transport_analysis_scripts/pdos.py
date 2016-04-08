from gpaw.transport.analysor import Transport_Plotter
import numpy as np
from pylab import *
import sys
from ase.data import chemical_symbols
if '*' in sys.argv[1]:
    fd=0
    bias_step = sys.argv[1].split('*')[0]
else:
    fd=1
    bias_step = sys.argv[1]

orbital_type = None
plotter = Transport_Plotter(fd)
plotter.plot_setup()
plotter.read_overhead()

legends = []
flags = ['b-', 'r-', 'g-' ,'c-', 'y-', 'm-', 'k-']

energies = np.linspace(-5, 5, 201)
for i, item in enumerate(sys.argv[2:]):
    group, orbital = item.split('_')
    elements = group.split('-')   

    if len(elements) == 1 and elements[0] == 'A':
        atom_indices = None
    else:
        atom_indices = []
        for element in elements:
            big = len(plotter.atoms)
            small = -1
            equal = None
            x, y, z = 999, 999, 999
            if ']' in element:
                tmp = element.split(']')
                x = element.index(']')
                tmp = tmp[1].split('[')[0]
                tmp = tmp.split('=')[0]
                small = int(tmp)
            if '[' in element:
                tmp = element.split('[')
                y = element.index('[')
                tmp = tmp[1].split(']')[0]
                tmp = tmp.split('=')[0]
                big = int(tmp)
            if '=' in element:
                tmp = element.split('=')
                z = element.index('=')
                tmp = tmp[1].split('[')[0]
                tmp = tmp.split(']')[0]
                equal = int(tmp)
            min = np.min([x,y,z])
            if min == 999:
                pass
            else:
                element = element[:min]
            for j, atom in enumerate(plotter.atoms):
                if atom.symbol == element and j > small and  j < big:
                    if equal is None or (equal is not None and j == equal):
                        atom_indices.append(j)
    if orbital == 'A':
        orbital_type = None
    else:
        orbital_type = orbital
    pdos = plotter.partial_dos(bias_step, 0, 0, None, atom_indices, orbital_type) 
    plot(energies, pdos, flags[i])
    dense_level=1
    if dense_level>1:
        from scipy import interpolate
        tck = interpolate.splrep(energies, pdos, s=0)
        numb = len(energies)
        newee = np.linspace(energies[0], energies[-1], numb * (dense_level))
        newpdos = interpolate.splev(newee, tck, der=0)
        ee = newee
        pdos = newpdos
        plot(ee, pdos, flags[i])

    legends.append(item)

legend(legends)
bias = plotter.get_info('bias', bias_step)
eye = np.zeros([10, 1]) + 1
f1 = bias[0] * eye
f2 = bias[1] * eye        
a1 = np.max(pdos)
l1 = np.linspace(0, a1, 10)

plot(f1, l1, 'r--')
plot(f2, l1, 'r--')

xlabel('Energy(eV)')
ylabel('Partial Density of States(Electron/eV)')
show()

