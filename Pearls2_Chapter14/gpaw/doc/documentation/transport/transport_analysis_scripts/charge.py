from __future__ import print_function
from gpaw.transport.analysor import Transport_Plotter
import numpy as np
from pylab import *
import sys
from ase.data import chemical_symbols

fd=0
bias_steps = sys.argv[1]

orbital_type = None
plotter = Transport_Plotter(fd)
plotter.plot_setup()
plotter.read_overhead()

legends = []
flags = ['b-o', 'r-o', 'g-o' ,'c-o', 'y-o', 'm-o', 'k-o']
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
    #print atom_indices
    if orbital == 'A':
        orbital_type = None
    else:
        orbital_type = orbital
    biases = []
    charges = []
    for bs in range(int(bias_steps)):
        try:
            charge = plotter.charge(bs, 0, atom_indices, orbital_type)    
            bias = plotter.get_info('bias', bs, 0)
            charges.append(charge)
            biases.append(bias[0]-bias[1])
        except IOError:
            print(' no file for bias_step', bs)
    biases = np.array(biases)       
    charges = np.array(charges)
    plot(biases, charges, flags[i])
    dense_level=1
    if dense_level>1:
        from scipy import interpolate
        tck = interpolate.splrep(biases, charges, s=0)
        numb = len(biases)
        newbiases = np.linspace(biases[0], biases[-1], numb * (dense_level))
        newcharges = interpolate.splev(newbiases, tck, der=0)
        biases = newbiases
        charges = newcharges
        plot(biases, charges, 'r-o')

    legends.append(item)

legend(legends)
xlabel('Bias(V)')
ylabel('charge(Electron)')
show()
            



