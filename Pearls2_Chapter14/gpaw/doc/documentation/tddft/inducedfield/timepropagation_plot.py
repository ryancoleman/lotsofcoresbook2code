# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt

from gpaw.inducedfield.inducedfield_base import read_data


# Helper function
def do_plot(d_g, ng, box, atom_a):
    # Take slice of data array
    d_yx = d_g[:, ng[1] // 2, :]
    y = np.linspace(0, box[0], ng[0])
    ylabel = u'x / Å'
    x = np.linspace(0, box[2], ng[2])
    xlabel = u'z / Å'
    
    # Plot
    plt.figure()
    ax = plt.subplot(1, 1, 1)
    X, Y = np.meshgrid(x, y)
    plt.contourf(X, Y, d_yx, 40)
    plt.colorbar()
    for atom in atom_a:
        pos = atom['pos']
        plt.scatter(pos[2], pos[0], s=50, c='k', marker='o')
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.xlim([x[0], x[-1]])
    plt.ylim([y[0], y[-1]])
    ax.set_aspect('equal')


# Read InducedField object as standard numpy arrays, so
#  - no GPAW's GridDescriptor is required
#  - this code is not parallel safe
data = read_data('na2_td_field.ind')

# Choose array
w = 1                           # Frequency index
freq = data['freq_w'][w]        # Frequency
box = np.diag(data['cell_cv'])  # Calculation box
d_g = data['Ffe_wg'][w]         # Data array
ng = d_g.shape                  # Size of grid
atom_a = data['atom_a']         # Atoms

do_plot(d_g, ng, box, atom_a)
plt.title('Field enhancement @ %.2f eV' % freq)
plt.savefig('na2_td_Ffe.png', bbox_inches='tight')

# Imaginary part of density
d_g = data['Frho_wg'][w].imag
ng = d_g.shape
do_plot(d_g, ng, box, atom_a)
plt.title('Imaginary part of induced charge density @ %.2f eV' % freq)
plt.savefig('na2_td_Frho.png', bbox_inches='tight')

# Imaginary part of potential
d_g = data['Fphi_wg'][w].imag
ng = d_g.shape
do_plot(d_g, ng, box, atom_a)
plt.title('Imaginary part of induced potential @ %.2f eV' % freq)
plt.savefig('na2_td_Fphi.png', bbox_inches='tight')
