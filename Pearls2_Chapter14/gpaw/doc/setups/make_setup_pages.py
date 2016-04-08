# -*- coding: utf-8 -*-
from __future__ import print_function
import os
import sys
import pickle

import matplotlib
matplotlib.use('Agg')

# Workaround for bug #562421 in python-matplotlib on Fedora
# https://bugzilla.redhat.com/show_bug.cgi?id=562421
matplotlib.rc('text', usetex=True)

import numpy as np
from ase.data import atomic_numbers, atomic_names
from ase.atoms import string2symbols
from ase.data.molecules import rest
from ase.data.g2_1 import data as molecule_data

from ase.data.g2_1_ref import atomization_vasp
from gpaw.atom.configurations import parameters

page = """.. Computer generated reST (make_setup_pages.py)
.. index:: %(name)s
.. _%(name)s:

================
%(name)s
================


Tests
=====

%(tests)s


Convergence tests
=================

The energy of %(aname)s dimer (`E_d`) and %(aname)s atom (`E_a`) is
calculated at diferent grid-spacings (`h`).

.. image:: ../_static/setups-data/%(symbol)s-dimer-eggbox.png


Setup details
=============

The setup has %(Nv)d valence electron%(plural)s and %(Nc)d electrons
in the frozen core.  It is based on a scalar-relativistic spin-paired
neutral all-electron PBE calculation.

Cutoffs and eigenvalues:

===  ===  =============  ==========
id   occ  eigenvals       cutoff
---  ---  -------------  ----------
%(table)s
===  ===  =============  ==========

Other cutoffs:

==========================  =====================
compensation charges        %(rcutcomp).2f Bohr
filtering                   %(rcutfilter).2f Bohr
core density                %(rcutcore).2f Bohr
==========================  =====================

Energy Contributions:

=========  ======================
Kinetic    %(Ekin).4f Ha
Potential  %(Epot).4f Ha
XC         %(Exc).4f Ha
---------  ----------------------
Total      %(Etot).4f Ha
=========  ======================


Wave functions, projectors, ...
-------------------------------

.. image:: ../_static/setups-data/%(symbol)s-setup.png



Back to :ref:`setups`.


"""

def make_page(symbol):
    filename = symbol + '.rst'
    try:
        data = pickle.load(open(
            '../_static/setups-data/%s.pckl' % symbol, 'rb'))
    except EOFError:
        print(symbol, 'missing!')
        return

    Z = atomic_numbers[symbol]
    name = atomic_names[Z]

    bulk = []
    if symbol in ('Ni Pd Pt La Na Nb Mg Li Pb Rb Rh Ta Ba Fe Mo C K Si W V ' +
                  'Zn Co Ag Ca Ir Al Cd Ge Au Cs Cr Cu').split():
        bulk.append(symbol)
        
    for alloy in ('LaN LiCl MgO NaCl GaN AlAs BP FeAl BN LiF NaF ' +
                  'SiC ZrC ZrN AlN VN NbC GaP AlP BAs GaAs MgS ' +
                  'ZnO NiAl CaO').split():
        if symbol in alloy:
            bulk.append(alloy)
            
    if len(bulk) > 0:
        if len(bulk) == 1:
            bulk = 'test for ' + bulk[0]
        else:
            bulk = 'tests for ' + ', '.join(bulk[:-1]) + ' and ' + bulk[-1]
        tests = 'See %s here: :ref:`bulk_tests`.  ' % bulk
    else:
        tests = ''
        
    molecules = []
    for x in atomization_vasp:
        if symbol in x:
            molecules.append(molecule_data[x]['name'])
    if molecules:
        names = [rest(m) for m in molecules]
        if len(molecules) == 1:
            mols = names[0]
        elif len(molecules) > 5:
            mols = ', '.join(names[:5]) + ' ...'
        else:
            mols = ', '.join(names[:-1]) + ' and ' + names[-1]
        tests += 'See tests for %s here: :ref:`molecule_tests`.' % mols

    if name[0] in 'AEIOUY':
        aname = 'an ' + name.lower()
    else:
        aname = 'a ' + name.lower()

    table = ''
    for n, l, f, e, rcut in data['nlfer']:
        if n == -1:
            n = '\*'
        table += '%2s%s  %3d  %10.3f Ha' % (n, 'spdf'[l], f, e)
        if rcut:
            table += '  %.2f Bohr\n' % rcut
        else:
            table += '\n'

    f = open(symbol + '.rst', 'w')
    f.write(page % {
        'name': name, 'aname': aname, 'tests': tests, 'symbol': symbol,
        'Nv': data['Nv'], 'Nc': data['Nc'], 'plural': 's'[:data['Nv'] > 1],
        'table': table,
        'rcutcomp': data['rcutcomp'],
        'rcutfilter': data['rcutfilter'],
        'rcutcore': data['rcore'],
        'Ekin': data['Ekin'], 'Epot': data['Epot'], 'Exc': data['Exc'],
        'Etot': data['Ekin'] + data['Epot'] + data['Exc']})
    f.close()

    # Make convergence test figure:
    d = data['d0'] * np.linspace(0.94, 1.06, 7)
    h = data['gridspacings']
    ng = len(h)
    Edimer0 = np.empty(ng)
    ddimer0 = np.empty(ng)
    Eegg = data['Eegg']
    Edimer = data['Edimer']
    for i in range(ng):
        E = Edimer[i]
        energy = np.polyfit(d**-1, E, 3)
        der = np.polyder(energy, 1)
        roots = np.roots(der)
        der2 = np.polyder(der, 1)
        if np.polyval(der2, roots[0]) > 0:
            root = roots[0]
        else:
            root = roots[1]
        if isinstance(root, complex):
            print('??????')
            root = root.real
        d0 = 1.0 / root.real
        E0 = np.polyval(energy, root)
        Edimer0[i] = E0
        ddimer0[i] = d0

    if not (d[0] < d0 < d[-1]):
        print(d, d0, symbol)

    print('%2s %.3f %+7.1f %%' % (symbol, d0, 100 * (d0 / d[3] - 1)))
    
    Ediss = 2 * Eegg[:, 0] - Edimer0

    import pylab as plt
    dpi = 80
    fig = plt.figure(figsize=(6, 11), dpi=dpi)
    fig.subplots_adjust(left=0.15, right=0.97, top=0.97, bottom=0.04)

    plt.subplot(311)
    plt.semilogy(h[:-1], 2 * abs(Eegg[:-1, 0] - Eegg[-1, 0]), '-o',
                 label=r'$2E_a$')
    plt.semilogy(h[:-1], abs(Edimer0[:-1] - Edimer0[-1]), '-o',
                 label=r'$E_d$')
    plt.semilogy(h[:-1], abs(Ediss[:-1] - Ediss[-1]), '-o',
                 label=r'$2E_a-E_d$')
    plt.semilogy(h, Eegg.ptp(axis=1), '-o', label=r'$E_{egg}$')
    #plt.title('Energy differences')
    plt.xlabel(r'h [\AA]')
    plt.ylabel('energy [eV]')
    plt.legend(loc='best')

    plt.subplot(312)
    plt.semilogy(h, np.abs(data['Fegg']).max(axis=1), '-o',
                 label=r'$|\mathbf{F}_{egg}|$')
    plt.semilogy(h, np.abs(data['Fdimer'].sum(axis=2)).max(axis=1), '-o',
                 label=r'$|\mathbf{F}_1 + \mathbf{F}_2|$')
    #plt.title('Forces')
    plt.xlabel(r'h [\AA]')
    plt.ylabel(r'force [eV/\AA]')
    plt.legend(loc='best')

    plt.subplot(313)
    d = ddimer0[-1]
    plt.plot(h, ddimer0, '-o')
    #plt.title('Bond length')
    plt.xlabel(r'h [\AA]')
    plt.axis(ymin=d * 0.98, ymax=d * 1.02)
    plt.ylabel(r'bond length [\AA]')

    plt.savefig('../_static/setups-data/%s-dimer-eggbox.png' % symbol, dpi=dpi)
    #plt.show()

args = sys.argv[1:]
if len(args) == 0:
    args = parameters.keys()
for symbol in args:
    make_page(symbol)
