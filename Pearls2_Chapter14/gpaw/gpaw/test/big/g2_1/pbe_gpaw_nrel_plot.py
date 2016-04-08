from __future__ import print_function
import os
import warnings
# silence matplotlib.use() warning
warnings.filterwarnings('ignore', '.*This call to matplotlib\.use.*',)

import csv

import numpy as np

import heapq

import matplotlib.pyplot as plt

from ase.data.molecules import latex

ann_fontsize = 'small'
label_fontsize = 12

from ase.data import atomic_numbers

from ase.structure import molecule
from ase.data.g2_1 import data

from gpaw.test.big.g2_1.pbe_nwchem_def2_qzvppd_analyse import ref

from gpaw.test.big.g2_1.pbe_gpaw_nrel_analyse import ref as calc
from gpaw.test.big.g2_1.pbe_gpaw_nrel_analyse import tag

formulas = ref['ea'].keys()
formulas.sort()

atoms = []
for f in formulas:
    m = molecule(f, data=data)
    atoms.extend(m.get_chemical_symbols())
atoms = list(set(atoms))
anum = [atomic_numbers[a] for a in atoms]

# sort according to anum
tmpsort = zip(anum, atoms)
tmpsort.sort()
a, atomssorted = zip(*tmpsort)

def containatom(atom, formulas, data):
    molecules = []
    for f in formulas:
        m = molecule(f, data=data)
        if atom in list(set(m.get_chemical_symbols())):
            molecules.append(f)
    return list(set(molecules))

def get_statistics(result, reference):
    skeys = result.keys()
    skeys.sort()
    res1 = []
    res2 = []
    for k in skeys:
        res1.append(k)
        res2.append(result[k])
    skeys = reference.keys()
    skeys.sort()
    ref1 = []
    ref2 = []
    for k in skeys:
        ref1.append(k)
        ref2.append(reference[k])
    resrm = []
    refrm = []
    # remove result without reference and reference without result
    toremove = list(set(res1).symmetric_difference(set(ref1)))
    for c in toremove:
        if c in res1:
            print(result + ' vs ' + reference + ': ' + result + ' removed ' + c)
            i = res1.index(c)
            res1.pop(i)
            res2.pop(i)
            resrm.append(c)
        if c in ref1:
            print(result + ' vs ' + reference + ': ' + reference + ' removed ' + c)
            i = ref1.index(c)
            ref1.pop(i)
            ref2.pop(i)
            refrm.append(c)
    # calculate result - reference
    diff = np.array(res2) - np.array(ref2)
    all = [(ref1[list(diff).index(e)], e) for e in diff]
    average = np.average(diff)
    absaverage = np.average(abs(diff))
    std = np.std(diff)
    nlargest = heapq.nlargest(4, abs(diff))
    # n largest errors
    largest = []
    for abse in nlargest:
        c = ref1[list(abs(diff)).index(abse)]
        e = diff[ref1.index(c)]
        largest.append((c, e))
    return average, absaverage, std, largest, all, resrm, refrm

# prepare plot

vs = [
    (calc['ea'].copy(), ref['ea'].copy(), 'G2_1'),
    ]

for atom in atoms:
    calc_atom = {}
    for a in containatom(atom, formulas, data):
        calc_atom.update({a: calc['ea'][a]})
    ref_atom = {}
    for a in containatom(atom, formulas, data):
        ref_atom.update({a: ref['ea'][a]})
    vs.append((calc_atom, ref_atom, atom))

no = []
average = []
absaverage = []
std = []
largest = []
all = []
resrm = []
refrm = []
labelsbase = []
labels = []
for n, results in enumerate(vs):
    label = results[2]
    labelsbase.append(label)
    labels.append(label)
    a, absa, s, l, al, rs, rf = get_statistics(results[0], results[1])
    no.append(n)
    average.append(a)
    absaverage.append(absa)
    std.append(s)
    largest.append(l)
    resrm.append(rs)
    refrm.append(rf)
    all.append(al)
#print no
#print average
#print absaverage
#print std
#print largest


def plot(xdata, ydata, std,
         title,
         xlabel, ylabel,
         label, color, alpha,
         miny, maxy,
         num=1,
         ):
    import matplotlib
    #matplotlib.use('Agg')
    import pylab
    import matplotlib.font_manager

    # all goes to figure num
    pylab.figure(num=num, figsize=(9.5, 9))
    pylab.gca().set_position([0.10, 0.20, 0.85, 0.60])
    # let the plot have fixed y-axis scale
    ywindow = maxy - miny
    #pylab.gca().set_ylim(miny, maxy+ywindow/5.0)
    pylab.gca().set_ylim(miny, maxy)
    #pylab.plot(xdata, ydata, 'b.', label=label, color=color)
    #pylab.plot(xdata, ydata, 'b-', label='_nolegend_', color=color)
    pylab.bar(xdata, ydata, 0.9, label=label, color=color, alpha=alpha)
    t = pylab.title(title)
    # http://old.nabble.com/More-space-between-title-and-secondary-x-axis-td31722298.html
    t.set_y(1.05)
    pylab.xlabel(xlabel)
    pylab.ylabel(ylabel)
    prop = matplotlib.font_manager.FontProperties(size=12)
    leg = pylab.legend(loc='upper right', fancybox=True, prop=prop)
    leg.get_frame().set_alpha(0.5)
    #pylab.savefig(directory_name + os.path.sep + out_prefix +'.eps', format='eps')

def plot_save(directory_name, out_prefix):
    from os.path import exists
    assert exists(directory_name)
    import pylab

    pylab.savefig(directory_name + os.path.sep + out_prefix +'.png', bbox_inches='tight')

import matplotlib
matplotlib.use('Agg')
from matplotlib import pylab, ticker
# print scaling results

num=1
miny = -0.3
maxy = 0.3

plot(
    no, average, std,
    'Atomization energy: E(atoms) - E(molecule)',
    'Subsets contain',
    'Calculated - Reference [eV]',
    'Average (avg)',
    'blue',
    0.5,
    miny, maxy,
    num=1,
    )
plot(
    no, absaverage, std,
    'Atomization energy: E(atoms) - E(molecule)',
    'Subsets contain',
    'Calculated - Reference [eV]',
    'Average absolute (abs)',
    'green',
    0.2,
    miny, maxy,
    num=1,
    )
zero = [0.0 for i in range(len(no))]
pylab.plot(no, zero, 'k-', label='_nolegend_')
ay1=pylab.gca()
ay1.xaxis.set_ticks([n + 0.5 for n in no])
ay1.xaxis.set_ticklabels(labels)
ay1.yaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())

for label in ay1.get_xticklabels() + ay1.get_yticklabels():
    label.set_fontsize(label_fontsize)

def errorslocation(n, n1):
    return (n + 0.1, 0.25 - 0.05 * n1)
def formulaslocation(n, n1):
    return (n + 0.3, - 0.10 - 0.05 * n1)

for n in range(len(no)):
    label = ''
    for n1, v in enumerate(['average', 'absaverage', 'std']):
        l = {
            'average': '\navg\n',
            'absaverage': '\nabs\n',
            'std': '\nstd\n',
            }[v]
        value = {
            'average': average,
            'absaverage': absaverage,
            'std': std,
            }[v]
        label += l + ' ' + str(round(value[n], 3)) + '\n'
    pylab.annotate(label,
                   xy=(n + 0.0, 0.0),
                   xytext=errorslocation(n, n1),
                   arrowprops=None,
                   horizontalalignment='left', verticalalignment='center',
                   fontsize=ann_fontsize,
                   )

# plot compounds with largest errors
for n, l in enumerate(largest):
    for n1, (c, e) in enumerate(l):
        name = latex(c) + '\n'
        # matplotlib.pyparsing.ParseFatalException: Expected end of math '$'
        # $\rm{SiH}2_\rm{s}3\rm{B}1\rm{d}$ (at char 0), (line:1, col:1)
        name = name.replace('\\rm', '')
        label = name + ' ' + str(round(e, 2))
        pylab.annotate(label,
                       xy=(n + 0.05, e),
                       xytext=formulaslocation(n, n1),
                       arrowprops=dict(width=0.05, headwidth=5.0, facecolor='black', shrink=1.00),
                       horizontalalignment='left', verticalalignment='center',
                       fontsize=ann_fontsize,
                       )
#pylab.show()
plot_save(".", tag + '_ea_vs')
#pylab.close(1)
