from __future__ import print_function
from ase.structure import molecule
from ase.parallel import paropen
from gpaw import GPAW
from gpaw.utilities.tools import split_formula
from gpaw.test import equal

cell = [10.,10.,10.]
data = paropen('data.txt', 'w')

##Reference from J. Chem. Phys. Vol 120 No. 15, 15 April 2004, page 6898
tpss_de = {
'Li2': 22.5,
}
tpss_old = {
'Li2': 22.7,
}

exp_bonds_dE = {
'Li2': (2.673,24.4),
}

niters_ref = {'Li2': 21, 'Li': 14}
niter_tolerance = 0

systems = ['Li2']

# Add atoms
for formula in systems:
    temp = split_formula(formula)
    for atom in temp:
        if atom not in systems:
            systems.append(atom)
energies = {}
niters = {}

# Calculate energies
for formula in systems:
    loa = molecule(formula)
    loa.set_cell(cell)
    loa.center()
    calc = GPAW(h=0.3,
                nbands=-2,
                xc='PBE',
                #fixmom=True,
                txt=formula + '.txt')
    if len(loa) == 1:
        calc.set(hund=True)
    else:
        pos = loa.get_positions()
        pos[1,:] = pos[0,:] + [0.0, 0.0, exp_bonds_dE[formula][0]]
        loa.set_positions(pos)
        loa.center()
    loa.set_calculator(calc)
    try:
        energy = loa.get_potential_energy()
        niters[formula] = calc.get_number_of_iterations()
        diff = calc.get_xc_difference('TPSS')
        energies[formula] = (energy, energy + diff)
    except:
        raise#print >> data, formula, 'Error'
    else:
        print(formula, energy, energy + diff, file=data)
    data.flush()

#calculate atomization energies
file = paropen('tpss.txt', 'w')
print('formula\tGPAW\tRef\tGPAW-Ref\tGPAW-exp', file=file)
mae_ref, mae_exp, mae_pbe, count = 0.0, 0.0, 0.0, 0
for formula in tpss_de.keys():
    try:
        atoms_formula = split_formula(formula)
        de_tpss = -1.0 * energies[formula][1]
        de_pbe = -1.0 * energies[formula][0]
        for atom_formula in atoms_formula:
            de_tpss += energies[atom_formula][1]
            de_pbe += energies[atom_formula][0]
    except:
        raise#print >>file, formula, 'Error'
    else:
        de_tpss *= 627.5/27.211
        de_pbe *= 627.5/27.211
        mae_ref += abs(de_tpss-tpss_de[formula])
        mae_exp += abs(de_tpss-exp_bonds_dE[formula][1])
        mae_pbe += abs(de_pbe-exp_bonds_dE[formula][1])
        count += 1
        out = "%s\t%.1f\t%.1f\t%.1f\t%.1f kcal/mol"%(formula,de_tpss,tpss_de[formula],
                                            de_tpss-tpss_de[formula],de_tpss-exp_bonds_dE[formula][1])
        print(out, file=file)
        file.flush()


#comparison to gpaw revision 5450 version value in kcal/mol (note the grid:0.3 Ang)
    equal(de_tpss, tpss_old[formula], 0.15)
