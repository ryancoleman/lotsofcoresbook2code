# Compare H adsorption energies and magnetic moments for a 12-ZGNR
# It needs 'ZGNR12.traj', 'ZGNR12_H.traj' and 'H2.traj'

from ase.io import read

reference = {'e_ads':-1.5,    #Reference from rev
             'magmom_C':0.5,
             'magmom_H':0.5}
en = {}
magmom = {}
for name in ['ZGNR12', 'ZGNR12_H', 'H2']:
    atoms = read(name+'.traj')
    en[name] = atoms.get_potential_energy()
    
e_ads = en['ZGNR12_H'] - en['ZGNR12'] - en['H2']
equal(e_ads ,reference['e_ads'], 0.001)

atoms =  read('ZGNR12_H.traj')
magmom_C = atoms.get_magnetic_moments()[1]
magmom_H = atoms.get_magnetic_moments()[25]
equal(magmom_C ,reference['magmom_C'], 0.001)
equal(magmom_H ,reference['magmom_H'], 0.001)
