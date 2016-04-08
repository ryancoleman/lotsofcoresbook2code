name= 'diamond333_hch'

a = 3.574

atoms = Atoms( [Atom('C', (0,0,0) ), Atom('C', (1,1,1)),
                Atom('C', (2,2,0) ), Atom('C', (3,3,1)),
                Atom('C', (2,0,2) ), Atom('C', (0,2,2)),
                Atom('C', (3,1,3) ), Atom('C', (1,3,3))],
               cell=(4,4,4),pbc=True)

atoms.set_cell((a,a,a), scale_atoms=True)
atoms *= (3,3,3)

calc = Calculator( h=0.2, txt = name +'.txt',
                    xc='PBE', setups={0:'hch1s'})

atoms.set_calculator(calc)

e = atoms.get_potential_energy()

calc.write(name + '.gpw')
