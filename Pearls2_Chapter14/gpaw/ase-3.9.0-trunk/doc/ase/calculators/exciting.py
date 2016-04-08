import os
from ase import Atoms
from ase.calculators.exciting import Exciting

# test structure, not real
a = Atoms('N3O', [(0, 0, 0), (1, 0, 0), (0, 0, 1), (0.5, 0.5, 0.5)], pbc=True)

calculator = Exciting(
    dir='excitingtestfiles',
    speciespath=os.environ['EXCITINGROOT']+'/species',
    paramdict={'title':{'text()':'N3O'},
               'groundstate':{'ngridk':'1 2 3','tforce':'true'},
               'relax':{},
               'properties':{'dos':{},
                             'bandstructure':
                                 {'plot1d':{'path':{'steps':'100',
                                                    'point':
                        [{'coord':'0.75000   0.50000   0.25000', 'label':'W'},
                         {'coord':'0.50000   0.50000   0.50000', 'label':'L'},
                         {'coord':'0.00000   0.00000   0.00000', 'label':'G'},  
                         {'coord':'0.50000   0.50000   0.00000', 'label':'X'},  
                         {'coord':'0.75000   0.50000   0.25000', 'label':'W'},  
                         {'coord':'0.75000   0.37500   0.37500', 'label':'K'}]
                                                    }}}}})

calculator.write(a)
