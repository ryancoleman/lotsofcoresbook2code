from __future__ import print_function
from ase import Atom, Atoms
from ase.optimize.bfgslinesearch import BFGSLineSearch
from ase.io import read
from gpaw import GPAW

#Relax the Au 55 atoms cluster with both lcao and paw modes
cluster = read('Au55.xyz')
cell = [(17.79365715,0,0),
        (0,19.60846479,0),
        (0,0, 19.84025464)]
cluster.set_cell(cell,scale_atoms=False)
cluster.center()

kwargs_lcao = dict(mode='lcao',
              #basis='dzp',
              convergence={'density':0.1, 'energy':0.1}
              )#poissonsolver=poissonsolver)

calc = GPAW(h=0.18,txt=None,**kwargs_lcao)
cluster.set_calculator(calc)

dyn1 = BFGSLineSearch(cluster, trajectory='Au_cluster_lcao.traj')
dyn1.run(fmax=0.02)
e_cluster_lcao = cluster.get_potential_energy()

#cluster.set_calculator(GPAW(h=0.18,tex=None))
dyn2 = BFGSLineSearch(cluster, trajectory='Au_cluster_paw.traj')
dyn2.run(fmax=0.02)
e_cluster_paw = cluster.get_potential_energy()

#Relax CO molecule with both lcao and paw modes
CO = Atoms([Atom('C',(1.0,1.0,1.0)),
            Atom('O',(1.0,1.0,2.3))],
            cell=(12,12.5,14.5))
CO.set_calculator(calc)
CO.center()
dyn3 = BFGSLineSearch(CO)
dyn3.run(fmax=0.02)
e_CO_lcao = CO.get_potential_energy()

CO.set_calculator(GPAW(h=0.18,tex=None))
dyn4 = BFGSLineSearch(CO)
dyn4.run(fmax=0.02)
e_CO_paw = CO.get_potential_energy()

CO_bond = CO.get_positions()[1][2] - CO.get_positions()[0][2] 

#Attach CO molecule onto the Au cluster
pos=[]
pos.append(cluster[34].position.copy())
pos.append(cluster[34].position.copy())
pos[0][2] += 1.8
pos[1][2] += 1.8 + CO_bond
CO = Atoms([Atom('C', pos[0]),
            Atom('O', pos[1])])
cluster.extend(CO)

#Relax the CO adsorbed Au cluster with both Lcao and paw modes
cluster.set_calculator(calc)
dyn5 = BFGSLineSearch(cluster, trajectory='CO_cluster_lcao.traj')
dyn5.run(fmax=0.02)
e_cocluster_lcao = cluster.get_potential_energy()

cluster.set_calculator(GPAW(h=0.18,tex=None))
dyn6 = BFGSLineSearch(cluster, trajectory='CO_cluster_paw.traj')
dyn6.run(fmax=0.02)
e_cocluster_paw = cluster.get_potential_energy()

#Print results
print('Adsorption energy of CO on Au cluster (lcao):', e_cocluster_lcao - e_CO_lcao - e_cluster_lcao)
print('Adsorption energy of CO on Au cluster (paw):', e_cocluster_paw - e_CO_paw - e_cluster_paw)
