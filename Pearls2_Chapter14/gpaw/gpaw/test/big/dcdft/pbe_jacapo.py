from ase.calculators.jacapo.pw91_psp import defaultpseudopotentials as psp

from ase.tasks.calcfactory import calculator_factory

from ase.test.tasks.dcdft import DeltaCodesDFTTask as Task

xc = 'PBE'

fit = (5, 0.02)

w = 0.06

pw = 700
dw = 1000
kd = 8.0

tag = 'dcdft_%s_jacapo' % xc.lower()

calcfactory = calculator_factory('jacapo',
                                 pw=pw,
                                 dw=dw,
                                 xc='PBE',
                                 #convergence={
                                 #'energy':0.00001,
                                 #'density':0.000001, # 0.0001
                                 #'occupation':0.001, # 0.001
                                 #'maxsteps':None,
                                 #'maxtime':None
                                 #},
                                 ft=w,
                                 symmetry=True,
                                 spinpol=True, # must set explicitly
                                 calculate_stress=False,
                                 deletenc=True,  # start fresh every time
                                 kptdensity=kd,
                                 )

taskopts = {}

task = Task(
    calcfactory=calcfactory,
    tag=tag,
    fit=fit,
    use_lock_files=True,
    **taskopts)

if __name__ == '__main__':
    keys = set(psp.keys()).intersection(set(task.collection.names))
    for m in ['Be', 'Cr', 'Kr', 'Xe']:
        keys.remove(m)  # EOS fails: maybe need higher ecut?
    for m in ['Zn', 'Zr']:
        keys.remove(m)  # do not converge
    # run just only system to check if scripts work
    keys = ['Si']
    task.run(keys)
