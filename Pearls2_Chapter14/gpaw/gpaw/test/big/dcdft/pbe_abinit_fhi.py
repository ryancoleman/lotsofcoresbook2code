from ase.tasks.calcfactory import CalculatorFactory

from ase.calculators.abinit import Abinit

from ase.test.tasks.dcdft import DeltaCodesDFTCollection
from ase.test.tasks.dcdft import DeltaCodesDFTTask as Task

xc = 'PBE'

fit = (5, 0.02)

w = 0.06

ecut = 1500
kd = 8.0

tag = 'dcdft_%s_abinit_fhi' % xc.lower()

class Factory(CalculatorFactory):
    def __init__(self, Class, name, label='label',
                 kpts=None, kptdensity=3.0,
                 **kwargs):
        CalculatorFactory.__init__(self, Class, name, label,
                                   kpts, kptdensity, **kwargs)
    def __call__(self, name, atoms):
        calculator = CalculatorFactory.__call__(self, name, atoms)
        # needed for cell optimization
        # http://www.abinit.org/documentation/helpfiles/for-v5.8/input_variables/varrlx.html#optcell
        # http://www.abinit.org/documentation/helpfiles/for-v5.8/input_variables/varrlx.html#ecutsm
        if 0:
            calculator.set({'ecutsm': 0.05}) # default 0.0
        # http://forum.abinit.org/viewtopic.php?f=8&t=1335
        if 0:
            calculator.set({'nsym': 1}) # default 0
        return calculator

calcopts = {'toldfe':1.0e-7,
            'ecut':ecut,
            'xc':'PBE',
            'pps':'fhi',
            'nstep':700,
            'diemac': 1000,
            'width':w,
            'kpts':kd,
            'chksymbreak':0,
            'ecutsm':0.05,
            'nsym':1,
}

calcfactory = Factory(Abinit, 'Abinit', 'label', **calcopts)

taskopts = {}

task = Task(
    calcfactory=calcfactory,
    tag=tag,
    fit=fit,
    use_lock_files=True,
    **taskopts)

if __name__ == '__main__':
    keys = task.collection.names
    for m in ['Mn']:
        keys.remove(m)  # EOS fails: reason unknown
    # run just only system to check if scripts work
    #task.run(keys)
    task.run(['Si'])
