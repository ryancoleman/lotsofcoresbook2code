import sys

from ase.calculators.jacapo.pw91_psp import defaultpseudopotentials as psp

from ase.test.tasks.dcdft import DeltaCodesDFTTask

from ase.tasks.calcfactory import calculator_factory

if len(sys.argv) == 1:
    xc = 'PBE'
    run = None
elif len(sys.argv) == 2:
    xc = 'PBE'
    run = sys.argv[1]
elif len(sys.argv) == 3:
    xc = sys.argv[1]
    run = sys.argv[2]

if run is not None:
    tag = 'scf_dcdft_%s_pw_%s' % (xc.lower(), run)
else:
    tag = 'scf_dcdft_%s_pw_jacapo' % (xc.lower())

if run is not None:
    tag += '_%s' % run

class Task(DeltaCodesDFTTask):
    def __init__(self, **kwargs):
        DeltaCodesDFTTask.__init__(self, **kwargs)

    def calculate(self, name, atoms):
        data = DeltaCodesDFTTask.calculate(self, name, atoms)
        try:
            steps = atoms.get_calculator().get_number_of_iterations()
        except (AttributeError, NotImplemented):
            steps = None
        data['calculator steps'] = steps
        return data

calcfactory = calculator_factory('jacapo',
                                 xc=xc,
                                 pw=340,
                                 dw=340,
                                 symmetry=True,
                                 ft=0.1,
                                 # need to be set explicitly
                                 spinpol=True,
                                 )

task = Task(
    calcfactory=calcfactory,
    tag=tag,
    use_lock_files=True,
    )

if __name__ == '__main__':
    keys = set(psp.keys()).intersection(set(task.collection.names))
    for m in ['Zn', 'Zr']:
        keys.remove(m)  # do not converge
    task.run(keys)
