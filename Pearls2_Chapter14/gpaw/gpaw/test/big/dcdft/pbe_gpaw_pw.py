import sys

from gpaw import PW
from gpaw.mixer import Mixer, MixerDif
from gpaw.factory import GPAWFactory
from gpaw.occupations import FermiDirac, MethfesselPaxton

from ase.test.tasks.dcdft import DeltaCodesDFTTask as Task

from gpaw.atom.configurations import parameters

elements_slow = ['Na', 'Sb', 'Bi', 'Sc', 'Mn', 'F', 'B']

class Factory(GPAWFactory):
    def __init__(self, show_text_output=False, write_gpw_file=None,
                 **kwargs):
        GPAWFactory.__init__(self, show_text_output=show_text_output,
                             write_gpw_file=write_gpw_file,
                             **kwargs)

    def __call__(self, name, atoms):
        calculator = GPAWFactory.__call__(self, name, atoms)
        calculator.set(nbands=-2)
        if name.split('-')[0] in ['Ne', 'Ar', 'Kr', 'Xe', 'Rn',
                                  'K', 'Ca', 'Cs', 'Ba', 'Sr', 'Rb']:
            # MDTMP: not enough orbitals - use dzp!
            calculator.set(nbands=-1)
            calculator.set(basis='dzp')
        if name.split('-')[0] in ['Rb']:
            calculator.set(nbands=-1)
        if name.split('-')[0] in ['In', 'Cs', 'Sb', 'Ni', 'Ta', 'As',
                                  'V', 'Ni', 'Li', 'Rb', 'Tl']:
            calculator.set(nbands=-3)
        if name.split('-')[0] in ['Ca', 'Zn', 'Zr', 'Pb', 'Sr',
                                  'Xe', 'Rn', 'Ru', 'N', 'Os',
                                  'Cl', 'Cd']:
            calculator.set(nbands=-4)
        if name.split('-')[0] in ['Sc']:
            calculator.set(nbands=-5)
        if name.split('-')[0] in ['Ba', 'Os']:
            calculator.set(nbands=-6)
        if name.split('-')[0] in ['K']:
            calculator.set(nbands=-7)
        if name.split('-')[0] in ['Ne', 'Ar', 'Kr']:
            calculator.set(nbands=-8)
        if name.split('-')[0] in ['Cr']:
            calculator.set(nbands=-10)
        if name.split('-')[0] in ['Ti', 'Cr', 'Fe', 'Y', 'Nb', 'Mo',
                                  'Te', 'Hf', 'Re', 'Hg', 'Sb', 'Ca',
                                  'Pd', 'Ni', 'Ta',
                                  'Ru', 'Rh', 'V', 'Ag', 'Ir',
                                  'W', 'Pt', 'Mn', 'Sn',
                                  'Zn', 'Os', 'Mg', 'Zn',
                                  'Na', 'Bi', 'Sc', 'Zr',
                                  ]:
            calculator.set(eigensolver='cg')
            calculator.set(parallel={'band': 1})
        if name.split('-')[0] in ['Li', 'Na']:
            # https://listserv.fysik.dtu.dk/pipermail/gpaw-developers/2012-May/002870.html
            calculator.set(h=0.05)
        if name.split('-')[0] in ['Cs']:
            calculator.set(symmetry='off')
            calculator.set(eigensolver='cg')
            calculator.set(parallel={'band': 1})
        if name.split('-')[0] in ['Mn']:
            calculator.set(mixer=MixerDif())
            calculator.set(maxiter=450)
        if name.split('-')[0] in ['Cr']:
            calculator.set(mixer=MixerDif())
            calculator.set(maxiter=650)
        return calculator

if len(sys.argv) == 1:
    element = None
elif sys.argv[1] in elements_slow:
    element = sys.argv[1]
else:
    element = None

xc = 'PBE'

fit = (5, 0.02)

w = 0.06

ecut = 1200
kd = 8.0

tag = 'dcdft_%s_gpaw_pw' % xc.lower()

task = Task(
    calcfactory=Factory(xc=xc,
                        mode=PW(ecut),
                        occupations=FermiDirac(w),
                        maxiter=250,
                        kptdensity=kd,
                        ),
    tag=tag,
    fit=fit,
    use_lock_files=True,
    )

if __name__ == '__main__':
    if element is None:
        keys = set(parameters.keys()).intersection(set(task.collection.names))
        for s in elements_slow:
            keys.remove(s)  # those are slow, run separately
    else:
        keys = [element]
    task.run(keys)
