from gpaw.atom.atompaw import AtomPAW
from gpaw.upf import UPFSetupData

def check(title, calc, epsref_n, threshold):
    eps_n = calc.wfs.kpt_u[0].eps_n
    print(title)
    for i, epsref in enumerate(epsref_n):
        err = abs(epsref - eps_n[i])
        ok = (err <= threshold)
        status = 'ok' if ok else 'FAILED'
        print('state %d | eps=%f | ref=%f | err=%f | tol=%s | %s' % (i,
                                                                     eps_n[i],
                                                                     epsref,
                                                                     err,
                                                                     threshold,
                                                                     status))
        assert ok

kwargs = dict(txt=None)

for setup in ['paw', 'hgh', UPFSetupData('H.pz-hgh.UPF')]:
    calc = AtomPAW('H', [[[1]]],
                   rcut=12.0, h=0.05,
                   setups={'H': setup}, **kwargs)
    tol = 5e-4 if setup in ['paw', 'hgh'] else 1e-3 # horrible UPF right now
    check('H %s' % setup, calc, [-0.233471], tol)


for setup in ['paw', 'hgh', UPFSetupData('O.pz-hgh.UPF')]:
    calc = AtomPAW('O', [[[2], [4]]],
                   rcut=10.0, h=0.035,
                   setups={'O': setup}, **kwargs)
    tol = 1e-3 if setup in ['paw', 'hgh'] else 5e-3 # horrible UPF right now
    check('O %s' % setup, calc, [-0.871362, -0.338381], tol)
