# http://www.answers.com/topic/hartree-units

from math import pi

from ase.units import _hbar, _eps0, _me, _e, Hartree

_autime = _hbar**3*(4*pi*_eps0)**2/(_me*_e**4) # 1 autime ~ 2.42e-17 s

# Conversion factors between ASE and GPAW units

attosec_to_autime = 1e-18/_autime # 1 as ~ 0.0413 autime
autime_to_attosec = _autime*1e18 # 1 autime ~ 24.2 as

eV_to_hartree = 1.0/Hartree # 1 eV ~ 0.0368 Eh
hartree_to_eV = Hartree # 1 Eh ~ 27.2 eV

hartree_to_aufrequency = 1.0 # Eh ~ autime^(-1) since hbar=1
aufrequency_to_hartree = 1.0 # autime^(-1) ~ Eh since hbar=1

eV_to_aufrequency = 1.0/Hartree # 1 eV ~ 0.0368 autime^(-1)
aufrequency_to_eV = Hartree # 1 autime^(-1) ~ 27.2 eV

assert eV_to_aufrequency == eV_to_hartree * hartree_to_aufrequency
assert aufrequency_to_eV == aufrequency_to_hartree * hartree_to_eV

