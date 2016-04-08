from ase.test import NotAvailable

from ase.calculators.mopac import Mopac

if Mopac().get_command() is None:
    raise NotAvailable('MOPAC required')

from ase.tasks.main import run

atoms, task = run("mopac molecule O2 O")
atoms, task = run('mopac molecule O2 O -s')
ae = 2 * task.data['O']['energy'] - task.data['O2']['energy']
print ae
assert abs(ae - 12.658) < 1e-3
