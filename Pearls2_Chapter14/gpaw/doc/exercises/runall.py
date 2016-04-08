import os
import sys
import pylab
n = 1
def show():
    global n
    pylab.savefig('x%d.png' % n)
    n += 1
pylab.show = show

if 0:
    def view(*args, **kwargs):
        return
    import ase
    ase.view = view

os.system('rm -rf test; mkdir test')
os.chdir('test')
for dir, script in [
    ('neb', 'neb1.py'),
    ('aluminium', 'Al_fcc.py'),
    ('aluminium', 'Al_fcc_convergence.py'),
    ('surface', 'Al100.py'),
    ('surface', 'work_function.py'),
    ('diffusion', 'initial.py'),
    ('diffusion', 'solution.py'),
    ('diffusion', 'densitydiff.py'),
    ('vibrations', 'h2o.py'),
    ('vibrations', 'H2O_vib.py'),
    ('iron', 'ferro.py'),
    ('iron', 'anti.py'),
    ('iron', 'non.py'),
    ('iron', 'PBE.py'),
    ('band_structure', 'Na_band.py'),
    ('band_structure', 'plot_band.py'),
    ('wannier', 'si.py'),
    ('wannier', 'wannier-si.py'),
    ('wannier', 'benzene.py'),
    ('wannier', 'wannier-benzene.py'),
    ('stm', 'HAl100.py'),
    ('wavefunctions', 'CO.py'),
    ('dos', 'pdos.py'),
    ('lrtddft', 'ground_state.py'),
    ('transport', 'pt_h2_tb_transport.py'),
    ('transport', 'pt_h2_lcao.py'),
    ('transport', 'pt_h2_lcao_transport.py')
    ]:
    execfile('../' + dir + '/' + script, {'k': 6, 'N': 5})
for dir, script, args in [
    ('stm', 'stm.py', ['HAl100.gpw']),
    ('dos', 'dos.py', ['Al-fcc.gpw', 'si.gpw', 'CO.gpw',
                       'ferro.gpw', 'anti.gpw', 'non.gpw'])]:
    for arg in args:
        sys.argv = ['', arg]
        execfile('../' + dir + '/' + script)
