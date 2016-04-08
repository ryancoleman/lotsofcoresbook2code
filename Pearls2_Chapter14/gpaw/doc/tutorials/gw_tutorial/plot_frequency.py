import pickle
# mathtext fails to create labels with matplotlib 0.99 on el6
import matplotlib
matplotlib.rc('text', usetex=True)
import matplotlib.pyplot as plt
from ase.units import Hartree

plt.figure(1)
plt.figure(figsize=(6.5, 4.5))
plt.figure(2)
plt.figure(figsize=(6.5, 4.5))

for wlin in [25.,50.,75.,100.]:

    e=[]
    GW_gap=[]

    for dw in [0.02,0.05,0.1,0.2,0.5]:

        data=pickle.load(open('Si_GW_wlin%s_dw%s.pckl' % (wlin, dw)))
        QP_skn=data['QP_skn']
        QPgap = (QP_skn[0,0,2] - QP_skn[0,0,1])*Hartree

        GW_gap.append(QPgap)
        e.append(dw)

    plt.plot(e, GW_gap, 'o-', label='$\omega_{\mathrm{lin}} = $ %s eV' %wlin)

plt.xlabel('$\Delta \omega$ (eV)')
plt.ylabel('Direct band gap (eV)')
plt.xlim([0.015, 0.6])
plt.xscale('log')
ax = plt.axes()
ax.set_xticks((0.02, 0.05, 0.1, 0.2, 0.5))
ax.set_xticklabels(('0.02', '0.05', '0.1', '0.2', '0.5'))
plt.ylim([3.2, 3.9])
plt.title('G$_0$W$_0$@LDA')
plt.legend(loc='upper left')
plt.savefig('Si_w.png')
