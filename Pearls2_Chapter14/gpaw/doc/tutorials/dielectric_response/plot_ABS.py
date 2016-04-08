import numpy as np
import matplotlib.pyplot as plt
    
plt.figure(figsize=(7, 5))
d = np.loadtxt('si_abs.csv', delimiter=',')
plt.plot(d[:, 0], d[:, 3], '-k', label='$\mathrm{Re}\epsilon(\omega)$')
plt.plot(d[:, 0], d[:, 4], '-r', label='$\mathrm{Im}\epsilon(\omega)$')

plt.title('Dielectric function of Si')
plt.legend()
plt.xlabel('Energy (eV)', fontsize=14)
plt.ylabel('$\epsilon$', fontsize=18)

ax = plt.gca()


def y(e):
    x = d[:, 0]
    i = (x < e).sum()
    return d[i, 4]

    
# data from G.Kresse, PRB 73, 045112 (2006)
for name, e in zip("E_0 E_1 E_2 E_0' E_1'".split(),
                   [2.53, 2.71, 3.72, 3.08, 4.50]):
    arr = plt.arrow(e, y(e) + 4.0, 0, -3,
                    width=0.01, head_width=0.1, head_length=1)
    ax.add_patch(arr)
    plt.text(e, y(e) + 4, '$' + name + '$')

plt.xlim(0, 6)
plt.ylim(0, 60)
plt.savefig('silicon_ABS.png', bbox_inches='tight')
plt.show()
