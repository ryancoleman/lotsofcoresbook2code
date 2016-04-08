import pickle
import matplotlib.pyplot as plt
point_names = ['W', 'L', '\Gamma', 'X', 'W', 'K']
x, X, e_kn = pickle.load(open('eigenvalues.pckl'))
e_kn -= e_kn[:, :4].max()
emin = e_kn.min() - 1
emax = 1
plt.figure(figsize=(6, 4))
for n in range(4):
    plt.plot(x, e_kn[:, n])
for p in X:
    plt.plot([p, p], [emin, emax], 'k-')
plt.xticks(X, ['$%s$' % n for n in point_names])
plt.axis(xmin=0, xmax=X[-1], ymin=emin, ymax=emax)
plt.xlabel('k-vector')
plt.ylabel('Energy (eV)')
plt.title('PBE bandstructure of Silicon')
plt.savefig('bs-PBE.png')
