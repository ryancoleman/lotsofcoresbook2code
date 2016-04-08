def agts(queue):
    al = queue.add('al.py', ncpus=8, walltime=12 * 60)
    queue.add('al.agts.py', deps=[al],
              creates=['Al_conv_ecut.png', 'Al_conv_k.png'])

if __name__ == '__main__':
    import pylab as plt
    from ase.utils.eos import EquationOfState
    from ase.io import read

    def fit(filename):
        configs = read(filename + '@:')
        volumes = [a.get_volume() for a in configs]
        energies = [a.get_potential_energy() for a in configs]
        eos = EquationOfState(volumes, energies)
        v0, e0, B = eos.fit()
        return (4 * v0)**(1 / 3.0)

    cutoffs = range(200, 501, 50)
    a = [fit('Al-%d.txt' % ecut) for ecut in cutoffs]
    plt.figure(figsize=(6, 4))
    plt.plot(cutoffs, a, 'o-')
    plt.axis(ymin=4.03, ymax=4.05)
    plt.xlabel('Plane-wave cutoff energy [eV]')
    plt.ylabel('lattice constant [Ang]')
    plt.savefig('Al_conv_ecut.png')

    kpoints = range(4, 17)
    plt.figure(figsize=(6, 4))
    a = [fit('Al-%02d.txt' % k) for k in kpoints]
    plt.plot(kpoints, a, '-')
    plt.xlabel('number of k-points')
    plt.ylabel('lattice constant [Ang]')
    plt.savefig('Al_conv_k.png')
    plt.show()
