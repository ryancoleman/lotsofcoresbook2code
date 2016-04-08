def agts(queue):
    co = queue.add('co.py', ncpus=4, walltime=5 * 60)
    queue.add('co.agts.py', deps=co)
 

if __name__ == '__main__':
    import numpy as np
    from ase.io import read

    configs = read('co.traj@:')
    e = [config.get_potential_energy() for config in configs]
    a = np.array([config.cell[0, 0] for config in configs])
    c = np.array([config.cell[2, 2] for config in configs])
    f = np.array([a**0, a, c, 0.5 * a**2, a * c, 0.5 * c**2])
    p = np.linalg.lstsq(f.T, e)[0]
    p0 = p[0]
    p1 = p[1:3]
    p2 = np.array([(p[3], p[4]),
                   (p[4], p[5])])
    a0, c0 = np.linalg.solve(p2.T, -p1)
    assert abs(a0 - 2.4937) < 0.001
    assert abs(c0 - 4.0424) < 0.001
    assert abs(a[4] - a0) < 0.001
    assert abs(c[4] - c0) < 0.001
