from gpaw.atom.generator2 import _generate


def derivatives(atoms):
    """Calculate derivatives of the energy with respect to PAW-setup
    parameters."""

    e00 = atoms.get_potential_energy()

    allsetups = atoms.calc.wfs.setups.setups

    parameters = {}
    for setup in allsetups.values():
        d = setup.data
        projectors = []
        radii = []
        ll = -1
        for n, l, e, r in zip(d.n_j, d.l_j, d.eps_j, d.rcut_j):
            if l > ll:
                radii.append(r)
                ll = l
            if n > 0:
                e = n
            projectors.append('%r%s' % (e, 'spdfg'[l]))

        kwargs = dict(symbol=d.symbol,
                      xc=d.xcname,
                      projectors=','.join(projectors),
                      radii=radii,
                      scalar_relativistic=(d.type == 'scalar-relativistic'),
                      alpha=1 / d.rcgauss**2,
                      gamma=d.gamma,
                      h=d.h,
                      l0=d.l0,
                      r0=d.r0,
                      nderiv0=d.nderiv0,
                      e0=d.e0)
        parameters[d.symbol] = kwargs

    
    for symbol, p in parameters.items():
        if atoms.calc.wfs.world.rank == 0:
            gen = _generate(**p)
            gen.make_paw_setup('derivative0').write_xml()
        atoms.calc.wfs.world.barrier()
        
    atoms.calc.set(setups='derivative0')
    e0 = atoms.get_potential_energy()
    
    def f(p, name):
        if atoms.calc.wfs.world.rank == 0:
            gen = _generate(**p)
            gen.make_paw_setup(name).write_xml()
        atoms.calc.wfs.world.barrier()
        atoms.calc.set(setups={None: 'derivative0', p['symbol']: name})
        e = atoms.get_potential_energy()
        return e

    results = {}
    for symbol, p in parameters.items():
        derivs = []
        for i, r in enumerate(p['radii']):
            p['radii'][i] = 1.01 * r
            e = f(p, 'derivative_rc%d' % i)
            derivs.append((e - e0) / (0.01 * r))
            p['radii'][i] = r

        alpha = p['alpha']
        p['alpha'] = 1.01 * alpha
        e = f(p, 'derivative_alpha')
        derivs.append((e - e0) / (0.01 * alpha))
        p['alpha'] = alpha

        r0 = p['r0']
        p['r0'] = 1.01 * r0
        e = f(p, 'derivative_r0')
        derivs.append((e - e0) / (0.01 * r0))
        p['r0'] = r0

        results[symbol] = derivs

    return results

