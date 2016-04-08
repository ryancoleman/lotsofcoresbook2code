from __future__ import print_function

import numpy as np

try:
    # Matplotlib is not a dependency
    import matplotlib as mpl
    mpl.use('Agg')  # force the antigrain backend
    from matplotlib.figure import Figure
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    import os
    import warnings
    # silence matplotlib.use() warning
    warnings.filterwarnings('ignore', '.*This call to matplotlib\.use.*',)
except (ImportError, RuntimeError):
    mpl = None

from scipy.optimize import leastsq
from ase.units import second
from ase.io.trajectory import PickleTrajectory

# Dimer oscillation model used for least-squares fit 
def f(p, t):
    return p[0] * np.cos(p[1] * (t - p[2])) + p[3]

# Jacobian of model with respect to its four parameters
def Df(p, t):
    return np.array([
        np.cos(p[1] * (t - p[2])),
        -p[0] * np.sin(p[1] * (t - p[2])) * (t - p[2]),
        p[0] * np.sin(p[1] * (t - p[2])) * p[1],
        np.ones_like(t)])


for name in ['h2_osc', 'n2_osc', 'na2_md', 'na2_osc']:
    print('\nAnalysing %s\n%s' % (name, '-'*32))

    # Import relevant test and make sure it has the prerequisite parameters
    m = __import__(name, {}, {})
    for attr in ['d_bond', 'd_disp', 'timestep', 'period', 'ndiv', 'niter']:
        if not hasattr(m, attr):
            raise ImportError('Module %s has no %s value' % (name, attr))

    # Read dimer bond length time series from trajectory file
    traj = PickleTrajectory(name + '_td.traj', 'r')
    nframes = len(traj)
    natoms = len(traj[0])
    symbol = traj[0].get_name()
    t_i = m.timestep * m.ndiv * np.arange(nframes)
    Ekin_i, Epot_i = np.empty(nframes), np.empty(nframes)
    R_iav = np.empty((nframes, natoms, 3))
    V_iav = np.empty((nframes, natoms, 3))
    A_iav = np.empty((nframes, natoms, 3))
    for i in range(nframes):
        Ekin_i[i] = traj[i].get_kinetic_energy()
        Epot_i[i] = traj[i].get_potential_energy()
        R_iav[i] = traj[i].get_positions()
        V_iav[i] = traj[i].get_velocities()
        A_iav[i] = traj[i].get_forces() / traj[i].get_masses()[:,np.newaxis]
    print('Read %d frames from trajectory...' % nframes)
    assert nframes * m.ndiv == m.niter, (nframes, m.ndiv, m.niter)
    traj.close()

    # Verify that energy was conserved
    dEstd = np.std(Ekin_i + Epot_i)
    dEmax = np.max(Ekin_i.ptp(), Epot_i.ptp())
    assert dEstd < 1e-2 * dEmax + 1e-6, (dEstd, dEmax)

    # Compare position, velocity and force time series using Velocity Verlet
    dt = m.timestep * m.ndiv * 1e-18 * second
    Rn_iav = R_iav[:-1] + V_iav[:-1] * dt + 0.5 * A_iav[:-1] * dt**2
    Vn_iav = V_iav[:-1] + 0.5 * (A_iav[:-1] + A_iav[1:]) * dt
    dRstd_av = (np.sum((Rn_iav - R_iav[1:])**2, axis=0) / len(Rn_iav))**0.5
    dVstd_av = (np.sum((Vn_iav - V_iav[1:])**2, axis=0) / len(Vn_iav))**0.5
    dRmax_av, dVmax_av = R_iav.ptp(axis=0), V_iav.ptp(axis=0)
    assert np.all(dRstd_av < 1e-4 * dRmax_av + 1e-12), (dRstd_av, dRmax_av)
    assert np.all(dVstd_av < 1e-4 * dVmax_av + 1e-12), (dVstd_av, dVmax_av)

    # Fit model to time series using imported parameters as an initial guess
    d_i = np.sum((R_iav[:,1] - R_iav[:,0])**2, axis=-1)**0.5
    p0 = (m.d_disp, 2 * np.pi / m.period, -m.timestep * m.ndiv, m.d_bond)
    p, cov, info, msg, status = leastsq(lambda p: f(p, t_i) - d_i, p0, \
        Dfun=lambda p: Df(p, t_i), col_deriv=True, full_output=True)
    print('leastsq returned %d: %s' % (status, msg.replace('\n ','')))
    print('p0=', np.asarray(p0))
    print('p =', p)
    assert status in range(1,4+1), (p, cov, info, msg, status)

    tol = 0.1 #TODO use m.reltol
    err = np.abs(2 * np.pi / p[1] - m.period) / m.period
    print('T=%13.9f fs, Tref=%13.9f fs, err=%5.2f %%, tol=%.1f %%' \
        % (2 * np.pi / p[1] * 1e-3, m.period * 1e-3, 1e2 * err, 1e2 * tol))

    if mpl:
        fig = Figure()
        ax = fig.add_axes([0.11, 0.1, 0.86, 0.83])
        raw = r',\;'.join([r'T=%.2f\mathrm{\,fs}',
                           r'T_\mathrm{ref}=%.2f\mathrm{\,fs}',
                           r'\eta=%.2f\,\%%'])
        mathmode = raw % (2 * np.pi / p[1] * 1e-3, m.period * 1e-3, 1e2 * err)
        ax.set_title(symbol + ' ($' + mathmode + '$)')
        ax.set_xlabel('Time [fs]')
        ax.set_ylabel('Dimer bond length [Ang]')
        ax.plot(t_i * 1e-3, d_i, '-b', t_i * 1e-3, f(p, t_i), '--k')
        ax.legend(('Ehrenfest data', r'$A\,\mathrm{cos}(\omega(t-t_0))+B$'))
        FigureCanvasAgg(fig).print_figure(name + '.png', dpi=90)

    if err > tol:
        print('Relative error %f %% > tolerance %f %%' % (1e2 * err, 1e2 * tol))
        raise SystemExit(1)

