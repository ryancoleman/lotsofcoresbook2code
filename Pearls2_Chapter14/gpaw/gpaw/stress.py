import numpy as np
import ase.units as units

from gpaw.utilities import unpack
from gpaw.wavefunctions.pw import PWWaveFunctions

def stress(calc):
    wfs = calc.wfs
    dens = calc.density
    ham = calc.hamiltonian

    if not isinstance(wfs, PWWaveFunctions):
        raise NotImplementedError('Calculation of stress tensor is only ' +
                                  'implemented for plane-wave mode.')
    if ham.xc.orbital_dependent:
        raise NotImplementedError('Calculation of stress tensor is not ' +
                                  'implemented for orbital-dependent ' +
                                  'XC functionals such as '+ham.xc.name)

    calc.timer.start('Stress tensor')

    s_vv = wfs.get_kinetic_stress().real

    s_vv += ham.xc.stress_tensor_contribution(dens.nt_sg)

    pd = dens.pd3
    p_G = 4 * np.pi * dens.rhot_q
    p_G[1:] /= pd.G2_qG[0][1:]**2
    G_Gv = pd.get_reciprocal_vectors()
    for v1 in range(3):
        s_vv[v1,v1] -= ham.epot
        for v2 in range(3):
            s_vv[v1,v2] += pd.integrate(p_G, dens.rhot_q *
                                        G_Gv[:, v1] * G_Gv[:, v2])
    s_vv += dens.ghat.stress_tensor_contribution(ham.vHt_q, dens.Q_aL)

    s_vv -= np.eye(3) * ham.ebar
    s_vv += ham.vbar.stress_tensor_contribution(dens.nt_sQ.sum(0))

    s_vv += dens.nct.stress_tensor_contribution(ham.vt_Q)

    s0 = 0.0
    s0_vv = 0.0
    for kpt in wfs.kpt_u:
        a_ani = {}
        for a, P_ni in kpt.P_ani.items():
            Pf_ni = P_ni * kpt.f_n[:, None]
            dH_ii = unpack(ham.dH_asp[a][kpt.s])
            dS_ii = ham.setups[a].dO_ii
            a_ni = (np.dot(Pf_ni, dH_ii) -
                    np.dot(Pf_ni * kpt.eps_n[:, None], dS_ii))
            s0 += np.vdot(P_ni, a_ni)
            a_ani[a] = 2 * a_ni.conj()
        s0_vv += wfs.pt.stress_tensor_contribution(kpt.psit_nG, a_ani,
                                               q=kpt.q)
    s0_vv -= s0.real * np.eye(3)
    wfs.bd.comm.sum(s0_vv)
    wfs.kd.comm.sum(s0_vv)
    s_vv += s0_vv

    vol = calc.atoms.get_volume() / units.Bohr**3
    s_vv = 0.5 / vol * (s_vv + s_vv.T)

    # Symmetrize:
    sigma_vv = np.zeros((3, 3))
    cell_cv = wfs.gd.cell_cv
    for U_cc in wfs.kd.symmetry.op_scc:
        M_vv = np.dot(np.linalg.inv(cell_cv),
                      np.dot(U_cc, cell_cv)).T
        sigma_vv += np.dot(np.dot(M_vv.T, s_vv), M_vv)
    sigma_vv /= len(wfs.kd.symmetry.op_scc)
    
    # Make sure all agree on the result (redundant calculation on
    # different cores involving BLAS might give slightly different
    # results):
    wfs.world.broadcast(sigma_vv, 0)

    calc.text('Stress tensor:')
    for sigma_v in sigma_vv:
        calc.text('%12.6f %12.6f %12.6f' %
                  tuple(units.Hartree / units.Bohr**3 * sigma_v))

    calc.timer.stop('Stress tensor')

    return sigma_vv
