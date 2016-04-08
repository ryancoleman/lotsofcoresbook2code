from math import sqrt, pi

import pylab as plt
import numpy as np

from ase.data import atomic_names as names


def analyse(generator, show=False):
    gen = generator
    symbol = gen.symbol

    colors = []
    id_j = []
    for j, n in enumerate(gen.vn_j):
        if n == -1:
            n = '*'
        id_j.append(str(n) + 'spdf'[gen.vl_j[j]])
        colors.append('krgbymc'[j])

    r_g = gen.r
    rmax = max(gen.rcut_l)

    # Construct logarithmic derivatives
    if len(gen.logd) > 0:
        rlog = gen.rlog
        logd = gen.logd
        elog = gen.elog
        ref = []
        for l, e in zip(gen.vl_j, gen.ve_j):
            i = elog.searchsorted(e)
            ref.append((elog[i], logd[l][1][i]))
        ref = np.array(ref)

    dpi = 80
    fig = plt.figure(figsize=(8.0, 12.0), dpi=dpi)
    fig.subplots_adjust(left=0.1, right=0.99, top=0.97, bottom=0.04,
                        hspace=0.25, wspace=0.26)

    plt.subplot(321)
    ymax = 0
    s = 1.0
    g = r_g.searchsorted(rmax * 2)
    for phi_g, phit_g, id, color in zip(gen.vu_j, gen.vs_j,
                                        id_j, colors):
        if id[0] != '*':
            lw = 2
            ymax = max(ymax, phi_g.max(), -phi_g.min())
        else:
            lw = 1
            s = ymax / max(np.abs(phi_g[:g]))
        plt.plot(r_g, s * phi_g, color + '-', label=id, lw=lw)
        plt.plot(r_g, s * phit_g, color + '--', label='_nolegend_', lw=lw)
    plt.legend(loc='lower right')
    plt.axis('tight')
    lim = plt.axis(xmin=0, xmax=rmax * 2, ymin=-ymax, ymax=ymax)
    plt.plot([rmax, rmax], lim[2:], 'k--', label='_nolegend_')
    plt.text(rmax, lim[2], r'$r_c$', ha='left', va='bottom', size=17)
    plt.title('Partial Waves')
    plt.xlabel('r [Bohr]')
    plt.ylabel(r'$r\ \phi,\ r\tilde{\phi}$, [Bohr$^{-1/2}$]')

    plt.subplot(322)
    ymax = 0
    s = 1.0
    for pt_g, id, color in zip(gen.vq_j, id_j, colors):
        if id[0] != '*':
            lw = 2
            ymax = max(ymax, pt_g.max(), -pt_g.min())
        else:
            lw = 1
            s = ymax / max(np.abs(pt_g))
        plt.plot(r_g, s * pt_g, color + '-', label=id, lw=lw)
    plt.axis('tight')
    lim = plt.axis(xmin=0, xmax=rmax * 1.2, ymin=-ymax, ymax=ymax)
    plt.plot([rmax, rmax], lim[2:], 'k--', label='_nolegend_')
    plt.text(rmax, lim[2], r'$r_c$', ha='left', va='bottom', size=17)
    plt.legend(loc='best')
    plt.title('Projectors')
    plt.xlabel('r [Bohr]')
    plt.ylabel(r'$r\ \tilde{p}$, [Bohr$^{-1/2}$]')

    plt.subplot(323)
    plt.plot(r_g, gen.nc, colors[0], label=r'$n_c$')
    plt.plot(r_g, gen.nct, colors[1], label=r'$\tilde{n}_c$')
    plt.axis('tight')
    lim = plt.axis(xmin=0, xmax=rmax * 1.2, ymax=max(gen.nct))
    plt.plot([rmax, rmax], lim[2:], 'k--', label='_nolegend_')
    plt.text(rmax, lim[2], r'$r_c$', ha='left', va='bottom', size=17)
    plt.legend(loc='best')
    plt.title('Densities')
    plt.xlabel('r [Bohr]')
    plt.ylabel(r'density [Bohr$^{-3}$]')

    plt.subplot(324)
    plt.plot(r_g[1:], gen.vr[1:] / r_g[1:] , label=r'$v$')
    plt.plot(r_g, gen.vt, label=r'$\tilde{v}$')
    plt.plot(r_g, gen.vt - gen.vbar, label=r'$\tilde{v}-\bar{v}$')
    plt.axis('tight')
    lim = plt.axis(xmin=0, xmax=rmax * 1.2,
                   ymin=(gen.vt - gen.vbar).min(), ymax=0)
    plt.plot([rmax, rmax], lim[2:], 'k--', label='_nolegend_')
    plt.text(rmax, lim[2], r'$r_c$', ha='left', va='bottom', size=17)
    plt.legend(loc='best')
    plt.title('Potentials')
    plt.xlabel('r [Bohr]')
    plt.ylabel('potential [Hartree]')

    plt.subplot(325)
    if len(gen.logd) > 0:
        plt.plot(ref[:, 0], ref[:, 1], 'ko', label='_nolegend_')
        for l, color in enumerate(colors[:4]):
            id = 'spdf'[l]
            plt.plot(elog, logd[l][0], linestyle='-', color=color, label=id)
            plt.plot(elog, logd[l][1], linestyle='--', color=color,
                     label='_nolegend_')
        plt.ylabel('log. deriv. at r=%.2f Bohr' % rlog)
        ymin = ref[:, 1].min()
        ymax = ref[:, 1].max()
        plt.axis(ymin=ymin - (ymax - ymin) * 0.1,
                 ymax=ymax + (ymax - ymin) * 0.1)
        plt.legend(loc='best')
    plt.title('Logarithmic Derivatives')
    plt.xlabel('Energy [Hartree]')

    plt.savefig('%s-setup.png' % symbol, dpi=dpi)
    
    if show:
        plt.show()
