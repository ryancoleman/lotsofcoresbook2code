import os
import gc
import sys
import time
import signal
import traceback

import numpy as np

from ase.utils import devnull

from gpaw.atom.generator import Generator
from gpaw.atom.configurations import parameters, tf_parameters
from gpaw.utilities import compiled_with_sl
from gpaw import setup_paths
from gpaw import mpi
import gpaw


def equal(x, y, tolerance=0, fail=True, msg=''):
    """Compare x and y."""

    if not np.isfinite(x - y).any() or (np.abs(x - y) > tolerance).any():
        msg = (msg + '%s != %s (error: |%s| > %.9g)' %
               (x, y, x - y, tolerance))
        if fail:
            raise AssertionError(msg)
        else:
            sys.stderr.write('WARNING: %s\n' % msg)


def findpeak(x, y):
    dx = x[1] - x[0]
    i = y.argmax()
    a, b, c = np.polyfit([-1, 0, 1], y[i - 1:i + 2], 2)
    assert a < 0
    x = -0.5 * b / a
    return dx * (i + x), a * x**2 + b * x + c

    
def gen(symbol, exx=False, name=None, **kwargs):
    if mpi.rank == 0:
        if 'scalarrel' not in kwargs:
            kwargs['scalarrel'] = True
        g = Generator(symbol, **kwargs)
        if 'orbital_free' in kwargs:
            g.run(exx=exx, name=name, use_restart_file=False,
                  **tf_parameters.get(symbol, {'rcut': 0.9}))
        else:
            g.run(exx=exx, name=name, use_restart_file=False,
                  **parameters[symbol])
    mpi.world.barrier()
    if setup_paths[0] != '.':
        setup_paths.insert(0, '.')


def wrap_pylab(names=[]):
    """Use Agg backend and prevent windows from popping up."""
    import matplotlib
    matplotlib.use('Agg')
    import pylab

    def show(names=names):
        if names:
            name = names.pop(0)
        else:
            name = 'fig.png'
        pylab.savefig(name)

    pylab.show = show


tests = [
    'gemm_complex.py',
    'ase3k_version.py',
    'kpt.py',
    'mpicomm.py',
    'numpy_core_multiarray_dot.py',
    'maxrss.py',  # verifies reported RAM allocation: fragile, don't move down
    'fileio/hdf5_noncontiguous.py',
    'cg2.py',
    'laplace.py',
    'lapack.py',
    'eigh.py',
    'parallel/submatrix_redist.py',
    'second_derivative.py',
    'parallel/parallel_eigh.py',
    'gp2.py',
    'blas.py',
    'Gauss.py',
    'nabla.py',
    'dot.py',
    'mmm.py',
    'lxc_fxc.py',
    'pbe_pw91.py',
    'gradient.py',
    'erf.py',
    'lf.py',
    'fsbt.py',
    'parallel/compare.py',
    'integral4.py',
    'zher.py',
    'gd.py',
    'pw/interpol.py',
    'screened_poisson.py',
    'xc.py',
    'XC2.py',
    'yukawa_radial.py',
    'dump_chi0.py',
    'vdw/potential.py',
    'lebedev.py',
    'fileio/hdf5_simple.py',
    'occupations.py',
    'derivatives.py',
    'parallel/realspace_blacs.py',
    'pw/reallfc.py',
    'parallel/pblas.py',
    'non_periodic.py',
    'spectrum.py',
    'pw/lfc.py',
    'gauss_func.py',
    'multipoletest.py',
    'noncollinear/xcgrid3d.py',
    'cluster.py',
    'poisson.py',
    'parallel/overlap.py',
    'parallel/scalapack.py',
    'gauss_wave.py',
    'transformations.py',
    'parallel/blacsdist.py',
    'ut_rsh.py',
    'pbc.py',
    'noncollinear/xccorr.py',
    'atoms_too_close.py',
    'harmonic.py',
    'proton.py',
    'timing.py',                            # ~1s
    'parallel/ut_parallel.py',              # ~1s
    'ut_csh.py',                            # ~1s
    'lcao_density.py',                      # ~1s
    'parallel/hamiltonian.py',              # ~1s
    'pw/stresstest.py',                     # ~1s
    'pw/fftmixer.py',                       # ~1s
    'usesymm.py',                           # ~1s
    'coulomb.py',                           # ~1s
    'xcatom.py',                            # ~1s
    'force_as_stop.py',                     # ~1s
    'vdwradii.py',                          # ~1s
    'ase3k.py',                             # ~1s
    'numpy_zdotc_graphite.py',              # ~1s
    'eed.py',                               # ~1s
    'gemv.py',                              # ~2s
    'fileio/idiotproof_setup.py',           # ~2s
    'ylexpand.py',                          # ~2s
    'keep_htpsit.py',                       # ~2s
    'gga_atom.py',                          # ~2s
    'hydrogen.py',                          # ~2s
    'restart2.py',                          # ~2s
    'aeatom.py',                            # ~2s
    'plt.py',                               # ~2s
    'ds_beta.py',                           # ~2s
    'multipoleH2O.py',                      # ~2s
    'noncollinear/h.py',                    # ~2s
    'stdout.py',                            # ~2s
    'lcao_largecellforce.py',               # ~2s
    'parallel/scalapack_diag_simple.py',    # ~2s
    'fixdensity.py',                        # ~2s
    'pseudopotential/ah.py',                # ~2s
    'lcao_restart.py',                      # ~2s
    'wfs_io.py',                            # ~3s
    'lrtddft2.py',                          # ~3s
    'fileio/file_reference.py',             # ~3s
    'cmrtest/cmr_test2.py',                 # ~3s
    'restart.py',                           # ~3s
    'broydenmixer.py',                      # ~3s
    'pw/fulldiagk.py',                      # ~3s
    'external_potential.py',                # ~3s
    'mixer.py',                             # ~3s
    'parallel/lcao_projections.py',         # ~3s
    'lcao_h2o.py',                          # ~3s
    'h2o_xas.py',                           # ~3s
    'wfs_auto.py',                          # ~3s
    'pw/fulldiag.py',                       # ~3s
    'symmetry_ft.py',                       # ~3s
    'aluminum_EELS_RPA.py',                 # ~3s
    'ewald.py',                             # ~4s
    'symmetry.py',                          # ~4s
    'revPBE.py',                            # ~4s
    'tf_mode_pbc.py',                       # ~4s
    'tf_mode.py',                           # ~4s
    'nonselfconsistentLDA.py',              # ~4s
    'aluminum_EELS_ALDA.py',                # ~4s
    'spin_contamination.py',                # ~4s
    'inducedfield_lrtddft.py',              # ~4s
    'H_force.py',                           # ~4s
    'usesymm2.py',                          # ~4s
    'mgga_restart.py',                      # ~4s
    'fixocc.py',                            # ~4s
    'spinFe3plus.py',                       # ~4s
    'fermisplit.py',                        # ~4s
    'Cl_minus.py',                          # ~4s
    'h2o_xas_recursion.py',                 # ~5s
    'nonselfconsistent.py',                 # ~5s
    'spinpol.py',                           # ~5s
    'exx_acdf.py',                          # ~5s
    'cg.py',                                # ~5s
    'kptpar.py',                            # ~5s
    'elf.py',                               # ~5s
    'blocked_rmm_diis.py',                  # ~5s
    'pw/slab.py',                           # ~5s
    'si.py',                                # ~5s
    'lcao_bsse.py',                         # ~5s
    'parallel/lcao_hamiltonian.py',         # ~5s
    'degeneracy.py',                        # ~5s
    'refine.py',                            # ~5s
    'gemm.py',                              # ~6s
    'al_chain.py',                          # ~6s
    'fileio/parallel.py',                   # ~6s
    'fixmom.py',                            # ~6s
    'exx_unocc.py',                         # ~6s
    'davidson.py',                          # ~6s
    'aedensity.py',                         # ~7s
    'pw/h.py',                              # ~7s
    'apmb.py',                              # ~7s
    'pseudopotential/hgh_h2o.py',           # ~7s
    'ed_wrapper.py',                        # ~7s
    'pw/bulk.py',                           # ~7s
    'ne_gllb.py',                           # ~7s
    'ed.py',                                # ~7s
    'lcao_force.py',                        # ~7s
    'fileio/restart_density.py',            # ~8s
    'rpa_energy_Ni.py',                     # ~8s
    'be_nltd_ip.py',                        # ~8s
    'test_ibzqpt.py',                       # ~8s
    'si_primitive.py',                      # ~9s
    'inducedfield_td.py',                   # ~9s
    'ehrenfest_nacl.py',                    # ~9s
    'fd2lcao_restart.py',                   # ~9s
    'gw_method.py',                         # ~9s
    'constant_electric_field.py',           # ~9s
    'complex.py',                           # ~9s
    'vdw/quick.py',                         # ~9s
    'bse_aluminum.py',                      # ~10s
    'Al2_lrtddft.py',                       # ~10s
    'ralda_energy_N2.py',                   # ~10s
    'gw_ppa.py',                            # ~10s
    'parallel/lcao_complicated.py',         # ~10s
    'bulk.py',                              # ~10s
    'scfsic_h2.py',                         # ~10s
    'lcao_bulk.py',                         # ~11s
    '2Al.py',                               # ~11s
    'kssingles_Be.py',                      # ~11s
    'relax.py',                             # ~11s
    'pw/mgo_hybrids.py',                    # ~11s
    'dscf_lcao.py',                         # ~12s
    '8Si.py',                               # ~12s
    'partitioning.py',                      # ~12s
    'lxc_xcatom.py',                        # ~12s
    'gllbatomic.py',                        # ~13s
    'guc_force.py',                         # ~13s
    'ralda_energy_Ni.py',                   # ~13s
    'simple_stm.py',                        # ~13s
    'ed_shapes.py',                         # ~14s
    'restart_band_structure.py',            # ~14s
    'exx.py',                               # ~14s
    'Hubbard_U.py',                         # ~15s
    'rpa_energy_Si.py',                     # ~15s
    'dipole.py',                            # ~15s
    'IP_oxygen.py',                         # ~15s
    'rpa_energy_Na.py',                     # ~15s
    'parallel/fd_parallel.py',              # ~15s
    'parallel/lcao_parallel.py',            # ~16s
    'atomize.py',                           # ~16s
    'excited_state.py',                     # ~16s
    'ne_disc.py',                           # ~16s
    'tpss.py',                              # ~18s
    'td_na2.py',                            # ~18s
    'exx_coarse.py',                        # ~18s
    'pplda.py',                             # ~18s
    'si_xas.py',                            # ~18s
    'mgga_sc.py',                           # ~19s
    'Hubbard_U_Zn.py',                      # ~20s
    # buildbot > 20 sec tests start here (add tests after lrtddft.py!)
    'lrtddft.py',                           # ~20s
    'parallel/fd_parallel_kpt.py',          # ~21s
    'pw/hyb.py',                            # ~21s
    'Cu.py',                                # ~21s
    'response_na_plasmon.py',               # ~22s
    'bse_diamond.py',                       # ~23s
    'fermilevel.py',                        # ~23s
    'parallel/ut_hsblacs.py',               # ~23s
    'ralda_energy_H2.py',                   # ~23s
    'diamond_absorption.py',                # ~24s
    'ralda_energy_Si.py',                   # ~24s
    'ldos.py',                              # ~25s
    'revPBE_Li.py',                         # ~26s
    'parallel/lcao_parallel_kpt.py',        # ~29s
    'h2o_dks.py',                           # ~30s
    'nsc_MGGA.py',                          # ~32s
    'diamond_gllb.py',                      # ~33s
    'MgO_exx_fd_vs_pw.py',                  # ~37s
    'vdw/quick_spin.py',                    # ~37s
    'bse_sym.py',                           # ~40s
    'parallel/ut_hsops.py',                 # ~41s
    'LDA_unstable.py',                      # ~42s
    'au02_absorption.py',                   # ~44s
    'wannierk.py',                          # ~45s
    'bse_vs_lrtddft.py',                    # ~45s
    'aluminum_testcell.py',                 # ~46s
    'pygga.py',                             # ~47s
    'ut_tddft.py',                          # ~49s
    'rpa_energy_N2.py',                     # ~52s
    'vdw/ar2.py',                           # ~53s
    'parallel/diamond_gllb.py',             # ~59s
    'beef.py',
    'pw/si_stress.py',                      # ~61s
    'chi0.py',                              # ~71s
    'scfsic_n2.py',                         # ~73s
    'transport.py',                         # ~73s
    'lrtddft3.py',                          # ~75s
    'nonlocalset.py',                       # ~82s
    # buildbot > 100 sec tests start here (add tests after lb94.py!)
    'lb94.py',                              # ~84s
    'AA_exx_enthalpy.py',                   # ~119s
    'lcao_tdgllbsc.py',                     # ~132s
    'bse_silicon.py',                       # ~143s
    'gwsi.py',                              # ~147s
    'pw/moleculecg.py',                     # duration unknown
    'potential.py',                         # duration unknown
    'pes.py',                               # duration unknown
    'lcao_pair_and_coulomb.py',             # duration unknown
    'asewannier.py',                        # duration unknown
    'exx_q.py',                             # duration unknown
    'pw/davidson_pw.py',                    # duration unknown
    'neb.py',                               # duration unknown
    'diamond_eps.py',                       # duration unknown
    'wannier_ethylene.py',                  # duration unknown
    'muffintinpot.py',                      # duration unknown
    'nscfsic.py',                           # duration unknown
    'coreeig.py',                           # duration unknown
    'bse_MoS2_cut.py',                      # duration unknown
    'parallel/scalapack_mpirecv_crash.py',  # duration unknown
    'cmrtest/cmr_test.py',                  # duration unknown
    'cmrtest/cmr_test3.py',                 # duration unknown
    'cmrtest/cmr_test4.py',                 # duration unknown
    'cmrtest/cmr_append.py',                # duration unknown
    'cmrtest/Li2_atomize.py']               # duration unknown

# 'fractional_translations.py',
# 'graphene_EELS.py', disabled while work is in progress on response code

# 'fractional_translations_med.py',
# 'fractional_translations_big.py',

# 'eigh_perf.py', # Requires LAPACK 3.2.1 or later
# XXX https://trac.fysik.dtu.dk/projects/gpaw/ticket/230
# 'parallel/scalapack_pdlasrt_hang.py',
# 'dscf_forces.py',
# 'stark_shift.py',


exclude = []

# not available on Windows
if os.name in ['ce', 'nt']:
    exclude += ['maxrss.py']

if mpi.size > 1:
    exclude += ['maxrss.py',
                'pes.py',
                'diamond_eps.py',
                'nscfsic.py',
                'coreeig.py',
                'asewannier.py',
                'wannier_ethylene.py',
                'muffintinpot.py',
                'stark_shift.py',
                'exx_q.py',
                'potential.py',
                # 'cmrtest/cmr_test3.py',
                # 'cmrtest/cmr_append.py',
                # 'cmrtest/Li2_atomize.py',  # started to hang May 2014
                'lcao_pair_and_coulomb.py',
                'bse_MoS2_cut.py',
                'pw/moleculecg.py',
                'pw/davidson_pw.py',
                # scipy.weave fails often in parallel due to
                # ~/.python*_compiled
                # https://github.com/scipy/scipy/issues/1895
                'scipy_test.py']

if mpi.size > 2:
    exclude += ['neb.py']

if mpi.size < 4:
    exclude += ['parallel/fd_parallel.py',
                'parallel/lcao_parallel.py',
                'parallel/pblas.py',
                'parallel/scalapack.py',
                'parallel/scalapack_diag_simple.py',
                'parallel/realspace_blacs.py',
                'AA_exx_enthalpy.py',
                'bse_aluminum.py',
                'bse_diamond.py',
                'bse_silicon.py',
                'bse_vs_lrtddft.py',
                'fileio/parallel.py',
                'parallel/diamond_gllb.py',
                'parallel/lcao_parallel_kpt.py',
                'parallel/fd_parallel_kpt.py']


if mpi.size != 4:
    exclude += ['parallel/scalapack_mpirecv_crash.py']
    exclude += ['parallel/scalapack_pdlasrt_hang.py']

if mpi.size == 1 or not compiled_with_sl():
    exclude += ['parallel/submatrix_redist.py']

if mpi.size != 1 and not compiled_with_sl():
    exclude += ['ralda_energy_H2.py',
                'ralda_energy_N2.py',
                'ralda_energy_Ni.py',
                'ralda_energy_Si.py',
                'bse_sym.py',
                'bse_silicon.py',
                'gwsi.py',
                'rpa_energy_N2.py',
                'pw/fulldiag.py',
                'pw/fulldiagk.py',
                'au02_absorption.py']

if sys.version_info < (2, 6):
    exclude.append('transport.py')
    
if np.__version__ < '1.6.0':
    exclude.append('chi0.py')

exclude = set(exclude)
    

class TestRunner:
    def __init__(self, tests, stream=sys.__stdout__, jobs=1,
                 show_output=False):
        if mpi.size > 1:
            assert jobs == 1
        self.jobs = jobs
        self.show_output = show_output
        self.tests = tests
        self.failed = []
        self.skipped = []
        self.garbage = []
        if mpi.rank == 0:
            self.log = stream
        else:
            self.log = devnull
        self.n = max([len(test) for test in tests])

    def run(self):
        self.log.write('=' * 77 + '\n')
        if not self.show_output:
            sys.stdout = devnull
        ntests = len(self.tests)
        t0 = time.time()
        if self.jobs == 1:
            self.run_single()
        else:
            # Run several processes using fork:
            self.run_forked()

        sys.stdout = sys.__stdout__
        self.log.write('=' * 77 + '\n')
        self.log.write('Ran %d tests out of %d in %.1f seconds\n' %
                       (ntests - len(self.tests) - len(self.skipped),
                        ntests, time.time() - t0))
        self.log.write('Tests skipped: %d\n' % len(self.skipped))
        if self.failed:
            self.log.write('Tests failed: %d\n' % len(self.failed))
        else:
            self.log.write('All tests passed!\n')
        self.log.write('=' * 77 + '\n')
        return self.failed

    def run_single(self):
        while self.tests:
            test = self.tests.pop(0)
            try:
                self.run_one(test)
            except KeyboardInterrupt:
                self.tests.append(test)
                break

    def run_forked(self):
        j = 0
        pids = {}
        while self.tests or j > 0:
            if self.tests and j < self.jobs:
                test = self.tests.pop(0)
                pid = os.fork()
                if pid == 0:
                    exitcode = self.run_one(test)
                    os._exit(exitcode)
                else:
                    j += 1
                    pids[pid] = test
            else:
                try:
                    while True:
                        pid, exitcode = os.wait()
                        if pid in pids:
                            break
                except KeyboardInterrupt:
                    for pid, test in pids.items():
                        os.kill(pid, signal.SIGHUP)
                        self.write_result(test, 'STOPPED', time.time())
                        self.tests.append(test)
                    break
                if exitcode == 512:
                    self.failed.append(pids[pid])
                elif exitcode == 256:
                    self.skipped.append(pids[pid])
                del pids[pid]
                j -= 1

    def run_one(self, test):
        exitcode_ok = 0
        exitcode_skip = 1
        exitcode_fail = 2

        if self.jobs == 1:
            self.log.write('%*s' % (-self.n, test))
            self.log.flush()

        t0 = time.time()
        filename = gpaw.__path__[0] + '/test/' + test

        failed = False
        skip = False

        if test in exclude:
            self.register_skipped(test, t0)
            return exitcode_skip

        try:
            loc = {}
            execfile(filename, loc)
            loc.clear()
            del loc
            self.check_garbage()
        except KeyboardInterrupt:
            self.write_result(test, 'STOPPED', t0)
            raise
        except ImportError, ex:
            module = ex.args[0].split()[-1].split('.')[0]
            if module in ['scipy', 'cmr', '_gpaw_hdf5']:
                skip = True
            else:
                failed = True
        except AttributeError, ex:
            if (ex.args[0] ==
                "'module' object has no attribute 'new_blacs_context'"):
                skip = True
            else:
                failed = True
        except Exception:
            failed = True

        mpi.ibarrier(timeout=60.0)  # guard against parallel hangs

        me = np.array(failed)
        everybody = np.empty(mpi.size, bool)
        mpi.world.all_gather(me, everybody)
        failed = everybody.any()
        skip = mpi.world.sum(int(skip))

        if failed:
            self.fail(test, np.argwhere(everybody).ravel(), t0)
            exitcode = exitcode_fail
        elif skip:
            self.register_skipped(test, t0)
            exitcode = exitcode_skip
        else:
            self.write_result(test, 'OK', t0)
            exitcode = exitcode_ok

        return exitcode

    def register_skipped(self, test, t0):
        self.write_result(test, 'SKIPPED', t0)
        self.skipped.append(test)
    
    def check_garbage(self):
        gc.collect()
        n = len(gc.garbage)
        self.garbage += gc.garbage
        del gc.garbage[:]
        assert n == 0, ('Leak: Uncollectable garbage (%d object%s) %s' %
                        (n, 's'[:n > 1], self.garbage))

    def fail(self, test, ranks, t0):
        if mpi.rank in ranks:
            if sys.version_info >= (2, 4, 0, 'final', 0):
                tb = traceback.format_exc()
            else:  # Python 2.3! XXX
                tb = ''
                traceback.print_exc()
        else:
            tb = ''
        if mpi.size == 1:
            text = 'FAILED!\n%s\n%s%s' % ('#' * 77, tb, '#' * 77)
            self.write_result(test, text, t0)
        else:
            tbs = {tb: [0]}
            for r in range(1, mpi.size):
                if mpi.rank == r:
                    mpi.send_string(tb, 0)
                elif mpi.rank == 0:
                    tb = mpi.receive_string(r)
                    if tb in tbs:
                        tbs[tb].append(r)
                    else:
                        tbs[tb] = [r]
            if mpi.rank == 0:
                text = ('FAILED! (rank %s)\n%s' %
                        (','.join([str(r) for r in ranks]), '#' * 77))
                for tb, ranks in tbs.items():
                    if tb:
                        text += ('\nRANK %s:\n' %
                                 ','.join([str(r) for r in ranks]))
                        text += '%s%s' % (tb, '#' * 77)
                self.write_result(test, text, t0)

        self.failed.append(test)

    def write_result(self, test, text, t0):
        t = time.time() - t0
        if self.jobs > 1:
            self.log.write('%*s' % (-self.n, test))
        self.log.write('%10.3f  %s\n' % (t, text))


if __name__ == '__main__':
    TestRunner(tests).run()
