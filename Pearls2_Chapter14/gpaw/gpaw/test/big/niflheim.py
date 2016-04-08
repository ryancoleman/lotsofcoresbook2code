import os
import subprocess

from gpaw.test.big.agts import Cluster


class NiflheimCluster(Cluster):
    def __init__(self, asepath='', setuppath='$GPAW_SETUP_PATH'):
        self.asepath = asepath
        self.setuppath = setuppath
        
    def submit(self, job):
        dir = os.getcwd()
        os.chdir(job.dir)

        self.write_pylab_wrapper(job)

        if job.queueopts is None:
            if job.ncpus < 4:
                ppn = '%d:opteron:ethernet' % job.ncpus
                nodes = 1
                arch = 'linux-x86_64-opteron-2.4'
            elif job.ncpus % 8 == 0:
                ppn = '8:xeon5570'
                nodes = job.ncpus // 8
                arch = 'linux-x86_64-xeon-2.4'
            else:
                assert job.ncpus % 4 == 0
                ppn = '4:opteron:ethernet'
                nodes = job.ncpus // 4
                arch = 'linux-x86_64-opteron-2.4'
            queueopts = '-l nodes=%d:ppn=%s' % (nodes, ppn)
        else:
            queueopts = job.queueopts
            arch = 'linux-x86_64-xeon-2.4'

        gpaw_python = os.path.join(dir, 'gpaw', 'build',
                                   'bin.' + arch, 'gpaw-python')

        submit_pythonpath = ':'.join([
            self.asepath,
            dir + '/gpaw',
            '%s/gpaw/build/lib.%s' % (dir, arch),
            '$PYTHONPATH'])

        run_command = '. /home/camp/modulefiles.sh&& '
        run_command += 'module load MATPLOTLIB&& '  # loads numpy, mpl, ...
        run_command += 'module load SCIPY&& '
        run_command += 'module load povray&& '
        run_command += 'module load ABINIT&& '
        run_command += 'module load DACAPO&& '
        run_command += 'module unload SCIENTIFICPYTHON&& '
        run_command += 'module load SCIENTIFICPYTHON/2.8&& '
        run_command += 'module use --append /home/niflheim/dulak/NWchem&& '
        run_command += 'module load NWCHEM/6.1-27.1.x86_64&& '

        if job.ncpus == 1:
            # don't use mpiexec here,
            # this allows one to start mpi inside the *.agts.py script
            run_command += 'module load '
            run_command += 'openmpi/1.3.3-1.el5.fys.gfortran43.4.3.2&& '
            run_command += ' export PYTHONPATH=' + submit_pythonpath + '&&'
            run_command += ' export GPAW_SETUP_PATH=' + self.setuppath + '&&'
            # we want to run GPAW/ASE scripts + gpaw-python with os.system!
            run_command += ' PATH=%s/ase/tools:$PATH' % dir + '&&'
            run_command += ' PATH=%s/gpaw/tools:$PATH' % dir + '&&'
            run_command += ' PATH=%s/gpaw/build/bin.%s:$PATH' % (dir, arch) + '&&'
            # we run other codes with asec
            run_command += ' export ASE_ABINIT_COMMAND="mpiexec abinip < PREFIX.files > PREFIX.log"&&'
        else:
            run_command += 'module load '
            run_command += 'openmpi/1.3.3-1.el5.fys.gfortran43.4.3.2&& '
            run_command += 'mpiexec --mca mpi_paffinity_alone 1'
            run_command += ' -x PYTHONPATH=' + submit_pythonpath
            run_command += ' -x GPAW_SETUP_PATH=' + self.setuppath
            run_command += ' -x OMP_NUM_THREADS=1'

        p = subprocess.Popen(
            ['/usr/local/bin/qsub',
             '-V',
             queueopts,
             '-l',
             'walltime=%d:%02d:00' %
             (job.walltime // 3600, job.walltime % 3600 // 60),
             '-N',
             job.name],
            stdin=subprocess.PIPE, stdout=subprocess.PIPE)

        out, err = p.communicate(
            'touch %s.start\n' % job.name +
            run_command +
            ' %s %s.py %s > %s.output\n' %
            (gpaw_python, job.script, job.args, job.name) +
            'echo $? > %s.done\n' % job.name)
        assert p.returncode == 0
        id = out.split('.')[0]
        job.pbsid = id
        os.chdir(dir)
