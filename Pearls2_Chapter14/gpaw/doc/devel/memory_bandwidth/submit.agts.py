import os

def agts(queue):
    if 0:  # disabled - this test use a mixture of shell and python
        # performs the "best performance" test from
        # https://wiki.fysik.dtu.dk/gpaw/devel/benchmarks.html#memory-benchmark
        #
        # Set the environment variables for your system before running this script!
        # Setting any variable in this script is ignored!
        machine = os.environ.get('MACHINE', 'TEST')
        ncores = int(os.environ.get('NCORES', 8))
        #
        prepare = queue.add('prepare.py',
                            queueopts='-l nodes=1:ppn=1',
                            ncpus=1, walltime=5, deps=[])
        run_numactl = queue.add('run_numactl.py',
                                queueopts='-l nodes=1:ppn=8:xeon8',
                                ncpus=1, walltime=1*60, deps=[prepare])
        analyse = queue.add('analyse.py',
                            queueopts='-l nodes=1:ppn=1',
                            ncpus=1, walltime=5, deps=[run_numactl],
                            creates=['memory_bandwidth_'+machine+'_py.png'],
                            show=['memory_bandwidth_'+machine+'_py.png'])
