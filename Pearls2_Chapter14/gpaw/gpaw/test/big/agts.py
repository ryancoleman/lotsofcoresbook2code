from __future__ import print_function
import os
import sys
import time
import random


jobstates = ['waiting', 'submitted', 'running', 'success', 'FAILED',
             'disabled', 'TIMEOUT']


class AGTSJob:
    def __init__(self, dir, script, queueopts=None,
                 ncpus=1, walltime=10 * 60, deps=None, creates=None,
                 show=None):
        """Advaced GPAW test system job.

        Example:

        >>> job = AGTSJob('doc/devel/256H2O', 'b256H2O.py --sl_default=4,4,16')
        >>> job.dir
        'doc/devel/256H2O'
        >>> job.script
        'b256H2O.py'
        >>> job.args
        '--sl_default=4,4,16'
        >>> job.name
        'b256H2O.py_--sl_default=4,4,16'
        >>> job.absname
        'doc/devel/256H2O/b256H2O.py_--sl_default=4,4,16'
        """

        if ' ' in script:
            script, self.args = script.split(' ', 1)
        else:
            self.args = ''
        pathname = os.path.normpath(os.path.join(dir, script))
        self.dir, self.script = os.path.split(pathname)
        if self.dir == '':
            self.dir = '.'
        self.absname = pathname
        if self.args:
            self.absname += '_' + self.args.replace(' ', '_')
        dir, self.name = os.path.split(self.absname)
        # any string valid for the batch system submit script, e.g.:
        # '-l nodes=2:ppn=4:opteron285'
        self.queueopts = queueopts
        self.ncpus = ncpus
        self.walltime = walltime

        if deps:
            if not isinstance(deps, (list, tuple)):
                deps = [deps]
            self.deps = deps
        else:
            self.deps = []

        if creates and not isinstance(creates, (list, tuple)):
            creates = [creates]
        self.creates = creates

        # Filenames to use for pylab.savefig() replacement of pylab.show():
        if not show:
            show = []
        self.show = show

        if os.path.exists('%s.status' % self.absname):
            self.status = open('%s.status' % self.absname).readline().strip()
        else:
            self.status = 'waiting'

        self.tstart = None
        self.tstop = None
        self.exitcode = None
        self.pbsid = None

    def set_status(self, status):
        self.status = status
        open('%s.status' % self.absname, 'w').write(status + '\n')

    def check_status(self):
        name = self.absname
        if self.status == 'running':
            if time.time() - self.tstart > self.walltime:
                self.set_status('TIMEOUT')
                return 'TIMEOUT'
            if os.path.exists('%s.done' % name):
                self.tstop = os.stat('%s.done' % name).st_mtime
                self.exitcode = int(open('%s.done' % name).readlines()[-1])
                if self.exitcode:
                    self.set_status('FAILED')
                else:
                    self.set_status('success')
                    if self.creates:
                        for filename in self.creates:
                            path = os.path.join(self.dir, filename)
                            if not os.path.isfile(path):
                                self.set_status('FAILED')
                                break
                return self.status

        elif self.status == 'submitted' and os.path.exists('%s.start' % name):
            self.tstart = os.stat('%s.start' % name).st_mtime
            self.set_status('running')
            return 'running'

        # Nothing happened:
        return None

    def clean(self):
        for name in ['start', 'done', 'status']:
            try:
                os.remove(self.absname + '.' + name)
            except OSError:
                pass

        self.status = 'waiting'


class Cluster:
    def write_pylab_wrapper(self, job):
        """Use Agg backend and prevent windows from popping up."""
        fd = open(job.script + '.py', 'w')
        fd.write('from gpaw.test import wrap_pylab\n')
        fd.write('wrap_pylab(%s)\n' % job.show)
        fd.write('execfile(%r)\n' % job.script)
        fd.close()
        
    def tick(self, nrunning):
        pass


class TestCluster(Cluster):
    def submit(self, job):
        if random.random() < 0.05:
            # randomly fail some of the jobs
            exitcode = 1
        else:
            exitcode = 0

        wait = random.randint(1, 12)
        cmd = 'sleep %s; touch %s.start; ' % (wait, job.absname)

        if random.random() < 0.05:
            # randomly time out some of the jobs
            pass
        else:
            duration = random.randint(4, 12)
            cmd += 'sleep %s; ' % duration
            if exitcode == 0 and job.creates:
                for filename in job.creates:
                    cmd += 'echo 42 > %s; ' % os.path.join(job.dir, filename)
            cmd += 'echo %d > %s.done' % (exitcode, job.absname)
        os.system('(%s)&' % cmd)


class LocalCluster(Cluster):
    def __init__(self):
        self.queue = []
        
    def submit(self, job):
        self.queue.append(job)
        
    def tick(self, nrunning):
        if self.queue and nrunning == 0:
            job = self.queue.pop(0)
            dir = os.getcwd()
            os.chdir(job.dir)
            self.write_pylab_wrapper(job)
            os.system('(touch %s.start;' % job.name +
                      'python %s.py %s > %s.output;' %
                      (job.script, job.args, job.name) +
                      'echo $? > %s.done)&' % job.name)
            os.chdir(dir)


class AGTSQueue:
    def __init__(self, sleeptime=60, fd=sys.stdout):
        self.sleeptime = sleeptime
        self.jobs = []

        if isinstance(fd, str):
            self.fd = open(fd, 'w')
        else:
            self.fd = fd

        # used by add() method:
        self._dir = None

    def count(self):
        self.N = dict([(state, 0) for state in jobstates])
        for j in self.jobs:
            self.N[j.status] += 1

    def log(self, job):
        self.count()
        self.fd.write('%s %3d %3d %3d %3d %3d %3d %3d %-49s %s\n' %
                      ((time.strftime('%H:%M:%S'),) +
                       tuple([self.N[state] for state in jobstates]) +
                       (job.absname, job.status)))
        self.fd.flush()
        self.status()

    def add(self, script, dir=None, queueopts=None, ncpus=1, walltime=15,
            deps=None, creates=None, show=None):
        """Add job.

        XXX move docs from doc/devel/testing to here and use Sphinx autodoc."""

        if dir is None:
            dir = self._dir
        job = AGTSJob(dir, script, queueopts, ncpus,
                      walltime * 60 + self.sleeptime,
                      deps, creates, show)
        self.jobs.append(job)
        return job

    def locate_tests(self):
        for root, dirs, files in os.walk('.'):
            if root.startswith('./_'):
                continue
            for fname in files:
                if fname.endswith('.agts.py'):
                    yield root, fname

    def collect(self):
        """Find agts.py files and collect jobs."""
        for dir, agtsfile in self.locate_tests():
            _global = {}
            execfile(os.path.join(dir, agtsfile), _global)
            self._dir = dir
            _global['agts'](self)
        self.normalize()

    def normalize(self):
        """Convert string dependencies to actual job objects."""
        for job in self.jobs:
            for i, dep in enumerate(job.deps):
                if not isinstance(dep, AGTSJob):
                    absname = os.path.normpath(os.path.join(job.dir, dep))
                    job.deps[i] = self.find(absname)

    def find(self, absname):
        """Find job with a particular name."""
        for job in self.jobs:
            if job.absname == absname:
                return job
        raise ValueError

    def find_job_and_dependencies(self, job):
        if isinstance(job, str):
            for j in self.jobs:
                if j.script == job:
                    job = j
                    break
            else:
                raise ValueError

        jobs = [job]
        for dep in job.deps:
            jobs += self.find_job_and_dependencies(dep)
        return jobs

    def run(self, cluster):
        """Run jobs and return the number of unsuccessful jobs."""
        self.fd.write('time       W   S   R   +   -   .   T job\n')
        jobs = self.jobs
        while True:
            done = True
            for job in jobs:
                if job.status == 'waiting':
                    done = False
                    ready = True
                    for dep in job.deps:
                        if dep.status != 'success':
                            ready = False
                            break
                    if ready:
                        cluster.submit(job)
                        job.set_status('submitted')
                        self.log(job)
                elif job.status in ['running', 'submitted']:
                    done = False

            if done:
                break

            time.sleep(self.sleeptime)

            nrunning = 0
            for job in jobs:
                newstatus = job.check_status()
                if newstatus:
                    self.log(job)
                    if newstatus in ['TIMEOUT', 'FAILED']:
                        self.fail(job)
                if job.status == 'running':
                    nrunning += 1

            cluster.tick(nrunning)
                        
        t = self.get_cpu_time()
        self.fd.write('CPU time: %d:%02d:%02d\n' %
                      (t // 3600, t // 60 % 60, t % 60))

        return len([None for job in self.jobs if job.status != 'success'])

    def status(self):
        fd = open('status.log', 'w')
        fd.write('# job                                              ' +
                 20 * ' ' +
                 'status      time   tmax ncpus  deps files id\n')
        for job in self.jobs:
            if job.tstop is not None:
                t = '%5d' % round(job.tstop - job.tstart)
            else:
                t = '     '
            if job.creates:
                c = len(job.creates)
            else:
                c = 0
            if job.pbsid is not None:
                id = job.pbsid
            else:
                id = ''
            fd.write('%-70s %-10s %s %6d %5d %5d %5d %s\n' %
                     (job.absname, job.status, t, job.walltime,
                      job.ncpus, len(job.deps), c, id))
        fd.close()

    def fail(self, dep):
        """Recursively disable jobs depending on failed job."""
        for job in self.jobs:
            if dep in job.deps:
                job.set_status('disabled')
                self.log(job)
                self.fail(job)

    def clean(self):
        for job in self.jobs:
            job.clean()

    def copy_created_files(self, dir):
        for job in self.jobs:
            if job.creates:
                for filename in job.creates:
                    path = os.path.join(job.dir, filename)
                    if os.path.isfile(path):
                        os.system('cp %s %s' %
                                  (path, os.path.join(dir, filename)))

    def get_cpu_time(self):
        """Calculate CPU time in seconds."""
        t = 0
        for job in self.jobs:
            if job.tstop is not None:
                t += job.ncpus * (job.tstop - job.tstart)
        return t


def main():
    from optparse import OptionParser

    parser = OptionParser(usage='%prog [options] [jobs]',
                          version='%prog 0.1')
    parser.add_option('-c', '--clean', action='store_true')
    parser.add_option('-r', '--run')

    opt, args = parser.parse_args()

    queue = AGTSQueue()

    # Find all jobs:
    queue.collect()

    if args:
        # Select specific job(s) and dependencies:
        newjobs = set()
        for arg in args:
            jobs = queue.find_job_and_dependencies(arg)
            for job in jobs:
                newjobs.add(job)
        queue.jobs = list(newjobs)

    if opt.clean:
        queue.clean()

    for job in queue.jobs:
        job.check_status()

    queue.count()
    for state in jobstates:
        print(('%9s %d' % (state, queue.N[state])))

    if opt.run:
        for job in queue.jobs:
            if job.status in ['FAILED', 'TIMEOUT', 'disabled']:
                job.clean()
                job.set_status('waiting')

        if opt.run == 'test':
            # Quick test using dummy cluster and timeout after only 20 seconds:
            queue.sleeptime = 2.0
            for job in queue.jobs:
                job.walltime = 20
            cluster = TestCluster()
        elif opt.run == 'niflheim':
            from gpaw.test.big.niflheim import NiflheimCluster
            cluster = NiflheimCluster()
        elif opt.run == 'local':
            cluster = LocalCluster()
        else:
            1 / 0

        queue.run(cluster)


if __name__ == '__main__':
    main()
