import os
import sys
import glob
import shutil
import subprocess

def cmd(c):
    x = os.system(c)
    assert x == 0, c

def fail(subject, email=None, filename='/dev/null', mailer='mail'):
    assert mailer in ['mailx', 'mail', 'mutt']
    import os
    if email is not None:
        if filename == '/dev/null':
            assert os.system('mail -s "%s" %s < %s' %
                             (subject, email, filename)) == 0
        else: # attachments
            filenames = filename.split()
            if mailer == 'mailx': # new mailx (12?)
                attach = ''
                for f in filenames:
                    attach += ' -a %s ' % f
                # send with empty body
                assert os.system('echo | mail %s -s "%s" %s' %
                                 (attach, subject, email)) == 0
            elif mailer == 'mail': # old mailx (8?)
                attach = '('
                for f in filenames:
                    ext = os.path.splitext(f)[-1]
                    if ext:
                        flog = os.path.basename(f).replace(ext, '.log')
                    else:
                        flog = f
                    attach += 'uuencode %s %s&&' % (f, flog)
                # remove final &&
                attach = attach[:-2]
                attach += ')'
                assert os.system('%s | mail -s "%s" %s' %
                                 (attach, subject, email)) == 0
            else: # mutt
                attach = ''
                for f in filenames:
                    attach += ' -a %s ' % f
                # send with empty body
                assert os.system('mutt %s -s "%s" -c %s < /dev/null' %
                                 (attach, subject, email)) == 0
    raise SystemExit

if '--dir' in sys.argv:
    i = sys.argv.index('--dir')
    dir = os.path.abspath(sys.argv[i+1])
else:
    dir = 'agts'

if '--email' in sys.argv:
    i = sys.argv.index('--email')
    email = sys.argv[i+1]
else:
    email = None

assert os.path.isdir(dir)

gpawdir = os.path.join(dir, 'gpaw')

# remove the old run directory
if os.path.isdir(dir):
    shutil.rmtree(dir)

os.mkdir(dir)
os.chdir(dir)

cmd('svn checkout https://svn.fysik.dtu.dk/projects/gpaw/trunk gpaw')

# a hack: link libxc for the common version built on surt statically
cmd('cp gpaw/customize.py gpaw/customize.py.libxca')
cmd('echo extra_link_args += [os.environ[\\\'LIBXC_HOME\\\'] + \\\'/lib/libxc.a\\\'] >> gpaw/customize.py.libxca')
cmd('echo libraries.remove\(\\\'xc\\\'\) >> gpaw/customize.py.libxca')

# a version of gpaw is needed for imports from within this script!
cmd("\
cd " + gpawdir + "&& \
source /home/opt/modulefiles/modulefiles_el6.sh&& \
module load libxc&& \
python setup.py build_ext \
--customize=customize.py.libxca 2>&1 > build_ext.log")

# import gpaw from where it was installed
sys.path.insert(0, gpawdir)

cmd("echo '\
cd '" + gpawdir + "'&& \
source /home/opt/modulefiles/modulefiles_el6.sh&& \
module load intel-compilers && \
python setup.py --remove-default-flags --customize=\
doc/install/Linux/Niflheim/el6-sl230s-tm-gfortran-openmpi-1.6.3-acml-4.4.0-sl-hdf5-1.8.10.py \
build_ext 2>&1 > surt.log' | ssh surt bash")

cmd("echo '\
cd '" + gpawdir + "'&& \
source /home/opt/modulefiles/modulefiles_el6.sh&& \
module load intel-compilers && \
source /home/camp/modulefiles.sh&& \
python setup.py --remove-default-flags --customize=\
doc/install/Linux/Niflheim/el6-dl160g6-tm-gfortran-openmpi-1.6.3-acml-4.4.0-sl-hdf5-1.8.10.py \
build_ext 2>&1 > muspel.log' | ssh muspel bash")

cmd("echo '\
cd '" + gpawdir + "'&& \
source /home/opt/modulefiles/modulefiles_el6.sh&& \
module load intel-compilers && \
source /home/camp/modulefiles.sh&& \
python setup.py --remove-default-flags --customize=\
doc/install/Linux/Niflheim/el6-x3455-tm-gfortran-openmpi-1.6.3-acml-4.4.0-sl-hdf5-1.8.10.py \
build_ext 2>&1 > slid.log' | ssh slid bash")

cmd("""wget --no-check-certificate --quiet \
http://wiki.fysik.dtu.dk/gpaw-files/gpaw-setups-latest.tar.gz && \
tar xzf gpaw-setups-latest.tar.gz && \
rm gpaw-setups-latest.tar.gz && \
mv gpaw-setups-[0-9]* gpaw/gpaw-setups""")

cmd('svn export https://svn.fysik.dtu.dk/projects/ase/trunk ase')

# ase needed
sys.path.insert(0, '%s/ase' % dir)

from gpaw.test.big.agts import AGTSQueue
from gpaw.test.big.niflheim6 import NiflheimCluster

queue = AGTSQueue()
queue.collect()
cluster = NiflheimCluster(asepath=os.path.join(dir, 'ase'),
                          setuppath=os.path.join(gpawdir, 'gpaw-setups'))
# Example below is confusing: job.script must NOT be the *.agts.py script,
# but the actual python script to be run!
# testsuite.agts.py does both: see gpaw/test/big/miscellaneous/testsuite.agts.py
#queue.jobs = [job for job in queue.jobs if job.script == 'testsuite.agts.py']

nfailed = queue.run(cluster)

gfiles = os.path.join(dir, 'agts-files')
if not os.path.isdir(gfiles):
    os.mkdir(gfiles)

queue.copy_created_files(gfiles)

# make files readable by go
files = glob.glob(gfiles + '/*')
for f in files:
    os.chmod(f, 0644)

from gpaw.version import version

subject = 'AGTS GPAW %s: ' % str(version)
# Send mail:
sfile = os.path.join(dir, 'status.log')
attach = sfile
if not nfailed:
    subject += ' succeeded'
    fail(subject, email, attach, mailer='mutt')
else:
    subject += ' failed'
    # attach failed tests error files
    ft = [l.split()[0] for l in open(sfile).readlines() if 'FAILED' in l]
    for t in ft:
        ef = glob.glob(os.path.join(dir, t) + '.e*')
        for f in ef:
            attach += ' ' + f
    fail(subject, email, attach, mailer='mutt')

if 0:
    # Analysis:
    import matplotlib
    matplotlib.use('Agg')
    from gpaw.test.big.analysis import analyse
    user = os.environ['USER']
    analyse(queue,
            '../analysis/analyse.pickle',  # file keeping history
            '../analysis',                 # Where to dump figures
            rev=niflheim.revision,
            #mailto='gpaw-developers@listserv.fysik.dtu.dk',
            mailserver='servfys.fysik.dtu.dk',
            attachment='status.log')
