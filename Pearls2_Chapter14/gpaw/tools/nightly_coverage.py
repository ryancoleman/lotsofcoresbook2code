#!/usr/bin/env python

import os
import sys
import re
import time
import glob
import tempfile
import numpy as np

from subprocess import Popen, PIPE
from trace import pickle
try:
    from cStringIO import StringIO
except ImportError:
    from StringIO import StringIO

def CoverageFilter(url, coverage, filtering, output, digits=4):
    """Filter a coverage file through 'grep -n >>>>>>' and an inlined
    stream editor script in order to convert to reStructuredText."""
    sedname = tempfile.mktemp(prefix='gpaw-sed4rst-')
    htmlfilt = [r's%([\*<._|]{1})%\\\1%g']
    #htmlfilt = ['s%\&%\&amp;%g', 's%>%\&gt;%g', 's%<%\&lt;%g']
    for i in range(1,digits):
        htmlfilt.append('s%^([0-9]){'+str(i)+'}[:-]{1}%'+' '*(digits-i)+'&%g')
    common = r's%^([ \t]*)([0-9]+)[:-]{1}'
    linefilt = ['s%^--$%| --%g',
        common + r'[ \t]*$%| `\2 <' + url + r'#L\2>`__:\1 %g',
        common + r'([ \t]*)([0-9]+):%| `\2 <' + url + r'#L\2>`__:\1 \3*\4* %g',
        common + r'([ \t]*)([>]{6})%| `\2 <' + url + r'#L\2>`__:\1 \3**\4**%g',
        common + r'%| `\2 <' + url + r'#L\2>`__:\1 %g']
    open(sedname, 'w').write('; '.join(htmlfilt + linefilt))
    la = lb = filtering
    if not isinstance(output, str):
        output.flush() # preserve temporal ordering by flushing pipe
        poi = Popen('grep -A%d -B%d -n ">>>>>>" "%s" | sed -rf "%s"' \
            % (la,lb,coverage,sedname), shell=True, stdout=output)
        assert poi.wait() == 0
    elif os.system('grep -A%d -B%d -n ">>>>>>" "%s" | sed -rf "%s" > "%s"' \
        % (la,lb,coverage,sedname,output)) != 0:
        raise RuntimeError('Could not parse coverage file "%s".' % coverage)
    os.system('rm -f "%s"' % sedname)


class CoverageData:
    """Contains coverage information in terms of line-by-line counters."""
    def __init__(self):
        self.nol, self.nos, self.nom, self.noc = {}, {}, {}, {}

    def get(self, devel='*'):
        _nol, _nos, _nom, _noc = [nox.get(devel,0) for nox in \
                                  [self.nol, self.nos, self.nom, self.noc]]
        if _nos == 0:
            ratio, avg = np.nan, np.nan
        else:
            ratio, avg = 1-_nom/float(_nos), _noc/float(_nos)
        return (_nol, _nos, _nom, 100*ratio, _noc, avg)

    def __iadd__(self, other):
        if not isinstance(other, CoverageData):
            raise ValueError('Only instances of CoverageData can be added.')
        for nox, other_nox in zip([self.nol, self.nos, self.nom, self.noc],
                                  [other.nol, other.nos, other.nom, other.noc]):
            for devel, count in other_nox.items():
                nox[devel] = nox.get(devel, 0) + count
        return self

    def ranking(self, reverse=False):
        devels = np.array([devel for devel in self.nol.keys() if devel != '*'])
        totals = np.array([self.nol[devel] for devel in devels])
        return sorted(zip(devels,totals), key=lambda x:x[1], reverse=reverse)

    def width(self):
        return max([len(devel) for devel in self.nol.keys() if devel != '*'])


class CoverageParser(CoverageData):
    """Parses any coverage file, optionally filtered through 'grep -n >>>>>>',
    and convert to reStructuredText while maintaining line-by-line counters."""
    def __init__(self, url, coverage, owners, filtering=None, digits=4):
        CoverageData.__init__(self)
        self.url = url
        if isinstance(coverage, str):
            self.wait, self.pin = self.create_pipe(coverage, filtering)
            self.direct = (filtering is None)
        else:
            self.pin = coverage
            self.direct = False #assumed to be numbered
        self.owners = owners
        self.digits = digits
        self.buf = None
        self.ln = 0

    def ranking(self, reverse=False): # overwrites, numpy is faster
        devels = np.unique(self.owners)
        totals = np.array([np.sum(self.owners==devel) for devel in devels])
        return sorted(zip(devels,totals), key=lambda x:x[1], reverse=reverse)

    def width(self):
        return self.owners.dtype.itemsize

    def create_pipe(self, covername, filtering=None):
        if filtering is None:
            return lambda: 0, open(covername, 'r')
        else:
            la = lb = filtering
            poi = Popen('grep -A%d -B%d -n ">>>>>>" "%s"' \
                % (la,lb,covername), shell=True, stdout=PIPE)
            return poi.wait, poi.stdout

    def increment(self, nox, amount=1):
        devel = self.owners[self.ln]
        nox[devel] = nox.get(devel, 0) + amount
        nox['*'] = nox.get('*', 0) + amount

    def digest(self, pattern, replace, add=False):
        assert pattern.startswith('^')
        m = re.match(pattern, self.buf)
        if m is None:
            return False
        args = (m.group(2),self.url,m.group(2),m.group(1)) + m.groups()[2:]
        self.buf = replace % args + self.buf[m.end():]
        if add:
            self.increment(self.noc, int(m.group(4)))
        return True

    def preprocess(self, line):
        for c in '\*<._|':
            line = line.replace(c, '\\'+c)
        m = re.match('^([0-9]+)([:-]{1}.*)$', line)
        self.ln = int(m.group(1))-1
        if self.digits is None:
            return line
        else:
            return '%*s%s\n' % ((self.digits,)+m.groups())

    def parse(self, line):
        self.buf = self.preprocess(line)
        if self.digest(r'^([ \t]*)([0-9]+)[:-]{1}[ \t]*$', \
                       r'| `%s <%s#L%s>`__:%s '):
            pass # empty line
        elif self.digest(r'^([ \t]*)([0-9]+)[:-]{1}([ \t]*)([>]{6})', \
                         r'| `%s <%s#L%s>`__:%s %s**%s**'):
            # statement with zero count
            self.increment(self.nos)
            self.increment(self.nom)
        elif self.digest(r'^([ \t]*)([0-9]+)[:-]{1}([ \t]*)([0-9]+):', \
                         r'| `%s <%s#L%s>`__:%s %s*%s* ', add=True):
            # statement with non-zero count
            self.increment(self.nos)
        elif self.digest(r'^([ \t]*)([0-9]+)[:-]{1}', r'| `%s <%s#L%s>`__:%s '):
            pass # no countable statement
        else:
            raise IOError('Could not parse "'+self.buf.strip('\n')+'"')
        self.increment(self.nol)
        return self.buf

    def __iter__(self):
        if self.direct:
            # Reading directly from coverage file so prefix line numbers
            # return ('%d:%s' % (l+1,s) for l,s in enumerate(self.pin))
            for ln,line in enumerate(self.pin):
                yield self.parse('%d:%s' % (ln+1, line))
        else:
            # Ignore grep contingency marks if encountered, else right-ajust
            for line in self.pin:
                if line.strip('\n') == '--':
                    yield '| --\n'
                else:
                    yield self.parse(line)
        assert self.wait() == 0
        raise StopIteration

    def write_to_stream(self, pout):
        for line in iter(self):
            pout.write(line)

    def write_to_file(self, filename):
        f = open(filename, 'w')
        self.write_to_stream(f)
        f.close()


class TableIO: # we can't subclass cStringIO.StringIO
    """Formats input into fixed-width text tables for reStructuredText."""
    def __init__(self, widths, simple=False, pipe=None):
        self.widths = widths
        if simple:
            self.seps = {'=':['',' ','\n'], 's':['',' ','\n'], '':['','','']}
            self.ws = ['=', 's', '=', '', '=']
        else:
            self.seps = {'=':['+=', '=+=', '=+\n'], '-':['+-', '-+-', '-+\n'],
                         's':['| ', ' | ', ' |\n'], '':['', '', '']}
            self.ws = ['-', 's', '=', '-', '']
        if pipe is None:
            pipe = StringIO()
        self.pipe = pipe
        self.set_formats(formats=('%-*s',)*len(self.widths), sizes=self.widths)

    def put(self, w, entries, func=None):
        if func is None:
            func = lambda l: w*l
        sep = self.seps[w]
        self.pipe.write(sep[0]+sep[1].join(map(func, entries))+sep[2])

    def set_formats(self, formats, sizes):
        self.formats, self.sizes = formats, sizes

    def add_section(self, w, title):
        self.pipe.write('\n\n%s\n%s\n%s\n\n' % (w*len(title),title,w*len(title)))

    def add_subtitle(self, w, title):
        self.pipe.write('\n%s\n%s\n\n' % (title,w*len(title)))

    def add_heading(self, labels=None):
        if labels is not None:
            self.put(self.ws[0], self.widths)
            self.put(self.ws[1], zip(self.widths,labels), lambda ls: '%-*s' % ls)
        self.put(self.ws[2], self.widths)

    def add_row(self, *args):
        formats, sizes = self.formats, self.sizes
        self.put(self.ws[1], zip(formats,sizes,args), lambda (f,s,a): f % (s,a))
        self.put(self.ws[3], self.widths)

    def write_to_stream(self, pout=None):
        self.put(self.ws[4], self.widths)
        if pout is not None:
            pout.write(self.pipe.getvalue())
            self.pipe.close()


class CoverageIO(TableIO):
    """Formats coverage into fixed-width text tables for reStructuredText."""

    _static_labels = ('NOL', 'NOS', 'NOM', 'Coverage', 'Executions', 'Average')
    _static_formats = ('%*d', '%*d', '%*d', '%*.2f %%', '%*d', '%*.2f')

    def __init__(self, widths, nol=4, nos=4, nom=4, rel=8, noc=11, avg=11, **kwargs):
        TableIO.__init__(self, widths + (nol, nos, nom, rel, noc, avg), **kwargs)
        self.set_formats( \
            formats=('%-*s',)*len(widths) + self._static_formats,
            sizes=widths + (nol, nos, nom, rel-2, noc, avg))

    def add_heading(self, labels=None):
        if labels is not None:
            TableIO.add_heading(self, labels + self._static_labels)
        else:
            TableIO.add_heading(self)

# -------------------------------------------------------------------

def fail(subject, filename='/dev/null'):
    assert os.system('mail -s "%s" christian.glinsvad@fysik.dtu.dk < %s' %
                     (subject, filename)) == 0
    raise SystemExit

def svninfo(url, target, rev='HEAD'):
    """Look up names of the developers who last altered lines of a target."""
    patn = r'^[ \t*]*[0-9]+[ \t]+([^ \t]+)[ \t]+.*$'
    poi = Popen('svn praise -r %s "%s" | sed -r "/%s/!d; s/%s/\\1/g"' \
        % (rev,os.path.join(url,target),patn,patn), shell=True, stdout=PIPE)
    if poi.wait() != 0:
        fail('Lookup of %s failed!' % target)
    return np.array(map(str.strip, poi.stdout.readlines()))

def svnexport(url, path):
    """Exports a clean directory tree from the repository specified by `url`
    into `path` and return the revision number (of the unversioned copy)."""
    poi = Popen('svn export "%s" %s' % (url,path), shell=True, stdout=PIPE)
    if poi.wait() != 0:
        fail('Checkout of %s failed!' % path.upper())
    lastline = poi.stdout.readlines()[-1]
    return re.match('^Exported revision ([0-9]+).$', lastline).group(1)

svnbase = 'https://svn.fysik.dtu.dk/projects/gpaw/trunk'
tmpdir = tempfile.mkdtemp(prefix='gpaw-coverage-')
hostname = os.getenv('HOSTNAME')
cpuruns = [1,2,4,8]
ignoredirs = ['gpaw/db', 'gpaw/sunday', 'gpaw/testing']

if '--rebuild' in sys.argv[1:]:
    # Build .rst files from .cover files without installing or running tests.
    assert os.path.isfile('counts.out')
    f = open('counts.out','r')
    rvs = dict([entry.strip('$ \n').split('=',1) for entry in 
                f.readline().split(';')])
    allfiles = f.readline().strip('$ \n').split(';')
    os.system('cp counts.out "%s"' % tmpdir)
    assert os.path.isdir('coverage')
    os.system('cp -r coverage "%s/coverage"' % tmpdir)
    os.chdir(tmpdir)
else:
    # Get SVN revision numbers, checkout a fresh version and install:
    rvs = {}
    os.chdir(tmpdir)
    rvs['gpaw'] = svnexport(svnbase, 'gpaw')
    rvs['ase'] = svnexport('https://svn.fysik.dtu.dk/projects/ase/trunk', 'ase')
    os.chdir('gpaw')
    f = open('counts.out','w')
    f.write('$ '+';'.join(map('='.join, rvs.items()))+'\n')
    allfiles = []
    for dirpath, dirnames, filenames in os.walk('gpaw'):
        if np.any([dirpath.startswith(d) for d in ignoredirs]):
            continue
        allfiles.extend([os.path.join(dirpath,filename) for filename in \
                         filenames if filename.endswith('.py')])
    f.write('$ '+';'.join(sorted(allfiles))+'\n')
    f.close()

    # Temporary installations of GPAW/ASE revisions and latest setups
    if hostname == 'thul.fysik.dtu.dk':
        customload = 'source /home/camp/modulefiles.sh; ' \
                     'module load NUMPY; ' \
                     'module load openmpi/1.3.3-1.el5.fys.gfortran43.4.3.2; '
        flags = ('--remove-default-flags ' +
                '--customize=' +
                 'doc/install/Linux/Niflheim/el5-xeon-gcc43-acml-4.3.0.py')
        pydir = 'lib64/python'
    elif hostname.endswith('.fysik.dtu.dk'):
        customload = ''
        flags = ''
        customrc = '%s/customize-%s.py' % (os.getenv('HOME'), hostname[:-13])
        if os.path.isfile(customrc):
            flags += ' --customize=' + customrc
        pydir = 'lib/python'
    else:
        raise EnvironmentError('Unknown host. Automatic installation failed.')

    if os.system(customload.replace(';', '&&') +
                 'python setup.py %s install --home=%s ' % (flags,tmpdir) +
                 '2>&1 | grep -v "c/libxc/src"') != 0 \
        or os.system('mv ../ase/ase ../%s' % pydir) != 0:
        fail('Installation failed!')

    os.system('wget --no-check-certificate --quiet ' +
              'http://wiki.fysik.dtu.dk/gpaw-files/gpaw-setups-latest.tar.gz')
    os.system('tar xvzf gpaw-setups-latest.tar.gz')
    setups = tmpdir + '/gpaw/' + glob.glob('gpaw-setups-[0-9]*')[0]

    # Repeatedly run test-suite in code coverage mode:
    args = '--coverage counts.pickle ' + ' '.join(sys.argv[1:])
    for cpus in cpuruns:
        tod = time.strftime('%d/%m-%Y %H:%M:%S')
        open('counts.out', 'a').write('\n\n%s - %d thread(s).\n' % (tod,cpus))
        if os.system(customload +
                     'export PYTHONPATH=%s/%s:$PYTHONPATH; ' % (tmpdir,pydir) +
                     'export GPAW_SETUP_PATH=%s; ' % setups +
                     'export IGNOREPATHS=%s; ' % os.getenv('HOME') +
                     'mpiexec -np %d ' % cpus +
                     tmpdir + '/bin/gpaw-python ' +
                     'tools/gpaw-test %s >>counts.out 2>&1' % args) != 0 \
            or not os.path.isfile('counts.pickle'):
            fail('Test coverage failed!', 'counts.out')

    # Convert pickled coverage information to .cover files with clear text
    if os.system('export PYTHONPATH=%s/%s:$PYTHONPATH; %s -m trace --report ' \
        '--missing --file counts.pickle --coverdir coverage >>counts.out 2>&1' \
        % (tmpdir,pydir,sys.executable)) != 0:
        fail('Coverage conversion failed!', 'counts.out')

    # Create tarball of the data and .cover files and copy it to shared folder
    if os.system('tar cvzf gpaw-counts-%s.tar.gz counts.* ' \
                 'coverage/gpaw.*.cover' % rvs['gpaw']) == 0:
        home = os.getenv('HOME')
        try:
            os.mkdir(home + '/sphinx')
        except OSError:
            pass
        os.system('cp -v --backup=existing gpaw-counts-%s.tar.gz ' \
        '"%s/sphinx/gpaw-counts-latest.tar.gz"' % (rvs['gpaw'],home))

# -------------------------------------------------------------------

# Parse output to test suite logfile and generate reStructuredText
f = open('testsuite.rst', 'w')
f.write('.. _testsuite:\n')
loginfo = TableIO((18, max(10,len(hostname))), simple=True, pipe=f)
loginfo.add_section('=', 'Test suite')
loginfo.add_heading()
for k,v in rvs.items():
    loginfo.add_row('%s revision:' % k.upper(), v)
loginfo.add_row('Build date:', time.strftime('%d/%m-%Y'))
loginfo.add_row('Ran on host:', hostname)
loginfo.write_to_stream()
assert os.system('tail -n+3 counts.out >testsuite.log') == 0
f.write('\n\n.. literalinclude:: testsuite.log\n\n')
f.close()

# Initialize pipes as tables for various categories
rlimits = np.array([0, 50, 90, 100, np.inf])
categories = ['Poor','Mediocre','Good', 'Complete']
pipes = []
l_filename = 42
l_filelink = len(':ref:` <>`') + 2*l_filename-4
for category in categories:
    pipe = CoverageIO((l_filelink,))
    pipe.add_subtitle('-', 'Files with %s coverage' % category.lower())
    pipe.add_heading(('Filename',))
    pipes.append(pipe)

urlbase = 'https://trac.fysik.dtu.dk/projects/gpaw/browser/trunk'
rstindex = 'coverage/index.rst'
tablename = 'coverage/summary.txt'
rankingname = 'coverage/ranking.txt'
devnull = open('/dev/null', 'w', buffering=0)

# Generate a toctree with an alphabetical list of files
f = open(rstindex, 'w')
f.write("""

-----------------------------------
List of files with missing coverage
-----------------------------------

.. toctree::
   :maxdepth: 1

""")

namefilt = re.compile('coverage/(gpaw\..*)\.cover')
cumcd = CoverageData()

for covername in sorted(glob.glob('coverage/gpaw.*.cover')):
    # Find out which .py file in GPAW this .cover file corresponds to
    refname = namefilt.match(covername).group(1)
    rstname = 'coverage/' + refname + '.rst'
    print('cover:', covername, '->', rstname, 'as :ref:`%s`' % refname)
    filename = refname.replace('.','/')+'.py' # unmangle paths
    fileurl = os.path.join(urlbase, '%s?rev=%s' % (filename,rvs['gpaw']))
    filelink = (':ref:`%s <%s>`' % (filename,refname)).ljust(l_filelink)
    if np.any([filename.startswith(d) for d in ignoredirs]):
         print('Ignored...')
         continue
    elif filename not in allfiles:
        raise RuntimeError('%s not in inventory of files.' % filename)
    allfiles.remove(filename)

    # Tally up how many developers have contributed to the file and by how much
    owners = svninfo(svnbase, filename, rvs['gpaw'])
    cp = CoverageParser(fileurl, covername, owners)
    cp.write_to_stream(devnull)
    t_nol, t_nos, t_nom, t_rel, t_noc, t_avg = cp.get()
    cumcd += cp

    # Skip files which are fully covered, i.e. have no ">>>>>>"
    if t_nos == 0:
        continue
    elif t_nom == 0:
        pipes[-1].add_row(filename, *cp.get())
        continue
    else:
        c = np.argwhere((rlimits[:-1]<=t_rel) & (t_rel<rlimits[1:])).item()
        pipes[c].add_row(filelink, *cp.get())

    f.write('   %s\n' % refname) # add to toctree in index file

    # Write reStructuredText file corresponding to coverage file
    g = open(rstname, 'w')
    g.write('\n\n.. _%s:\n' % refname)

    # Add header with a simple table of date and versions
    header = TableIO((22, 10), simple=True, pipe=g)
    header.add_section('=', 'Coverage of %s' % filename)
    header.add_heading()
    for k,v in rvs.items():
        header.add_row('%s revision:' % k.upper(), v)
    header.add_row('Date of compilation:', time.strftime('%d/%m-%Y'))
    header.add_row('Coverage category:', categories[c])
    header.write_to_stream()

    # Add summary table with current coverage information listed by developer
    g.write("""
-------------------------------------
Distribution of coverage by developer
-------------------------------------
The following table summarizes the coverage of the file ``%s`` by the test suite.
For each developer who has made a contribution to the Python file, the coverage
summary specifies the number of lines (NOL), number of executable statements (NOS)
and number of coverage misses (NOM) for the subset of lines he/she has committed.

"""  % filename)
    l_devel = max(len('Developer'), cp.width())
    td = CoverageIO((l_devel,), simple=False, pipe=g)
    td.add_heading(('Developer',))
    for devel,total in cp.ranking(reverse=True):
        td.add_row(devel, *cp.get(devel))
    td.put('-', td.widths)
    td.add_row('*Total:*', *cp.get())
    td.write_to_stream()

    # Add relevant excerpt of the .cover file in nice reStructuredText format
    g.write("""
Line ownership is determined from the second column of::

   svn praise -r %s "%s/%s"

-------------------
Test suite coverage
-------------------
Below is an excerpt of the coverage file ``%s``, which was generated by running
the test suite as detailed in the :ref:`code coverage <coverage>` section. Only
lines marked by `>>>>>>` were not covered, thereby giving rise to concern, but
a few adjacent lines are also listed in order to show the context of the code.
Blue numbers specify the number of times an executable statement has been run.

""" % (rvs['gpaw'],svnbase,filename,covername))
    #CoverageParser(fileurl, covername, owners, 3).write_to_stream(g)
    CoverageFilter(fileurl, covername, 3, g)
    g.close()

# Include the summary and ranking tables in toctree
f.write("""

Back to :ref:`code coverage <coverage>`.

""")
f.close()



# Build summary with tables for the various categories
f = open(tablename, 'w')
f.write("""
-------
Summary
-------

The following tables summarize the coverage of Python files in GPAW
by the test suite. Files are divided into the following %d categories
based on the amount of test suite coverage:

""" % len(categories))
end = 'of the executable statements were covered.\n'
f.write('- :ref:`%s`: Less than %.0f %% %s' % (categories[0].lower(),rlimits[1],end))
for c,category in tuple(enumerate(categories))[1:-1]:
    f.write('- :ref:`%s`: %.0f %% - %.0f %% %s' \
         % (category.lower(),rlimits[c],rlimits[c+1],end))
f.write('- :ref:`%s`: %.0f %% %s' % (categories[-1].lower(),rlimits[-2],end))
f.write("""

For each file, the coverage summary specifies the number of lines (NOL), number
of executable statements (NOS) and number of coverage misses (NOM). Additionaly,
the total number of statement executions is presented, as well as the average
number of executions for each statement. Note that the test suite has been run
multiple times as invidual Python threads, thus causing a total of %d Python
imports for each file. Naturally, Python files which are not imported by the
test suite will generate no coverage information, so these are listed under
:ref:`files ignored <ignored>`.

""" % np.sum(cpuruns))
for category,pipe in zip(categories,pipes):
    f.write('\n\n.. _%s:\n' % category.lower())
    pipe.write_to_stream(f)

f.write("""

.. _ignored:

Files without coverage
----------------------
The following files were ignored because they were not imported:

""")
for filename in allfiles:
    f.write('- %s\n' % filename)
f.close()



# Add another table with cumulated coverage information listed by developer
f = open(rankingname, 'w')
f.write("""
-------------------------------------
Distribution of coverage by developer
-------------------------------------
The following table summarizes the cumulative coverage of the GPAW code base by
the test suite. For each developer who has made a contribution to any of these
Python files, the coverage summary specifies the total number of lines (NOL),
total number of executable statements (NOS) and total number of coverage misses
(NOM) for the subset of lines he/she has committed.

""")
l_devel = max(len('Developer'), cumcd.width())
td = CoverageIO((l_devel,), nol=5, nos=5, nom=5, noc=12, simple=False, pipe=f)
td.add_heading(('Developer',))
for devel,total in cumcd.ranking(reverse=True):
    td.add_row(devel, *cumcd.get(devel))
td.put('-', td.widths)
td.add_row('*Total:*', *cumcd.get())
td.write_to_stream()
f.close()


# Create tarball of the generated .rst files and copy it to shared folder
if os.system('tar cvzf gpaw-coverage-%s.tar.gz testsuite.* ' \
             'coverage/*.rst coverage/*.txt' % rvs['gpaw']) == 0:
    home = os.getenv('HOME')
    try:
        os.mkdir(home + '/sphinx')
    except OSError:
        pass
    assert os.system('cp -v --backup=existing gpaw-coverage-%s.tar.gz ' \
        '"%s/sphinx/gpaw-coverage-latest.tar.gz"' % (rvs['gpaw'],home)) == 0

os.system('cd; rm -r ' + tmpdir)
