from __future__ import print_function
import os
import cPickle as pickle
import numpy as np
import matplotlib.pyplot as pl
import tempfile
import smtplib

try:
    from email.mime.text import MIMEText
    from email.mime.multipart import MIMEMultipart
except ImportError:
    from email.MIMEMultipart import MIMEMultipart
    from email.MIMEText import MIMEText

"""Database structure:
dict(testname: [(rev, runtime, info), (rev, runtime, info), ...])
    rev: SVN revision
    runtime: Run time in seconds. Negative for crashed jobs!
    info: A string describing the outcome
"""

class DatabaseHandler:
    """Database class for keeping timings and info for long tests"""
    def __init__(self, filename):
        self.filename = filename
        self.data = dict()

    def read(self):
        if os.path.isfile(self.filename):
            self.data = pickle.load(file(self.filename))
        else:
            print('File does not exist, starting from scratch')

    def write(self, filename=None):
        if filename is None:
            filename = self.filename
        if os.path.isfile(filename):
            os.rename(filename, filename + '.old')
        pickle.dump(self.data, open(filename, 'wb'))

    def add_data(self, name, rev, runtime, info):
        if not self.data.has_key(name):
            self.data[name] = []
        self.data[name].append((rev, runtime, info))

    def get_data(self, name):
        """Return rev, time_array"""
        revs, runtimes = [], []
        if self.data.has_key(name):
            for datapoint in self.data[name]:
                revs.append(datapoint[0])
                runtimes.append(datapoint[1])

        return np.asarray(revs), np.asarray(runtimes)

    def update(self, queue, rev):
        """Add all new data to database"""
        for job in queue.jobs:
            absname = job.absname

            tstart = job.tstart
            if tstart is None:
                tstart = np.nan
            tstop = job.tstop
            if tstop is None:
                tstop = np.nan

            info = job.status

            self.add_data(absname, rev, tstop - tstart, info)

class TestAnalyzer:
    def __init__(self, name, revs, runtimes):
        self.name = name
        self.revs = revs
        self.runtimes = runtimes
        self.better = []
        self.worse = []
        self.relchange = None
        self.abschange = None

    def analyze(self, reltol=0.1, abstol=5.0):
        """Analyze timings

        When looking at a point, attention is needed if it deviates more than
        10\% from the median of previous points. If such a point occurs the
        analysis is restarted.
        """
        self.better = []
        self.worse = []
        abschange = 0.0
        relchange = 0.0
        status = 0
        current_first = 0   # Point to start analysis from
        for i in range(1, len(self.runtimes)):
            tmpruntimes = self.runtimes[current_first:i]
            median = np.median(tmpruntimes[np.isfinite(tmpruntimes)])
            if np.isnan(median):
                current_first = i
            elif np.isfinite(self.runtimes[i]):
                abschange = self.runtimes[i] - median
                relchange = abschange / median
                if relchange < -reltol and abschange < -abstol:
                    # Improvement
                    current_first = i
                    self.better.append(i)
                    status = -1
                elif relchange > reltol and abschange > abstol:
                    # Regression 
                    current_first = i
                    self.worse.append(i)
                    status = 1
                else:
                    status = 0

        self.status = status
        self.abschange = abschange
        self.relchange = relchange * 100

    def plot(self, outputdir=None):
        if outputdir is None:
            return
        fig = pl.figure()
        ax = fig.add_subplot(1, 1, 1)
        ax.plot(self.revs, self.runtimes, 'ko-')
        ax.plot(self.revs[self.better],
                self.runtimes[self.better],
                'go', markersize=8)
        ax.plot(self.revs[self.worse],
                self.runtimes[self.worse],
                'ro', markersize=8)
        ax.set_title(self.name)
        if not outputdir.endswith('/'):
            outputdir += '/'
        figname = self.name.replace('/','_')
        fig.savefig(outputdir + figname + '.png')

class MailGenerator:
    def __init__(self, queue):
        self.better = []
        self.worse = []

        self.FAILED = []
        self.TIMEOUT = []

    def add_test(self, ta):
        if ta.abschange < 0.0:
            self.add_better(ta)
        else:
            self.add_worse(ta)

    def add_better(self, ta):
        self.better.append((ta.name, ta.runtimes[-1],
                            ta.abschange, ta.relchange))

    def add_worse(self, ta):
        self.worse.append((ta.name, ta.runtimes[-1],
                            ta.abschange, ta.relchange))

    def add_failed(self, name):
        self.FAILED.append(name)

    def add_timeout(self, name):
        self.TIMEOUT.append(name)

    def generate_mail(self):
        mail = ''
        if len(self.FAILED):
            mail += 'The following %i tests failed:\n' % len(self.FAILED)
            for name in self.FAILED:
                mail += '%s\n' % name
            mail += '\n'

        if len(self.TIMEOUT):
            mail += 'The following %i tests timed out:\n' % len(self.TIMEOUT)
            for name in self.TIMEOUT:
                mail += '%s\n' % name
            mail += '\n'

        if len(self.FAILED) or len(self.TIMEOUT):
            mail += 'See attached log for details\n\n'

        mail += 'Benchmark results from weekly tests:\n\n'

        if len(self.better):
            mail += 'The following %i tests improved:\n' % len(self.better)
            for test in self.better:
                mail += '%-55s %7.1fs:%+8.1fs (%+7.2f%%)\n' % test
        else:
            mail += 'No tests improved!\n'

        mail += '\n'

        if len(self.worse):
            mail += 'The following %i tests regressed:\n' % len(self.worse)
            for test in self.worse:
                mail += '%-55s %7.1fs:%+8.1fs (%+7.2f%%)\n' % test
        else:
            mail += 'No tests regressed!\n'

        return mail

    def generate_subject(self):
        subject = 'Weekly tests: '
        if len(self.FAILED):
            subject += '%i FAILURES, ' % len(self.FAILED)
        if len(self.TIMEOUT):
            subject += '%i TIMEOUTS, ' % len(self.TIMEOUT)
        subject += '%i regressions, ' % len(self.worse)
        subject += '%i improvements.' % len(self.better)
        return subject

    def send_mail(self, address, server, attachment=None):
        msg = MIMEMultipart()
        msg['Subject'] = self.generate_subject()
        me = os.environ['USER'] + '@' + server
        msg['From'] = me
        msg['To'] = address
        msg.attach(MIMEText(self.generate_mail()))

        if attachment:
            a = MIMEText(open(attachment).read())
            a.add_header('Content-Disposition', 'attachment',
                         filename='status.log')
            msg.attach(a)
        
        s = smtplib.SMTP(server)
        s.sendmail(me, [address], msg.as_string())
        s.quit()


#def csv2database(infile, outfile):
#    """Use this file once to import the old data from csv"""
#    csvdata = np.recfromcsv(infile)
#    db = DatabaseHandler(outfile)
#    for test in csvdata:
#        name = test[0]
#        for i in range(1, len(test) - 1):
#            runtime = float(test[i])
#            info = ''
#            db.add_data(name, 0, runtime, info)
#    db.write()

def analyse(queue, dbpath, outputdir=None, rev=None,
            mailto=None, mailserver=None, attachment=None):
    """Analyse runtimes from testsuite

    Parameters:
        queue: AGTSQueue
            Que to analuze
        dbpath: str
            Path to file storing previous results
        outputdir: str|None
            If str, figures will be put in this dir
        rev: int|None
            GPAW revision. If None time.time() is used
        mailto: str|None
            Mailaddres to send results to. If None, results will be printed to
            stdout.
    """
    if rev is None:
        import time
        rev = time.time()
    db = DatabaseHandler(dbpath)
    db.read()
    db.update(queue, rev)
    db.write()
    mg = MailGenerator(queue)
    for job in queue.jobs:
        name = job.absname
        if job.status == 'success':
            revs, runtimes = db.get_data(name)
            ta = TestAnalyzer(name, revs, runtimes)
            ta.analyze(abstol=25, reltol=0.04)
            if ta.status:
                mg.add_test(ta)
            ta.plot(outputdir)
        elif job.status == 'FAILED':
            mg.add_failed(name)
        elif job.status == 'TIMEOUT':
            mg.add_timeout(name)

    if mailto is not None:
        mg.send_mail(mailto, mailserver, attachment)
    else:
        print(mg.generate_subject())
        print(mg.generate_mail())
