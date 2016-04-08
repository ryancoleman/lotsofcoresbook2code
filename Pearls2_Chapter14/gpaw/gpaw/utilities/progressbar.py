from __future__ import print_function, division
import sys
import functools
from time import time, sleep

from ase.utils import devnull

from gpaw.utilities.memory import maxrss


class ProgressBar:
    def __init__(self, fd=sys.stdout, nobar=False):
        """Progress-bar.
        
        Usage::
            
            pb = ProgressBar()
            for i in range(10):
                pb.update(i / 10.0)
                do_stuff()
            pb.finish()
        """
        self.fd = fd
        
        try:
            self.tty = fd.isatty()
        except AttributeError:
            self.tty = False
        
        self.done = False
        self.n = None
        self.nobar = nobar
        self.t = time()
        self.symbols = ['-', '*']
        fd.flush()
        
    def update(self, x):
        """Update progress-bar (0 <= x <= 1)."""
        if x == 0 or self.done:
            return

        if self.tty:
            N = 35
        elif self.nobar:
            N = 10
        else:
            N = 40

        n = int(N * x)
        t = time() - self.t
        t_dt = self.format_time(t)
        est = self.format_time(t / x)
        p = functools.partial(print, file=self.fd)

        if self.tty:
            bar = '-' * (n - 1) + self.symbols[int(t % len(self.symbols))]
            p(('\r{0} / {1} ({2:.0f}%) |{3:' + str(N) + '}| ')
              .format(t_dt, est, x * 100, bar), end='')
            print(' {0:.0f} MB/core'.format(maxrss() / 1024**2), end='')
            if x == 1:
                p()
                self.done = True
            self.fd.flush()
        elif self.nobar:
            if self.n is None:
                p(('Started: {0:.0f} MB/core')
                  .format(maxrss() / 1024**2))
                self.n = 0
            if n > self.n:
                p(('{0} of {1} ({2:.0f}%) {3:.0f} MB/core')
                  .format(t_dt, est, x * 100, maxrss() / 1024**2))
                self.fd.flush()
                self.n = n
            if x == 1:
                p('Finished in {0}'.format(t_dt))
                self.fd.flush()
                self.done = True
        else:
            if self.n is None:
                p('{0}s |'.format(t / x), end='')
                self.n = 0
            if n > self.n:
                p('-' * (n - self.n), end='')
                self.fd.flush()
                self.n = n
            if x == 1:
                p('| Time: {0:.3f}s'.format(t))
                self.fd.flush()
                self.done = True
        
    def format_time(self, seconds):
        m, s = divmod(seconds, 60)
        h, m = divmod(m, 60)
        return "{0:.0f}h{1:.0f}m{2:.0f}s".format(h, m, s)
                
    def finish(self):
        self.update(1)


def test():
    for fd in [sys.stdout, devnull, open('pb.txt', 'w')]:
        print(fd)
        pb = ProgressBar(fd)
        for i in range(20):
            pb.update(i / 20)
            sleep(0.03)
            pb.update((i + 1) / 20)
        pb.finish()


if __name__ == '__main__':
    test()
