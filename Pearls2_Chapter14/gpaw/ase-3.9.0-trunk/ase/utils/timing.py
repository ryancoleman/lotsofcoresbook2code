from __future__ import print_function

# Copyright (C) 2003  CAMP
# Please see the accompanying LICENSE file for further information.

"""A replacement for the ``time.clock()`` function.

From the clock man page::

       Note that the time can wrap around. On a 32bit system
       where CLOCKS_PER_SEC equals 1000000 this function will
       return the same value approximately every 72 minutes.

The ``clock()`` function defined below tries to fix this problem.
However, if the ``clock()`` function is not called often enough (more
than 72 minutes between two calls), then there is no way of knowing
how many times the ``time.clock()`` function has wrapped arround! - in
this case a huge number is returned (1.0e100).  This problem can be
avoided by calling the ``update()`` function at intervals smaller than
72 minutes."""

import sys
import time
import math
import functools

import numpy as np

wrap = 1e-6 * 2**32

# Global variables:
c0 = time.clock()
t0 = time.time()
cputime = 0.0
trouble = False


def clock():
    """clock() -> floating point number

    Return the CPU time in seconds since the start of the process."""

    update()
    if trouble:
        return 1.0e100
    return cputime

    
def update():
    global trouble, t0, c0, cputime
    if trouble:
        return
    t = time.time()
    c = time.clock()
    if t - t0 >= wrap:
        trouble = True
        return
    dc = c - c0
    if dc < 0.0:
        dc += wrap
    cputime += dc
    t0 = t
    c0 = c

    
def function_timer(func, *args, **kwargs):
    out = kwargs.pop('timeout', sys.stdout)
    t1 = time.time()
    r = func(*args, **kwargs)
    t2 = time.time()
    print(t2 - t1, file=out)
    return r

    
class Timer:
    def __init__(self, print_levels=1000):
        self.timers = {}
        self.t0 = time.time()
        self.running = []
        self.print_levels = print_levels
    
    def print_info(self, calc):
        """Override to get to write info during calculator's initialize()."""
        pass

    def start(self, name):
        names = tuple(self.running + [name])
        self.timers[names] = self.timers.get(names, 0.0) - time.time()
        self.running.append(name)
        
    def stop(self, name=None):
        if name is None:
            name = self.running[-1]
        names = tuple(self.running)
        running = self.running.pop()
        if name != running:
            raise RuntimeError('Must stop timers by stack order.  '
                               'Requested stopping of %s but topmost is %s'
                               % (name, running))
        self.timers[names] += time.time()
        return names

    def __call__(self, name):
        """Context manager for timing a block of code.
        
        Example (t is a timer object)::

            with t('Add two numbers'):
                x = 2 + 2
                
            # same as this:
            t.start(Add two numbers')
            x = 2 + 2
            t.stop()
        """
        self.start(name)
        return self
        
    def __enter__(self):
        pass
        
    def __exit__(self, *args):
        self.stop()
            
    def get_time(self, *names):
        return self.timers[names]
                
    def write(self, out=sys.stdout):
        were_running = list(self.running)
        while self.running:
            self.stop()
        if len(self.timers) == 0:
            return

        t0 = time.time()
        tot = t0 - self.t0

        n = max([len(names[-1]) + len(names) for names in self.timers]) + 1
        line = '=' * (n + 26) + '\n'
        out.write(line)
        out.write('%-*s    incl.     excl.\n' % (n, 'Timing:'))
        out.write(line)
        tother = tot
        
        inclusive = self.timers.copy()
        exclusive = self.timers.copy()
        keys = exclusive.keys()
        keys.sort()
        for names in keys:
            t = exclusive[names]
            if len(names) > 1:
                if len(names) < self.print_levels + 1:
                    exclusive[names[:-1]] -= t
            else:
                tother -= t
        exclusive[('Other',)] = tother
        inclusive[('Other',)] = tother
        keys.append(('Other',))
        for names in keys:
            t = exclusive[names]
            tinclusive = inclusive[names]
            r = t / tot
            p = 100 * r
            i = int(40 * r + 0.5)
            if i == 0:
                bar = '|'
            else:
                bar = '|%s|' % ('-' * (i - 1))
            level = len(names)
            if level > self.print_levels:
                continue
            name = (level - 1) * ' ' + names[-1] + ':'
            out.write('%-*s%9.3f %9.3f %5.1f%% %s\n' %
                      (n, name, tinclusive, t, p, bar))
        out.write(line)
        out.write('%-*s%9.3f %5.1f%%\n' % (n + 10, 'Total:', tot, 100.0))
        out.write(line)
        out.write('date: %s\n' % time.asctime())
        
        for name in were_running:
            self.start(name)

    def add(self, timer):
        for name, t in timer.timers.items():
            self.timers[name] = self.timers.get(name, 0.0) + t

            
class timer:
    """Decorator for timing a method call.
    
    Example::

        from ase.utils.timing import timer, Timer

        class A:
            def __init__(self):
                self.timer = Timer()
                
            @timer('Add two numbers')
            def add(self, x, y):
                return x + y
                
        """
    def __init__(self, name):
        self.name = name
    
    def __call__(self, method):
        @functools.wraps(method)
        def new_method(slf, *args, **kwargs):
            slf.timer.start(self.name)
            x = method(slf, *args, **kwargs)
            try:
                slf.timer.stop()
            except IndexError:
                pass
            return x
        return new_method

