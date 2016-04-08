#!/usr/bin/env python

label = 'oscillator'

from os import environ
octopus = environ['OCTOPUS_SCRIPT']
locals = {'label': label}
execfile(octopus, {}, locals)
exitcode = locals['exitcode']
if exitcode != 0:
    raise RuntimeError(('Octopus exited with exit code: %d.  ' +
                        'Check %s.log for more information.') %
                       (exitcode, label))

