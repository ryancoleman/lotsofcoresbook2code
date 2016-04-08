#!/usr/bin/env python
from __future__ import print_function
import os
import sys
import xml.sax
from optparse import OptionParser
from math import exp, log


class Reader(xml.sax.handler.ContentHandler):
    def __init__(self, list, name, state):
        self.list = list
        self.name = name
        self.state = state
        self.grid = {}
        xml.sax.handler.ContentHandler.__init__(self)

    def read(self, filename):
        if filename.endswith('.gz'):
            import gzip
            source = gzip.open(filename)
        else:
            source = open(filename)
        xml.sax.parse(source, self)

    def startElement(self, name, attrs):
        if name == 'state' and self.list:
            print(attrs['id'], file=sys.stderr)
        elif name == 'radial_grid':
            istart = int(attrs['istart'])
            iend = int(attrs['iend'])
            ii = range(istart, iend + 1)
            id = attrs['id']
            eq = attrs['eq']
            if eq == 'r=a*exp(d*i)':
                a = float(attrs['a'])
                d = float(attrs['d'])
                self.grid[id] = [a * exp(d * i) for i in ii]
            elif eq == 'r=a*i/(n-i)':
                a = float(attrs['a'])
                n = float(attrs['n'])
                self.grid[id] = [a * i / (n - i) for i in ii]
            elif eq == 'r=a*(exp(d*i)-1)':
                a = float(attrs['a'])
                d = float(attrs['d'])
                self.grid[id] = [a * (exp(d * i) - 1.0) for i in ii]
            elif eq == 'r=d*i':
                d = float(attrs['d'])
                self.grid[id] = [d * i for i in ii]
            elif eq == 'r=(i/n+a)^5/a-a^4':
                a = float(attrs['a'])
                n = float(attrs['n'])
                self.grid[id] = [(i / n + a)**5 / a - a**4 for i in ii]
            else:
                raise RuntimeError('Unknown grid type: ' + attrs['eq'])
        elif name == self.name and attrs.get('state', None) == self.state:
            self.r = self.grid[attrs['grid']]
            self.data = ''
        else:
            self.data = None

    def characters(self, data):
        if self.data is not None:
            self.data += data

    def endElement(self, name):
        if self.data is not None:
            for r, x in zip(self.r, self.data.split()):
                print(r, x)
                

op = OptionParser(usage='%prog [options] setup[.gz]',
                      version='%prog 0.2')
op.add_option('-x', '--extract',
                  help='Function to extract.',
                  metavar='<name>')
op.add_option('-s', '--state',
                  help='Select valence state.',
                  metavar='<channel>')
op.add_option('-l', '--list', action='store_true',
                  help='List valence states.')
options, args = op.parse_args()
if len(args) != 1:
        op.error('incorrect number of arguments')


reader = Reader(options.list, options.extract, options.state)
reader.read(args[0])
