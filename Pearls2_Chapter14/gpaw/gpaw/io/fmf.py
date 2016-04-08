import os
import platform
import time

from gpaw.version import version

class FMF:
    """Full-Metadata Format

    Full-Metadata Format after 
    http://www.sir.uni-freiburg.de/repository/2009/SI20090302a/SI20090302a.pdf"""
    def __init__(self, 
                 title='-',
                 creator=None,
                 place=None,
                 escape='#'):
        self.escape = escape
        self.estimate_creator(creator)
        self.title = title
        self.place = place
 
    def header(self):
        ec = self.escape
        header =  ec + '; -*- fmf version: 1.0 -*-\n'
        header += ec + '[*reference]\n'
        header += ec + 'creator: ' + self.creator + '\n'
        header += ec + 'created: ' + time.strftime('%Y-%m-%d %H:%M') + '\n'

        title = self.title
        if isinstance(title, str):
            title = title.split('\n')
        for i, line in enumerate(title):
            if i == 0:
                header += ec + 'title: '
            else:
                header += ec + '       '
            header += line + '\n'

        header += ec + 'gpaw-version: ' + version + '\n'
        try:
            import socket
            header += ec + 'hostname: ' + socket.gethostname() + '\n'
        except:
            pass
        header += ec + 'architecture: ' + platform.uname()[4]  + '\n'

        return header

    def data(self, definitions):
        ec = self.escape
        data = ec + '[* data definitions]\n'
        for definition in definitions:
            data += ec + definition + '\n'
        data += ec + '[* data]\n'
        return data

    def field(self, title, entries):
        ec = self.escape
        res = ec + '[' + title + ']\n'
        for entry in entries:
            res += ec + entry + '\n'
        return res

    def estimate_creator(self, creator=None):
        if creator is not None:
            self.creator = creator
            return

        try:
            # get 
            import getpass
            username = getpass.getuser()
            
            try:
                import pwd
                gecos = pwd.getpwnam(username).pw_gecos
                hasalpha = False
                for letter in gecos:
                    if letter.isalpha():
                        hasalpha = True
                if hasalpha:
                    creator = gecos
                else:
                    creator = username
            except:
                creator = username

        except:
            creator = 'unknown'

        self.creator = creator

