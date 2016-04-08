def agts(queue):
    queue.add('dos.agts.py', ncpus=1,
              deps=['../iron/PBE.py',
                    '../wavefunctions/CO.py',
                    '../wannier/si.py',
                    '../aluminium/Al_fcc.py'])

if __name__ == '__main__':
    import os
    import sys
    for filename in ['aluminium/Al-fcc.gpw',
                     'wannier/si.gpw',
                     'wavefunctions/CO.gpw',
                     'iron/ferro.gpw',
                     'iron/anti.gpw',
                     'iron/non.gpw']:
        sys.argv = ['', '../' + filename]
        execfile('dos.py')
    if not os.path.islink('ferro.gpw'):
        os.symlink('../iron/ferro.gpw', 'ferro.gpw')
    execfile('pdos.py')
