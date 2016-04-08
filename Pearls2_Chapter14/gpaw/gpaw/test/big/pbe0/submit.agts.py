from __future__ import print_function
def agts(queue):
    gaps = queue.add('gaps.py', ncpus=16, walltime=5 * 60)
    queue.add('submit.agts.py', deps=gaps)
 

if __name__ == '__main__':
    import ase.db
    c = ase.db.connect('gaps.db')
    for d in sorted(c.select(), key=lambda d: d.name):
        print((d.name))
        for k, e in d.data.items():
            r = e[:2].tolist() + (e[:2]-e[2:4]).tolist() + (e[:2]-e[4:]).tolist()
            print((k, ' '.join(['%6.3f' % x for x in r])))
