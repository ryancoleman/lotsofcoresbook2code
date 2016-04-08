def agts(queue):
    #generate = queue.add('generate.py', ncpus=1, walltime=20)
    G = [queue.add('g21gpaw.py %d' % i, walltime=40 * 60)
         for i in range(4)]
    N = [queue.add('g21nwchem.py %d' % i, walltime=10 * 60,
                   queueopts='-l nodes=1:ppn=4:opteron4', ncpus=1)
         for i in range(4)]
    queue.add('analyse.py', deps=G + N, creates='g2-1.csv')
