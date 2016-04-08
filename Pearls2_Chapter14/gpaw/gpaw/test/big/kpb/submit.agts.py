def agts(queue):
    molecules = [queue.add('molecules.py %d' % i,
                           ncpus=1,
                           walltime=2*60)
                 for i in range(2)]
    queue.add('check.py', deps=molecules)
 
