def agts(queue):
    a = queue.add('dcdft_gpaw.py', ncpus=4, walltime=20)
    queue.add('testdb.py', deps=a)
    queue.add('extract.py dcdft.db', creates='dcdft.db_raw.txt')
