import threading

import numpy as np


def run(task, nthreads, *args):
    """Parallelize task over first index of args."""
    
    N = len(args[0])

    def target(*args):
        for x in zip(*args):
            task(*x)

    if nthreads == 1:
        target(*args)
        return [0, N]
        
    n = np.linspace(0, N, nthreads + 1).round().astype(int)
    
    threads = []
    for n1, n2 in zip(n[:-1], n[1:]):
        thread = threading.Thread(target=target, args=[x[n1:n2] for x in args])
        thread.start()
        threads.append(thread)
    
    for thread in threads:
        thread.join()

    return n
