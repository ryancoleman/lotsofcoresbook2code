import sys

def agts(queue):
    if sys.version_info > (2, 5):
        # needs itertools.combinations introduced in 2.6
        calc = queue.add('fc_butadiene.py', walltime=30)
    else:
        pass
    
 
