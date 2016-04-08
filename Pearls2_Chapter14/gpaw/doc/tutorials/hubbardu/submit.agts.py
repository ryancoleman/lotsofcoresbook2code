def agts(queue):
    queue.add('nio.py')
    n = queue.add('n.py')
    queue.add('check.py', deps=n, creates='gaps.csv')

