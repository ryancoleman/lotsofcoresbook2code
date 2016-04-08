def agts(queue):
    queue.add('CO_Au111.agts.py',
              creates=['CO_Au111.csv'])

if __name__ == "__main__":
    from ase.optimize.test.CO_Au111 import *
