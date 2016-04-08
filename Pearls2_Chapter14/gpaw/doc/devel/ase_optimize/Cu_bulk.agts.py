def agts(queue):
    queue.add('Cu_bulk.agts.py',
              creates=['Cu_bulk.csv'])

if __name__ == "__main__":
    from ase.optimize.test.Cu_bulk import *
