def agts(queue):
    queue.add('N2Cu_relax.agts.py',
              creates=['N2Cu-N2.csv', 'N2Cu-surf.csv'])

if __name__ == "__main__":
    from ase.optimize.test.N2Cu_relax import *
