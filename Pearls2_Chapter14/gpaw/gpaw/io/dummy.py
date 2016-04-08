from gpaw.mpi import world

class DummyWriter:
    def __init__(self, name, comm=world):
        self.comm = comm # for possible future use
        self.master = (self.comm == 0)
        # Only master does meaningful I/O
        if self.master:
            raise RuntimeError('Dummy writer should not have master!')
        
    def dimension(self, name, value):
        pass

    def __setitem__(self, name, value):
        pass
 
    def add(self, name, shape, array=None, dtype=None, units=None,
            parallel=False, write=True):
        pass

    def fill(self, array, *indices, **kwargs):
        pass
 
    def close(self):
        pass



