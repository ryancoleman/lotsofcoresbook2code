import numpy as np
from ase.io.pupynere import NetCDFFile


class Reader:
    def __init__(self, filename, comm):
        self.nc = NetCDFFile(filename)

    def dimension(self, name):
        return self.nc.dimensions[name]
    
    def __getitem__(self, name):
        value = getattr(self.nc, name)
        if isinstance(value, str):
            try:
                value = eval(value)
            except (SyntaxError, NameError):
                pass
            return value
        else:
            return value[0]

    def has_array(self, name):
        return name in self.nc.variables
    
    def get(self, name, *indices):
        var = self.nc.variables[name]
        if var.shape == ():
            return var.getValue()
        else:
            if var.dimensions[-1] == 'two':
                x = var[indices]
                array = np.empty(x.shape[:-1], complex)
                array.real = x[..., 0]
                array.imag = x[..., 1]
                return array
            else:
                return var[indices]

    def get_reference(self, name, indices):
        return NetCDFReference(self.nc.variables[name], indices)
    
    def close(self):
        self.nc.close()


class NetCDFReference:
    def __init__(self, var, indices):
        self.var = var
        self.indices = indices
        self.cmplx = (var.dimensions[-1] == 'two')
        n = len(indices)
        if self.cmplx:
            self.shape = var.shape[n:-1]
        else:
            self.shape = var.shape[n:]

    def __len__(self):
        return self.shape[0]

    def __getitem__(self, indices):
        if not isinstance(indices, tuple):
            indices = (indices,)
        if self.cmplx:
            x = self.var[self.indices + indices]
            array = np.empty(x.shape[:-1], complex)
            array.real = x[..., 0]
            array.imag = x[..., 1]
            return array
        else:
            return self.var[self.indices + indices]
