from gpaw.grid_descriptor import GridDescriptor
gd = GridDescriptor([4, 4, 4])
a = gd.empty(dtype=complex)
a[:] = 1.0
assert gd.integrate(a.real, a.real) == 1.0
