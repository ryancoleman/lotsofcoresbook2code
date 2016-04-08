import ase

tol = 1e-9

a = ase.Atoms('CC', [[9.5,5,5],[0.5,5,5]], cell=[10,10,10], pbc=True)

assert abs(a.get_distance(0,1)-9.0) < tol
assert abs(a.get_distance(0,1,mic=True)-1.0) < tol

a.set_distance(0,1, 1.5, mic=True)

assert abs(a.get_distance(0,1)-8.5) < tol
assert abs(a.get_distance(0,1,mic=True)-1.5) < tol

a.set_distance(0,1, 1.5)

assert abs(a.get_distance(0,1)-1.5) < tol
assert abs(a.get_distance(0,1,mic=True)-1.5) < tol
