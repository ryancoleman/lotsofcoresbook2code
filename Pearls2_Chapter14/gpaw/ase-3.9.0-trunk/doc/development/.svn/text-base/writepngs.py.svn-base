from ase import io

t = io.read('mytrajectory.traj')

for i, s in enumerate(t):
    # rotate to the desired direction
    s.rotate('z', 'x', rotate_cell=True)

    # repeat with keeping old cell
    cell = s.get_cell()
    s = s.repeat((1, 3, 3))
    s.set_cell(cell)

    ofname = str(i) + '.png'
    print('writing', ofname)
    io.write(ofname, s,
             show_unit_cell=True,
             bbox=[-3, -5, 50, 22])  # set bbox by hand, try and error
