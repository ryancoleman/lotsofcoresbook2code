from ase.test import cli
from ase.db import connect

cmd = """
ase-build H | ase-run emt -d y.json &&
ase-build H2O | ase-run emt -d y.json &&
ase-build O2 | ase-run emt -d y.json &&
ase-build H2 | ase-run emt -f 0.02 -d y.json &&
ase-build O2 | ase-run emt -f 0.02 -d y.json &&
ase-build -x fcc Cu | ase-run emt -E 5 -d y.json &&
ase-db y.json natoms=1,Cu=1 --delete --yes &&
ase-db y.json "H>0" -k hydro=1,abc=42,foo=bar &&
ase-db y.json "H>0" --delete-keys foo"""

for name in ['y.json', 'y.db']:  # , 'postgresql://localhost', 'mysql://localhost']:
    cli(cmd.replace('y.json', name))
    con = connect(name)
    assert len(list(con.select())) == 5
    assert len(list(con.select('hydro'))) == 3
    assert con.get_atoms(H=1)[0].magmom == 1
    assert len(list(con.select('foo'))) == 0
    assert len(list(con.select(abc=42))) == 3
    assert len(list(con.select('abc'))) == 3
    assert len(list(con.select('abc,foo'))) == 0
    assert len(list(con.select('abc,hydro'))) == 3
    assert len(list(con.select(foo='bar'))) == 0
    assert len(list(con.select(formula='H2'))) == 1
    assert len(list(con.select(formula='H2O'))) == 1

    id = con.reserve(abc=7)
    assert con[id].abc == 7
