def h_define_filter(s):
    return s.startswith('#define  XC')

def define_functionals_from_h_file(version,
                                   hfile='c/libxc/src/xc_funcs.h',
                                   pyfile='gpaw/xc/libxc_functionals.py'):
    """Extract the values of #define from an c include hfile."""
    from os.path import abspath
    hf = open(abspath(hfile), 'r')
    lines = filter(h_define_filter, hf.readlines())
    hf.close()
    assert len(lines) > 0
    pf = open(abspath(pyfile), 'w')
    pf.write('# Computer generated code! Hands off!\n')
    pf.write('# libxc: version ' + str(version) + '\n')
    pf.write('# http://www.tddft.org/programs/octopus/wiki/index.php/Libxc\n')
    # open dictionary
    pf.write('libxc_functionals = {\n')
    # Put all the defines into the dictionary
    for line in lines:
        # extract the define
        splitline = line.split(None, 2)
        # extract the values stripping c comments
        (name, number, descr) = (
            splitline[1],
            int(splitline[2].split(None, 1)[0]),
            splitline[2].split(None, 1)[1].strip())
        name = name.replace('XC_', '', 1) # remove first "XC_"
        pf.write("'" + str(name) + "'" + ": " + str(number) + "," +
                 "# " + str(descr) + "\n")
    # close dictionary
    pf.write('}\n')
    pf.close()

define_functionals_from_h_file(version='1.2.0')
