def installed():
    import os
    from ase.test import NotAvailable
    # check if CASTEP_COMMAND is set a environment variable
    if not os.environ.has_key('CASTEP_COMMAND'):
        print("WARNING: Environment variable CASTEP_COMMAND is not set")
        print("Will set CASTEP_COMMAND  = castep for the sake of this test")
        print("Please change it if this does not run castep in your environment")
        os.environ['CASTEP_COMMAND'] = 'castep'

    if not (os.system('which %s > /dev/null 2>&1' % os.environ['CASTEP_COMMAND']) == 0):
        raise NotAvailable("""Could not find CASTEP. If you have it
                              installed make sure, you set the CASTEP_COMMAND
                              environment variable correctly""")
    return True
