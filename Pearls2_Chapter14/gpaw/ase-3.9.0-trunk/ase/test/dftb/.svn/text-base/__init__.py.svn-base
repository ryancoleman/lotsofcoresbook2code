def installed():
    import os
    from ase.test import NotAvailable
    dcmd = os.getenv('DFTB_COMMAND')
    dpre = os.getenv('DFTB_PREFIX')
    if dcmd == None:
        raise NotAvailable('DFTB_COMMAND not defined')
    if dpre == None:
        raise NotAvailable('DFTB_PREFIX not defined (for slater koster files)')
    return True
