def installed():
    import os
    from ase.test import NotAvailable
    vcmd = os.getenv('VASP_COMMAND')
    vscr = os.getenv('VASP_SCRIPT')
    if vcmd == None and vscr == None:
        raise NotAvailable('Neither VASP_COMMAND nor VASP_SCRIPT defined')
    return True
