from __future__ import print_function
import numpy as np

def Bootstrap(chi0_wGG, Nw, Kc_GG, printtxt, print_bootstrap, world):

    Nw_local = chi0_wGG.shape[0]
    npw = chi0_wGG.shape[1]
    
    # arxiv 1107.0199
    fxc_GG = np.zeros((npw, npw), dtype=complex)
    tmp_GG = np.eye(npw, npw)
    dminv_wGG = np.zeros((Nw_local, npw, npw), dtype=complex)

    dflocal_w = np.zeros(Nw_local, dtype=complex)
    df_w = np.zeros(Nw, dtype=complex)
                
    for iscf in range(120):
        dminvold_wGG = dminv_wGG.copy()
        Kxc_GG = Kc_GG + fxc_GG

        for iw in range(Nw_local):
            chi_GG = np.dot(chi0_wGG[iw], np.linalg.inv(tmp_GG - np.dot(Kxc_GG, chi0_wGG[iw])))
            dminv_wGG[iw] = tmp_GG + np.dot(Kc_GG, chi_GG)

        if world.rank == 0:
            alpha = dminv_wGG[0,0,0] / (Kc_GG[0,0] * chi0_wGG[0,0,0])
            fxc_GG = alpha * Kc_GG
        world.broadcast(fxc_GG, 0)
    
        error = np.abs(dminvold_wGG - dminv_wGG).sum()
        if world.sum(error) < 0.1:
            printtxt('Self consistent fxc finished in %d iterations ! ' %(iscf))
            break
        if iscf > 100:
            printtxt('Too many fxc scf steps !')
    
        if print_bootstrap:
            for iw in range(Nw_local):
                dflocal_w[iw] = np.linalg.inv(dminv_wGG[iw])[0,0]
            world.all_gather(dflocal_w, df_w)
            if world.rank == 0:
                f = open('df_scf%d' %(iscf), 'w')
                for iw in range(Nw):
                    print(np.real(df_w[iw]), np.imag(df_w[iw]), file=f)
                f.close()
            world.barrier()
        
    for iw in range(Nw_local):
        dflocal_w[iw] = np.linalg.inv(dminv_wGG[iw])[0,0]
    world.all_gather(dflocal_w, df_w)
    
    return df_w
