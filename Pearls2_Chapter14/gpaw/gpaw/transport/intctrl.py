import numpy as np

class PathInfo:
    def __init__(self, type, nlead):
        self.type = type
        self.num = 0
        self.lead_num = nlead
        self.energy = []
        self.weight = []
        self.nres = 0
        self.sigma = []
        for i in range(nlead):
            self.sigma.append([])
        if type == 'eq':
            self.fermi_factor = []
        elif type == 'ne':
            self.fermi_factor = []
            for i in range(nlead):
                self.fermi_factor.append([[], []])
        else:
            raise TypeError('unkown PathInfo type')

    def add(self, elist, wlist, flist, siglist):
        self.num += len(elist)
        self.energy += elist
        self.weight += wlist
        if self.type == 'eq':
            self.fermi_factor += flist
        elif self.type == 'ne':
            for i in range(self.lead_num):
                for j in [0, 1]:
                    self.fermi_factor[i][j] += flist[i][j]
        else:
            raise TypeError('unkown PathInfo type')
        for i in range(self.lead_num):
            self.sigma[i] += siglist[i]

    def set_nres(self, nres):
        self.nres = nres
        
class IntCtrl:
    """
    Parameters:
        kt             // Temperature : default = 0 K
        leadfermi(nLead)  // Fermi Energy of each Lead with biasV
        maxfermi, minfermi   // default = max(Fermi) & min(Fermi)    
        bias               // >0 for ul>u_r
        eqintpath,  eqinttol,  eqdelta,  eqresz     // eqInt Ctrl
        locintpath, lcointtol, locdelta, locresz    // locInt Ctrl

        neintmethod  // 0 - Manual Method : Linear (default)// 1 - Auto Method
        neintstep    // used in Manual Method, defalut = 1e-2
        neinttol     // used in Auto Method, default = 1e-3
        neintpath    // [ minfermi leadfermi maxfermi ] + eta ( 1e-8 )}
    """
    
    def __init__(self, kt, efermi, bias, min_energy=-100,
                                                               neintmethod=0,
                                                               neintstep=0.02,
                                                               eqinttol=1e-4,
                                                               verbose=False):
        #if u_l>u_r,bias>0
        self.kt = kt
        self.leadfermi = []
        for i in range(len(bias)):
            self.leadfermi.append(efermi[i] + bias[i])
        self.minfermi = min(self.leadfermi)
        self.maxfermi = max(self.leadfermi)
        self.eqinttol = eqinttol
        self.biastol = 1e-10
        
        #eq-Integral Path : 
        #default = [-100  -100+20i  10i  1e-8i] + minfermi
        #        = [-100  -100+20i  10i-20kt  5i*2*pi*kt-20kt
        #           5i*2*pi*kt+20kt] + minfermi
        
        nkt = 4 * self.kt
        dkt = 8 * np.pi * self.kt
        #self.eqintpath = [-20.0, -20.0 + dkt * 1.j, -nkt + dkt * 1.j,
        #                  dkt * 1.j + nkt]
        self.eqintpath = [ min_energy, 
                               min_energy + (10 + dkt)*1.j, 
                              #-nkt + (10 + dkt) * 1.j, 
                              -nkt + dkt * 1.j, 
                               dkt *1.j +nkt]
        self.eqdelta = dkt
        nRes = 8
        if abs( nRes - (np.round((nRes - 1) // 2) * 2 + 1)) < 1e-3 :
            print('Warning: Residue Point too close to IntPath!')
        self.eqresz = range(1, nRes, 2)
        for i in range(len(self.eqresz)):
            self.eqresz[i] *=  1.j * np.pi * self.kt
        if verbose:
            print('--eqIntCtrl: Tol = ', self.eqinttol, 'Delta =', \
                                          self.eqdelta, ' nRes =', self.eqresz)

        for i in range(len(self.eqintpath)):
            self.eqintpath[i] = self.eqintpath[i] + self.minfermi
            
        for i in range(len(self.eqresz)):
            self.eqresz[i] = self.eqresz[i] + self.minfermi        

        # read bound-Integral Path : 
        # default = [ minfermi+1e-8i   minfermi+1i
        #           maxfermi+1i   maxfermi+1e-8i ] 
        #         = [ minfermi-20kt+5i*2*pi*kt   maxfermi+20kt+5i*2*pi*kt ]
        #         = [ minfermi-20kt+5i*2*pi*kt   Mid+1i
        #            maxfermi+20kt+5i*2*pi*kt ]

        self.locinttol = self.eqinttol
        if (self.maxfermi - self.minfermi)< self.biastol:
            self.locintpath = None
            self.locdelta = 0
            self.locresz = 0
            if verbose:
                print('--locInt: None')
        else:
            nkt = 4 * self.kt
            dkt = 8 * np.pi * self.kt

            if self.maxfermi-self.minfermi < 0.2 + 2 * nkt or dkt > 0.5:
                self.locintpath = [self.minfermi - nkt + dkt * 1.j,
                                   self.maxfermi + nkt + dkt * 1.j]
            else:
                self.locintpath = [self.minfermi - nkt + dkt * 1.j,
                                   self.minfermi + nkt + dkt * 1.j,
                                   (self.maxfermi + self.minfermi) / 2 + 1.j,
                                   self.maxfermi - nkt + dkt * 1.j,
                                   self.maxfermi + nkt + dkt * 1.j]
            self.locdelta = dkt
            nRes = 8
            self.locresz = np.array(range(1, nRes, 2)
                                                   ) * 1.j * np.pi * self.kt
            tmp = len(range(1, nRes, 2))
            self.locresz = np.resize(self.locresz, [2, tmp])
            for i in range(tmp):
                self.locresz[0][i] += self.minfermi
                self.locresz[1][i] += self.maxfermi
            if verbose:
                print('--locInt: Tol =', self.locinttol, 'Delta =', \
                               self.locdelta, 'nRes=', len(self.locresz[0]))

        #ne-Integral Path : 
        # default = [ minfermi  leadfermi  maxfermi ] 
        # IntMethod = Manual Method
        # -------------------------------------------------- 
        # -- Integral Method -- 
        self.neinttol = 1e-3        
        self.neintmethod= neintmethod # 0: Linear 1: Auto

        # -- Integral Step--

        self.neintstep = neintstep                    
        
        # -- Integral Path --

        if self.maxfermi -self.minfermi < self.biastol:
            self.neintpath = []
        else :
            nkt = 4 * kt
            self.neintpath = [self.minfermi - nkt, self.maxfermi + nkt]

        # -- Integral eta --

        for i in range(len(self.neintpath)):                         
            self.neintpath[i] += 1e-2j 

        if len(self.neintpath) == 0:
            if verbose:
                print(' --neInt: None')
        elif self.neintmethod == 0:
            if verbose:
                print(' --neInt: ManualEp -> Step=', self.neintstep, 'Eta =',\
                                              np.imag(self.neintpath[0]))
        elif self.neintmethod == 1:
            if verbose:
                print(' --neInt: AutoEp   -> Tol =', self.neinttol,  'Eta =',\
                                              np.imag(self.neintpath[0]))

        

