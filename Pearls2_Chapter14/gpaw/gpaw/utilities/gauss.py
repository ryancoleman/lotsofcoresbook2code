
import numpy as np
from numpy import sqrt, pi, exp

import _gpaw
from gpaw import debug
from gpaw.utilities.tools import coordinates
from gpaw.utilities import erf, is_contiguous

# computer generated code:
# use c/bmgs/sharmonic.py::construct_gauss_code(lmax) to generate more
Y_L = [
  '0.28209479177387814',
  '0.48860251190291992 * y',
  '0.48860251190291992 * z',
  '0.48860251190291992 * x',
  '1.0925484305920792 * x*y',
  '1.0925484305920792 * y*z',
  '0.31539156525252005 * (-r2+3*z*z)',
  '1.0925484305920792 * x*z',
  '0.54627421529603959 * (-y*y+x*x)',
  '0.59004358992664352 * (-y*y*y+3*x*x*y)',
  '2.8906114426405538 * x*y*z',
  '0.45704579946446577 * (5*y*z*z-y*r2)',
  '0.3731763325901154 * (-3*z*r2+5*z*z*z)',
  '0.45704579946446577 * (-x*r2+5*x*z*z)',
  '1.4453057213202769 * (-y*y*z+x*x*z)',
  '0.59004358992664352 * (x*x*x-3*x*y*y)',
  '2.5033429417967046 * (x*x*x*y-x*y*y*y)',
  '1.7701307697799307 * (3*x*x*y*z-y*y*y*z)',
  '0.94617469575756008 * (-x*y*r2+7*x*y*z*z)',
  '0.66904654355728921 * (-3*y*z*r2+7*y*z*z*z)',
  '0.10578554691520431 * (3*r2*r2-30*z*z*r2+35*z*z*z*z)',
  '0.66904654355728921 * (7*x*z*z*z-3*x*z*r2)',
  '0.47308734787878004 * (y*y*r2+7*x*x*z*z-x*x*r2-7*y*y*z*z)',
  '1.7701307697799307 * (x*x*x*z-3*x*y*y*z)',
  '0.62583573544917614 * (-6*x*x*y*y+x*x*x*x+y*y*y*y)',
  '0.65638205684017015 * (y*y*y*y*y+5*x*x*x*x*y-10*x*x*y*y*y)',
  '8.3026492595241645 * (x*x*x*y*z-x*y*y*y*z)',
  '0.48923829943525038 * (y*y*y*r2-9*y*y*y*z*z-3*x*x*y*r2+27*x*x*y*z*z)',
  '4.7935367849733241 * (3*x*y*z*z*z-x*y*z*r2)',
  '0.45294665119569694 * (-14*y*z*z*r2+y*r2*r2+21*y*z*z*z*z)',
  '0.1169503224534236 * (63*z*z*z*z*z+15*z*r2*r2-70*z*z*z*r2)',
  '0.45294665119569694 * (x*r2*r2-14*x*z*z*r2+21*x*z*z*z*z)',
  '2.3967683924866621 * (-3*y*y*z*z*z+y*y*z*r2+3*x*x*z*z*z-x*x*z*r2)',
  '0.48923829943525038 * (9*x*x*x*z*z-27*x*y*y*z*z-x*x*x*r2+3*x*y*y*r2)',
  '2.0756623148810411 * (y*y*y*y*z-6*x*x*y*y*z+x*x*x*x*z)',
  '0.65638205684017015 * (-10*x*x*x*y*y+5*x*y*y*y*y+x*x*x*x*x)',
]
gauss_L = [
  'sqrt(a**3*4)/pi * exp(-a*r2)',
  'sqrt(a**5*5.333333333333333)/pi * y * exp(-a*r2)',
  'sqrt(a**5*5.333333333333333)/pi * z * exp(-a*r2)',
  'sqrt(a**5*5.333333333333333)/pi * x * exp(-a*r2)',
  'sqrt(a**7*4.266666666666667)/pi * x*y * exp(-a*r2)',
  'sqrt(a**7*4.266666666666667)/pi * y*z * exp(-a*r2)',
  'sqrt(a**7*0.35555555555555557)/pi * (-r2+3*z*z) * exp(-a*r2)',
  'sqrt(a**7*4.266666666666667)/pi * x*z * exp(-a*r2)',
  'sqrt(a**7*1.0666666666666667)/pi * (-y*y+x*x) * exp(-a*r2)',
  'sqrt(a**9*0.10158730158730159)/pi * (-y*y*y+3*x*x*y) * exp(-a*r2)',
  'sqrt(a**9*2.4380952380952383)/pi * x*y*z * exp(-a*r2)',
  'sqrt(a**9*0.06095238095238095)/pi * (5*y*z*z-y*r2) * exp(-a*r2)',
  'sqrt(a**9*0.040634920634920635)/pi * (-3*z*r2+5*z*z*z) * exp(-a*r2)',
  'sqrt(a**9*0.06095238095238095)/pi * (-x*r2+5*x*z*z) * exp(-a*r2)',
  'sqrt(a**9*0.6095238095238096)/pi * (-y*y*z+x*x*z) * exp(-a*r2)',
  'sqrt(a**9*0.10158730158730159)/pi * (x*x*x-3*x*y*y) * exp(-a*r2)',
  'sqrt(a**11*0.09029982363315697)/pi * (x*x*x*y-x*y*y*y) * exp(-a*r2)',
  'sqrt(a**11*0.045149911816578486)/pi * (3*x*x*y*z-y*y*y*z) * exp(-a*r2)',
  'sqrt(a**11*0.01289997480473671)/pi * (-x*y*r2+7*x*y*z*z) * exp(-a*r2)',
  'sqrt(a**11*0.006449987402368355)/pi * (-3*y*z*r2+7*y*z*z*z) * exp(-a*r2)',
  'sqrt(a**11*0.00016124968505920888)/pi * (3*r2*r2-30*z*z*r2+35*z*z*z*z) * exp(-a*r2)',
  'sqrt(a**11*0.006449987402368355)/pi * (7*x*z*z*z-3*x*z*r2) * exp(-a*r2)',
  'sqrt(a**11*0.0032249937011841773)/pi * (y*y*r2+7*x*x*z*z-x*x*r2-7*y*y*z*z) * exp(-a*r2)',
  'sqrt(a**11*0.045149911816578486)/pi * (x*x*x*z-3*x*y*y*z) * exp(-a*r2)',
  'sqrt(a**11*0.005643738977072311)/pi * (-6*x*x*y*y+x*x*x*x+y*y*y*y) * exp(-a*r2)',
  'sqrt(a**13*0.00020522687189353857)/pi * (y*y*y*y*y+5*x*x*x*x*y-10*x*x*y*y*y) * exp(-a*r2)',
  'sqrt(a**13*0.03283629950296617)/pi * (x*x*x*y*z-x*y*y*y*z) * exp(-a*r2)',
  'sqrt(a**13*0.00011401492882974365)/pi * (y*y*y*r2-9*y*y*y*z*z-3*x*x*y*r2+27*x*x*y*z*z) * exp(-a*r2)',
  'sqrt(a**13*0.01094543316765539)/pi * (3*x*y*z*z*z-x*y*z*r2) * exp(-a*r2)',
  'sqrt(a**13*9.772708185406599e-05)/pi * (-14*y*z*z*r2+y*r2*r2+21*y*z*z*z*z) * exp(-a*r2)',
  'sqrt(a**13*6.515138790271066e-06)/pi * (63*z*z*z*z*z+15*z*r2*r2-70*z*z*z*r2) * exp(-a*r2)',
  'sqrt(a**13*9.772708185406599e-05)/pi * (x*r2*r2-14*x*z*z*r2+21*x*z*z*z*z) * exp(-a*r2)',
  'sqrt(a**13*0.0027363582919138476)/pi * (-3*y*y*z*z*z+y*y*z*r2+3*x*x*z*z*z-x*x*z*r2) * exp(-a*r2)',
  'sqrt(a**13*0.00011401492882974365)/pi * (9*x*x*x*z*z-27*x*y*y*z*z-x*x*x*r2+3*x*y*y*r2) * exp(-a*r2)',
  'sqrt(a**13*0.0020522687189353855)/pi * (y*y*y*y*z-6*x*x*y*y*z+x*x*x*x*z) * exp(-a*r2)',
  'sqrt(a**13*0.00020522687189353857)/pi * (-10*x*x*x*y*y+5*x*y*y*y*y+x*x*x*x*x) * exp(-a*r2)',
]
gausspot_L = [
  '2.0*1.7724538509055159*erf(sqrt(a)*r)/r',
  '1.1547005383792515*(1.7724538509055159*erf(sqrt(a)*r)-(2)*sqrt(a)*r*exp(-a*r2))/r/r2*y',
  '1.1547005383792515*(1.7724538509055159*erf(sqrt(a)*r)-(2)*sqrt(a)*r*exp(-a*r2))/r/r2*z',
  '1.1547005383792515*(1.7724538509055159*erf(sqrt(a)*r)-(2)*sqrt(a)*r*exp(-a*r2))/r/r2*x',
  '0.5163977794943222*(5.3173615527165481*erf(sqrt(a)*r)-(6+4*(sqrt(a)*r)**2)*sqrt(a)*r*exp(-a*r2))/r/r2**2*x*y',
  '0.5163977794943222*(5.3173615527165481*erf(sqrt(a)*r)-(6+4*(sqrt(a)*r)**2)*sqrt(a)*r*exp(-a*r2))/r/r2**2*y*z',
  '0.14907119849998599*(5.3173615527165481*erf(sqrt(a)*r)-(6+4*(sqrt(a)*r)**2)*sqrt(a)*r*exp(-a*r2))/r/r2**2*(-r2+3*z*z)',
  '0.5163977794943222*(5.3173615527165481*erf(sqrt(a)*r)-(6+4*(sqrt(a)*r)**2)*sqrt(a)*r*exp(-a*r2))/r/r2**2*x*z',
  '0.2581988897471611*(5.3173615527165481*erf(sqrt(a)*r)-(6+4*(sqrt(a)*r)**2)*sqrt(a)*r*exp(-a*r2))/r/r2**2*(-y*y+x*x)',
  '0.03984095364447979*(26.586807763582737*erf(sqrt(a)*r)-(30+20*(sqrt(a)*r)**2+8*(sqrt(a)*r)**4)*sqrt(a)*r*exp(-a*r2))/r/r2**3*(-y*y*y+3*x*x*y)',
  '0.19518001458970666*(26.586807763582737*erf(sqrt(a)*r)-(30+20*(sqrt(a)*r)**2+8*(sqrt(a)*r)**4)*sqrt(a)*r*exp(-a*r2))/r/r2**3*x*y*z',
  '0.03086066999241838*(26.586807763582737*erf(sqrt(a)*r)-(30+20*(sqrt(a)*r)**2+8*(sqrt(a)*r)**4)*sqrt(a)*r*exp(-a*r2))/r/r2**3*(5*y*z*z-y*r2)',
  '0.02519763153394848*(26.586807763582737*erf(sqrt(a)*r)-(30+20*(sqrt(a)*r)**2+8*(sqrt(a)*r)**4)*sqrt(a)*r*exp(-a*r2))/r/r2**3*(-3*z*r2+5*z*z*z)',
  '0.03086066999241838*(26.586807763582737*erf(sqrt(a)*r)-(30+20*(sqrt(a)*r)**2+8*(sqrt(a)*r)**4)*sqrt(a)*r*exp(-a*r2))/r/r2**3*(-x*r2+5*x*z*z)',
  '0.09759000729485333*(26.586807763582737*erf(sqrt(a)*r)-(30+20*(sqrt(a)*r)**2+8*(sqrt(a)*r)**4)*sqrt(a)*r*exp(-a*r2))/r/r2**3*(-y*y*z+x*x*z)',
  '0.03984095364447979*(26.586807763582737*erf(sqrt(a)*r)-(30+20*(sqrt(a)*r)**2+8*(sqrt(a)*r)**4)*sqrt(a)*r*exp(-a*r2))/r/r2**3*(x*x*x-3*x*y*y)',
  '0.018781205660633703*(186.10765434507917*erf(sqrt(a)*r)-(210+140*(sqrt(a)*r)**2+56*(sqrt(a)*r)**4+16*(sqrt(a)*r)**6)*sqrt(a)*r*exp(-a*r2))/r/r2**4*(x*x*x*y-x*y*y*y)',
  '0.013280317881493264*(186.10765434507917*erf(sqrt(a)*r)-(210+140*(sqrt(a)*r)**2+56*(sqrt(a)*r)**4+16*(sqrt(a)*r)**6)*sqrt(a)*r*exp(-a*r2))/r/r2**4*(3*x*x*y*z-y*y*y*z)',
  '0.007098628499999332*(186.10765434507917*erf(sqrt(a)*r)-(210+140*(sqrt(a)*r)**2+56*(sqrt(a)*r)**4+16*(sqrt(a)*r)**6)*sqrt(a)*r*exp(-a*r2))/r/r2**4*(-x*y*r2+7*x*y*z*z)',
  '0.005019488349473618*(186.10765434507917*erf(sqrt(a)*r)-(210+140*(sqrt(a)*r)**2+56*(sqrt(a)*r)**4+16*(sqrt(a)*r)**6)*sqrt(a)*r*exp(-a*r2))/r/r2**4*(-3*y*z*r2+7*y*z*z*z)',
  '0.0007936507936507937*(186.10765434507917*erf(sqrt(a)*r)-(210+140*(sqrt(a)*r)**2+56*(sqrt(a)*r)**4+16*(sqrt(a)*r)**6)*sqrt(a)*r*exp(-a*r2))/r/r2**4*(3*r2*r2-30*z*z*r2+35*z*z*z*z)',
  '0.005019488349473618*(186.10765434507917*erf(sqrt(a)*r)-(210+140*(sqrt(a)*r)**2+56*(sqrt(a)*r)**4+16*(sqrt(a)*r)**6)*sqrt(a)*r*exp(-a*r2))/r/r2**4*(7*x*z*z*z-3*x*z*r2)',
  '0.003549314249999666*(186.10765434507917*erf(sqrt(a)*r)-(210+140*(sqrt(a)*r)**2+56*(sqrt(a)*r)**4+16*(sqrt(a)*r)**6)*sqrt(a)*r*exp(-a*r2))/r/r2**4*(y*y*r2+7*x*x*z*z-x*x*r2-7*y*y*z*z)',
  '0.013280317881493264*(186.10765434507917*erf(sqrt(a)*r)-(210+140*(sqrt(a)*r)**2+56*(sqrt(a)*r)**4+16*(sqrt(a)*r)**6)*sqrt(a)*r*exp(-a*r2))/r/r2**4*(x*x*x*z-3*x*y*y*z)',
  '0.004695301415158426*(186.10765434507917*erf(sqrt(a)*r)-(210+140*(sqrt(a)*r)**2+56*(sqrt(a)*r)**4+16*(sqrt(a)*r)**6)*sqrt(a)*r*exp(-a*r2))/r/r2**4*(-6*x*x*y*y+x*x*x*x+y*y*y*y)',
  '0.00044767942445854463*(1674.9688891057126*erf(sqrt(a)*r)-(1890+1260*(sqrt(a)*r)**2+504*(sqrt(a)*r)**4+144*(sqrt(a)*r)**6+32*(sqrt(a)*r)**8)*sqrt(a)*r*exp(-a*r2))/r/r2**5*(y*y*y*y*y+5*x*x*x*x*y-10*x*x*y*y*y)',
  '0.005662746571529173*(1674.9688891057126*erf(sqrt(a)*r)-(1890+1260*(sqrt(a)*r)**2+504*(sqrt(a)*r)**4+144*(sqrt(a)*r)**6+32*(sqrt(a)*r)**8)*sqrt(a)*r*exp(-a*r2))/r/r2**5*(x*x*x*y*z-x*y*y*y*z)',
  '0.0003336805417390959*(1674.9688891057126*erf(sqrt(a)*r)-(1890+1260*(sqrt(a)*r)**2+504*(sqrt(a)*r)**4+144*(sqrt(a)*r)**6+32*(sqrt(a)*r)**8)*sqrt(a)*r*exp(-a*r2))/r/r2**5*(y*y*y*r2-9*y*y*y*z*z-3*x*x*y*r2+27*x*x*y*z*z)',
  '0.0032693882574249982*(1674.9688891057126*erf(sqrt(a)*r)-(1890+1260*(sqrt(a)*r)**2+504*(sqrt(a)*r)**4+144*(sqrt(a)*r)**6+32*(sqrt(a)*r)**8)*sqrt(a)*r*exp(-a*r2))/r/r2**5*(3*x*y*z*z*z-x*y*z*r2)',
  '0.0003089281524450488*(1674.9688891057126*erf(sqrt(a)*r)-(1890+1260*(sqrt(a)*r)**2+504*(sqrt(a)*r)**4+144*(sqrt(a)*r)**6+32*(sqrt(a)*r)**8)*sqrt(a)*r*exp(-a*r2))/r/r2**5*(-14*y*z*z*r2+y*r2*r2+21*y*z*z*z*z)',
  '7.976490597295335e-05*(1674.9688891057126*erf(sqrt(a)*r)-(1890+1260*(sqrt(a)*r)**2+504*(sqrt(a)*r)**4+144*(sqrt(a)*r)**6+32*(sqrt(a)*r)**8)*sqrt(a)*r*exp(-a*r2))/r/r2**5*(63*z*z*z*z*z+15*z*r2*r2-70*z*z*z*r2)',
  '0.0003089281524450488*(1674.9688891057126*erf(sqrt(a)*r)-(1890+1260*(sqrt(a)*r)**2+504*(sqrt(a)*r)**4+144*(sqrt(a)*r)**6+32*(sqrt(a)*r)**8)*sqrt(a)*r*exp(-a*r2))/r/r2**5*(x*r2*r2-14*x*z*z*r2+21*x*z*z*z*z)',
  '0.0016346941287124991*(1674.9688891057126*erf(sqrt(a)*r)-(1890+1260*(sqrt(a)*r)**2+504*(sqrt(a)*r)**4+144*(sqrt(a)*r)**6+32*(sqrt(a)*r)**8)*sqrt(a)*r*exp(-a*r2))/r/r2**5*(-3*y*y*z*z*z+y*y*z*r2+3*x*x*z*z*z-x*x*z*r2)',
  '0.0003336805417390959*(1674.9688891057126*erf(sqrt(a)*r)-(1890+1260*(sqrt(a)*r)**2+504*(sqrt(a)*r)**4+144*(sqrt(a)*r)**6+32*(sqrt(a)*r)**8)*sqrt(a)*r*exp(-a*r2))/r/r2**5*(9*x*x*x*z*z-27*x*y*y*z*z-x*x*x*r2+3*x*y*y*r2)',
  '0.0014156866428822932*(1674.9688891057126*erf(sqrt(a)*r)-(1890+1260*(sqrt(a)*r)**2+504*(sqrt(a)*r)**4+144*(sqrt(a)*r)**6+32*(sqrt(a)*r)**8)*sqrt(a)*r*exp(-a*r2))/r/r2**5*(y*y*y*y*z-6*x*x*y*y*z+x*x*x*x*z)',
  '0.00044767942445854463*(1674.9688891057126*erf(sqrt(a)*r)-(1890+1260*(sqrt(a)*r)**2+504*(sqrt(a)*r)**4+144*(sqrt(a)*r)**6+32*(sqrt(a)*r)**8)*sqrt(a)*r*exp(-a*r2))/r/r2**5*(-10*x*x*x*y*y+5*x*y*y*y*y+x*x*x*x*x)',
]

def Y_L(L, x, y, z, r2):
    if L == 0:
        return 0.28209479177387814
    elif L == 1:
        return 0.48860251190291992 * y
    elif L == 2:
        return 0.48860251190291992 * z
    elif L == 3:
        return 0.48860251190291992 * x
    elif L == 4:
        return 1.0925484305920792 * x*y
    elif L == 5:
        return 1.0925484305920792 * y*z
    elif L == 6:
        return 0.31539156525252005 * (3*z*z-r2)
    elif L == 7:
        return 1.0925484305920792 * x*z
    elif L == 8:
        return 0.54627421529603959 * (x*x-y*y)

def gauss_L(a, L, x, y, z, r2, exp_ar2):
    if L == 0:
        return sqrt(a**3*4)/pi * exp_ar2
    elif L == 1:
        return sqrt(a**5*5.333333333333333)/pi * y * exp_ar2
    elif L == 2:
        return sqrt(a**5*5.333333333333333)/pi * z * exp_ar2
    elif L == 3:
        return sqrt(a**5*5.333333333333333)/pi * x * exp_ar2
    elif L == 4:
        return sqrt(a**7*4.2666666666666666)/pi * x*y * exp_ar2
    elif L == 5:
        return sqrt(a**7*4.2666666666666666)/pi * y*z * exp_ar2
    elif L == 6:
        return sqrt(a**7*0.35555555555555557)/pi * (3*z*z-r2) * exp_ar2
    elif L == 7:
        return sqrt(a**7*4.2666666666666666)/pi * x*z * exp_ar2
    elif L == 8:
        return sqrt(a**7*1.0666666666666667)/pi * (x*x-y*y) * exp_ar2

def gausspot_L(a, L, x, y, z, r, r2, erf_sar, exp_ar2):
    if L == 0:
        return 2.0*1.7724538509055159*erf_sar/r
    elif L == 1:
        return 1.1547005383792515*(1.7724538509055159*erf_sar-2*sqrt(a)*r*exp_ar2)/r/r2*y
    elif L == 2:
        return 1.1547005383792515*(1.7724538509055159*erf_sar-2*sqrt(a)*r*exp_ar2)/r/r2*z
    elif L == 3:
        return 1.1547005383792515*(1.7724538509055159*erf_sar-2*sqrt(a)*r*exp_ar2)/r/r2*x
    elif L == 4:
        return 0.5163977794943222*(5.3173615527165481*erf_sar-(6+4*(sqrt(a)*r)**2)*sqrt(a)*r*exp_ar2)/r/r2**2*x*y
    elif L == 5:
        return 0.5163977794943222*(5.3173615527165481*erf_sar-(6+4*(sqrt(a)*r)**2)*sqrt(a)*r*exp_ar2)/r/r2**2*y*z
    elif L == 6:
        return 0.14907119849998599*(5.3173615527165481*erf_sar-(6+4*(sqrt(a)*r)**2)*sqrt(a)*r*exp_ar2)/r/r2**2*(3*z*z-r2)
    elif L == 7:
        return 0.5163977794943222*(5.3173615527165481*erf_sar-(6+4*(sqrt(a)*r)**2)*sqrt(a)*r*exp_ar2)/r/r2**2*x*z
    elif L == 8:
        return 0.2581988897471611*(5.3173615527165481*erf_sar-(6+4*(sqrt(a)*r)**2)*sqrt(a)*r*exp_ar2)/r/r2**2*(x*x-y*y)

# end of computer generated code

class Gaussian:
    """Class offering several utilities related to the generalized gaussians.

    Generalized gaussians are defined by::
    
                       _____                           2  
                      /  1       l!         l+3/2  -a r   l  m
       g (x,y,z) =   / ----- --------- (4 a)      e      r  Y (x,y,z),
        L          \/  4 pi  (2l + 1)!                       l

    where a is the inverse width of the gaussian, and Y_l^m is a real
    spherical harmonic.
    The gaussians are centered in the middle of input grid-descriptor."""
    
    def __init__(self, gd, a=19., center=None):
        self.gd = gd
        self.xyz, self.r2 = coordinates(gd, center)
        self.r = np.sqrt(self.r2)
        self.set_width(a)
        self.exp_ar2 = exp(-self.a * self.r2) 
        self.erf_sar = erf(sqrt(self.a) * self.r)


    def set_width(self, a):
        """Set exponent of exp-function to -a on the boundary."""
        self.a = 4 * a * (self.gd.icell_cv**2).sum(1).max()
        
    def get_gauss(self, L):
        a = self.a
        x, y, z  = tuple(self.xyz)
        r2 = self.r2
        exp_ar2 = self.exp_ar2
        return gauss_L(a, L, x, y, z, r2, exp_ar2)
    
    def get_gauss_pot(self, L):
        a = self. a
        x, y, z  = tuple(self.xyz)
        r2 = self.r2
        r = self.r
        erf_sar = self.erf_sar
        exp_ar2 = self.exp_ar2
        return gausspot_L(a, L, x, y, z, r, r2, erf_sar, exp_ar2)

    def get_moment(self, n, L):
        r2 = self.r2
        x, y, z = tuple(self.xyz)
        return self.gd.integrate(n * Y_L(L, x, y, z, r2))

    def remove_moment(self, n, L, q=None):
        # Determine multipole moment
        if q is None:
            q = self.get_moment(n, L)

        # Don't do anything if moment is less than the tolerance
        if abs(q) < 1e-7:
            return 0.

        # Remove moment from input density
        n -= q * self.get_gauss(L)

        # Return correction
        return q * self.get_gauss_pot(L)


def gaussian_wave(r_vG, r0_v, sigma, k_v=None, A=None, dtype=float,
                  out_G=None):
    """Generates function values for atom-centered Gaussian waves.

    ::
    
                         _ _
        _            / -|r-r0|^2 \           _ _
      f(r) = A * exp( ----------- ) * exp( i k.r )
                     \ 2 sigma^2 /

    If the parameter A is not specified, the Gaussian wave is normalized::

                                                  oo
           /    ____        \ -3/2               /       _  2  2
      A = (    /    '        )        =>    4 pi | dr |f(r)|  r  = 1
           \ \/  pi   sigma /                    /
                                                   0

    Parameters:

    r_vG: ndarray
        Set of coordinates defining the grid positions.
    r0_v: ndarray
        Set of coordinates defining the center of the Gaussian envelope.
    sigma: float
        Specifies the spatial width of the Gaussian envelope.
    k_v: ndarray or None
        Set of reciprocal lattice coordinates defining the wave vector.
        An argument of None is interpreted as the gamma point i.e. k_v=0.
    A: float, complex or None
        Specifies the amplitude of the Gaussian wave. Normalizes if None.
    dtype: type, defaults to float
        Specifies the output data type. Only returns the real-part if float.
    out_G: ndarray or None
        Optional pre-allocated buffer to fill in values. Allocates if None.

    """
    if k_v is None:
        k_v = np.zeros(r0_v.shape)

    if A is None:
        # 4*pi*int(exp(-r^2/(2*sigma^2))^2 * r^2, r=0...infinity)
        # = sigma^3*pi^(3/2) = 1/A^2 -> A = (sqrt(Pi)*sigma)^(-3/2)
        A = 1/(sigma*np.pi**0.5)**1.5

    if debug:
        assert is_contiguous(r_vG, float)
        assert is_contiguous(r0_v, float)
        assert is_contiguous(k_v, float)
        assert r_vG.ndim >= 2 and r_vG.shape[0] > 0
        assert r0_v.ndim == 1 and r0_v.shape[0] > 0
        assert k_v.ndim == 1 and k_v.shape[0] > 0
        assert (r_vG.shape[0],) == r0_v.shape == k_v.shape
        assert sigma > 0

    if out_G is None:
        out_G = np.empty(r_vG.shape[1:], dtype=dtype)
    elif debug:
        assert is_contiguous(out_G)
        assert out_G.shape == r_vG.shape[1:]

    # slice_v2vG = [slice(None)] + [np.newaxis]*3
    # gw = lambda r_vG, r0_v, sigma, k_v, A=1/(sigma*np.pi**0.5)**1.5: \
    #    * np.exp(-np.sum((r_vG-r0_v[slice_v2vG])**2, axis=0)/(2*sigma**2)) \
    #    * np.exp(1j*np.sum(np.r_vG*k_v[slice_v2vG], axis=0)) * A
    _gpaw.utilities_gaussian_wave(A, r_vG, r0_v, sigma, k_v, out_G)
    return out_G

