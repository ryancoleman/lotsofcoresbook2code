/************************************************************************
 Implements Perdew, Tao, Staroverov & Scuseria 
   meta-Generalized Gradient Approximation.

  Exchange part
************************************************************************/

#include <math.h>
#include <stdlib.h>
#include <xc.h>
#include "xc_mgga.h"

typedef struct tpss_params {
  common_params common; // needs to be at the beginning of every functional_params
  XC(func_type) *x_aux;
  XC(func_type) *c_aux1;
  XC(func_type) *c_aux2;
} tpss_params;
  
/* some parameters */
static double b=0.40, c=1.59096, e=1.537, kappa=0.804, mu=0.21951;

/* This is Equation (7) from the paper and its derivatives */
static void 
x_tpss_7(double p, double alpha, 
	 double *qb, double *dqbdp, double *dqbdalpha)
{
  /* Eq. (7) */
  double a = sqrt(1.0 + b*alpha*(alpha-1.0)), h = 9.0/20.0;

  *qb    = h*(alpha - 1.0)/a + 2.0*p/3.0;
  *dqbdp = 2.0/3.0;
  *dqbdalpha = h*(1.0 + 0.5*b*(alpha-1.0))/pow(a, 3);
}


/* Equation (10) in all it's glory */
static 
void x_tpss_10(double p, double alpha,
	       double *x, double *dxdp, double *dxdalpha)
{
  double x1, dxdp1, dxdalpha1;
  double aux1, ap, apsr, p2;
  double qb, dqbdp, dqbdalpha;
  
  /* Equation 7 */
  x_tpss_7(p, alpha, &qb, &dqbdp, &dqbdalpha);

  p2   = p*p; 
  aux1 = 10.0/81.0;
  ap = (3*alpha + 5*p)*(3*alpha + 5*p);
  apsr = (3*alpha + 5*p);
  
  /* first we handle the numerator */
  x1    = 0.0;
  dxdp1 = 0.0;
  dxdalpha1 = 0.0;

  { /* first term */
    double a = (9*alpha*alpha+30*alpha*p+50*p2), a2 = a*a;
    x1    += aux1*p + 25*c*p2*p*ap/a2;
    dxdp1 += aux1 + ((3*225*c*p2*alpha*alpha+ 4*750*c*p*p2*alpha + 5*625*c*p2*p2)*a2 - 25*c*p2*p*ap*2*a*(30*alpha+50*2*p))/(a2*a2);
    dxdalpha1 += ((225*c*p*p2*2*alpha + 750*c*p2*p2)*a2 - 25*c*p2*p*ap*2*a*(9*2*alpha+30*p))/(a2*a2);
  }
  
  { /* second term */
    double a = 146.0/2025.0*qb;
    x1    += a*qb;
    dxdp1 += 2.0*a*dqbdp;
    dxdalpha1 += 2.0*a*dqbdalpha;
  }
  
  { /* third term */
    double h = 73.0/(405*sqrt(2.0));
    x1    += -h*qb*p/apsr * sqrt(ap+9);
    dxdp1 += -h * qb *((3*alpha)/ap * sqrt(ap+9) + p/apsr * 1./2. * pow(ap+9,-1./2.)* 2*apsr*5) - h*p/apsr*sqrt(ap+9)*dqbdp; 
    dxdalpha1 += -h*qb*( (-1)*p*3/ap * sqrt(ap+9) + p/apsr * 1./2. * pow(ap+9,-1./2.)* 2*apsr*3) - h*p/apsr*sqrt(ap+9)*dqbdalpha;
  }
  

  { /* forth term */
    double a = aux1*aux1/kappa;
    x1    += a*p2;
    dxdp1 += a*2.0*p;
    dxdalpha1 += 0.0;
  }
  
  { /* fifth term */
    x1    += 20*sqrt(e)*p2/(9*ap);
    dxdp1 += 20*sqrt(e)/9*(2*p*ap-p2*2*(3*alpha + 5*p)*5)/(ap*ap);
    dxdalpha1 +=-20*2*sqrt(e)/3*p2/(ap*(3*alpha + 5*p));
  }
  
  { /* sixth term */
    double a = e*mu;
    x1    += a*p*p2;
    dxdp1 += a*3.0*p2;
    dxdalpha1 += 0.0;
  }
  
  /* and now the denominator */
  {
    double a = 1.0+sqrt(e)*p, a2 = a*a;
    *x    = x1/a2;
    *dxdp = (dxdp1*a - 2.0*sqrt(e)*x1)/(a2*a);
    *dxdalpha = dxdalpha1/a2;
  }
}

static void 
x_tpss_para(XC(func_type) *lda_aux, const double *rho, const double sigma, const double tau_,
	    double *energy, double *dedd, double *vsigma, double *dedtau)
{

  double gdms, p, tau, tauw;
  double x, dxdp, dxdalpha, Fx, dFxdx;
  double tau_lsda, exunif, vxunif, dtau_lsdadd;
  double dpdd, dpdsigma;
  double alpha, dalphadd, dalphadsigma, dalphadtau; 
  double aux =  (3./10.) * pow((3*M_PI*M_PI),2./3.); 


  /* get the uniform gas energy and potential */
  const int np = 1;
  XC(lda_exc_vxc)(lda_aux, np, rho, &exunif, &vxunif);

  /* calculate |nabla rho|^2 */
  gdms = max(MIN_GRAD*MIN_GRAD, sigma);
  
  /* Eq. (4) */
  p = gdms/(4.0*pow(3*M_PI*M_PI, 2.0/3.0)*pow(rho[0], 8.0/3.0));
  dpdd = -(8.0/3.0)*p/rho[0];
  dpdsigma= 1/(4.0*pow(3*M_PI*M_PI, 2.0/3.0)*pow(rho[0], 8.0/3.0));

  /* von Weisaecker kinetic energy density */
  tauw = max(gdms/(8.0*rho[0]), 1.0e-12);
  tau = max(tau_, tauw);

  tau_lsda = aux * pow(rho[0],5./3.); 
  dtau_lsdadd = aux * 5./3.* pow(rho[0],2./3.);
  
  alpha = (tau - tauw)/tau_lsda;

  if(fabs(tauw-tau_)< 1.0e-10){
    dalphadsigma = 0.0;
    dalphadtau = 0.0;
    dalphadd = 0.0; 
  }else{
    dalphadtau = 1./tau_lsda;
    dalphadsigma = -1./(tau_lsda*8.0*rho[0]);
    dalphadd = (tauw/rho[0]* tau_lsda - (tau - tauw) * dtau_lsdadd)/ pow(tau_lsda,2.); 
  }

  /* get Eq. (10) */
  x_tpss_10(p, alpha, &x, &dxdp, &dxdalpha);

  { /* Eq. (5) */
    double a = kappa/(kappa + x);
    Fx    = 1.0 + kappa*(1.0 - a);
    dFxdx = a*a;
  }
  
  { /* Eq. (3) */

    *energy = exunif*Fx*rho[0];

    /* exunif is en per particle already so we multiply by n the terms with exunif*/

    *dedd   = vxunif*Fx + exunif*dFxdx*rho[0]*(dxdp*dpdd + dxdalpha*dalphadd);

    *vsigma = exunif*dFxdx*rho[0]*(dxdp*dpdsigma + dxdalpha*dalphadsigma);

    *dedtau = exunif*dFxdx*rho[0]*(dxdalpha*dalphadtau);


  }
}

static void 
XC(mgga_x_tpss)(void *p, const double *rho, const double *sigma, const double *tau,
                double *e, double *dedd, double *vsigma, double *dedtau)
{
  tpss_params *par = (tpss_params*)p;
  if(par->common.nspin == XC_UNPOLARIZED){
    double en;
    x_tpss_para(par->x_aux, rho, sigma[0], tau[0], &en, dedd, vsigma, dedtau);
    *e = en/(rho[0]+rho[1]);
  }else{ 
    /* The spin polarized version is handle using the exact spin scaling
       Ex[n1, n2] = (Ex[2*n1] + Ex[2*n2])/2
    */

    *e = 0.0;

    double e2na, e2nb, rhoa[2], rhob[2];

    double vsigmapart[3]; 
	  
    rhoa[0]=2*rho[0];
    rhoa[1]=0.0;
    rhob[0]=2*rho[1];
    rhob[1]=0.0;
		  
    x_tpss_para(par->x_aux, rhoa, 4*sigma[0], 2.0*tau[0], &e2na, &(dedd[0]), &(vsigmapart[0]), &(dedtau[0]));

    x_tpss_para(par->x_aux, rhob, 4*sigma[2], 2.0*tau[1], &e2nb, &(dedd[1]), &(vsigmapart[2]), &(dedtau[1]));
		 
    *e = (e2na + e2nb )/(2.*(rho[0]+rho[1]));
    vsigma[0] = 2*vsigmapart[0];
    vsigma[2] = 2*vsigmapart[2];
  }
}

/************************************************************************
 Implements Perdew, Tao, Staroverov & Scuseria 
   meta-Generalized Gradient Approximation.
   J. Chem. Phys. 120, 6898 (2004)
   http://dx.doi.org/10.1063/1.1665298

  Correlation part
************************************************************************/

/* some parameters */
static double d = 2.8;

/* Equation (14) */
static void
c_tpss_14(double csi, double zeta, double *C, double *dCdcsi, double *dCdzeta)
{
  double fz, C0, dC0dz, dfzdz;
  double z2 = zeta*zeta;
    
  /* Equation (13) */
  C0    = 0.53 + z2*(0.87 + z2*(0.50 + z2*2.26));
  dC0dz = zeta*(2.0*0.87 + z2*(4.0*0.5 + z2*6.0*2.26));  /*OK*/
  
  fz    = 0.5*(pow(1.0 + zeta, -4.0/3.0) + pow(1.0 - zeta, -4.0/3.0));
  dfzdz = 0.5*(-4.0/3.0)*(pow(1.0 + zeta, -7.0/3.0) - pow(1.0 - zeta, -7.0/3.0)); /*OK*/
  
  { /* Equation (14) */
    double csi2 = csi*csi;
    double a = 1.0 + csi2*fz, a4 = pow(a, 4);
    
    *C      =  C0 / a4;
    *dCdcsi = -8.0*C0*csi*fz/(a*a4);  /*added C OK*/
    *dCdzeta = (dC0dz*a - C0*4.0*csi2*dfzdz)/(a*a4);  /*OK*/
  }
}

/* Equation (12) */
static void c_tpss_12(XC(func_type) *aux1, XC(func_type) *aux2, int nspin, const double *rho, const double *sigma, 
                      double dens, double zeta, double z,
                      double *e_PKZB, double *de_PKZBdd, double *de_PKZBdsigma, double *de_PKZBdz)
{
  /*some incoming variables:  
    dens = rho[0] + rho[1]
    z = tau_w/tau
    zeta = (rho[0] - rho[1])/dens*/

  double e_PBE, e_PBEup, e_PBEdn;
  double de_PBEdd[2], de_PBEdsigma[3], de_PBEddup[2], de_PBEdsigmaup[3], de_PBEdddn[2], de_PBEdsigmadn[3] ;
  double aux, zsq;
  double dzetadd[2], dcsidd[2], dcsidsigma[3];  

  double C, dCdcsi, dCdzeta;
  double densp[2], densp2[2], sigmatot[3], sigmaup[3], sigmadn[3];
  int i;
  /*initialize dCdcsi and dCdzeta and the energy*/
  dCdcsi = dCdzeta = 0.0;  
  e_PBE = 0.0;
  e_PBEup = 0.0;
  e_PBEdn = 0.0;

  /* get the PBE stuff */
  if(nspin== XC_UNPOLARIZED)
    { densp[0]=rho[0]/2.;
      densp[1]=rho[0]/2.;
      sigmatot[0] = sigma[0]/4.;
      sigmatot[1] = sigma[0]/4.;
      sigmatot[2] = sigma[0]/4.;
    }else{
    densp[0] = rho[0];
    densp[1] = rho[1];
    sigmatot[0] = sigma[0];
    sigmatot[1] = sigma[1];
    sigmatot[2] = sigma[2];
  }

  /* e_PBE */
  XC(func_type) *auxfunc = (nspin == XC_UNPOLARIZED) ? aux2 : aux1;
  const int np = 1;
  XC(gga_exc_vxc)(auxfunc, np, densp, sigmatot, &e_PBE, de_PBEdd, de_PBEdsigma); 

  densp2[0]=densp[0];
  densp2[1]=0.0;

  if(nspin== XC_UNPOLARIZED)
    {
      sigmaup[0] = sigma[0]/4.;
      sigmaup[1] = 0.;
      sigmaup[2] = 0.;
    }else{
    sigmaup[0] = sigma[0];
    sigmaup[1] = 0.;
    sigmaup[2] = 0.;
  }
  /* e_PBE spin up */
  XC(gga_exc_vxc)(auxfunc, np, densp2, sigmaup, &e_PBEup, de_PBEddup, de_PBEdsigmaup); 
  
  densp2[0]=densp[1];
  densp2[1]=0.0;

  if(nspin== XC_UNPOLARIZED)
    {
      sigmadn[0] = sigma[0]/4.;
      sigmadn[1] = 0.;
      sigmadn[2] = 0.;
    }else{
    sigmadn[0] = sigma[2];
    sigmadn[1] = 0.;
    sigmadn[2] = 0.;
  }

  /* e_PBE spin down */
  XC(gga_exc_vxc)(auxfunc, np,  densp2, sigmadn, &e_PBEdn, de_PBEdddn, de_PBEdsigmadn); 
  
  /*get Eq. (13) and (14) for the polarized case*/
  if(nspin == XC_UNPOLARIZED){   
    C          = 0.53;
    dzetadd[0] = 0.0;
    dcsidd [0] = 0.0;
    dzetadd[1] = 0.0;
    dcsidd [1] = 0.0;
    for(i=0; i<3; i++) dcsidsigma[i] = 0.0;
  }else{
    // initialize derivatives
    for(i=0; i<2; i++){
      dzetadd[i] = 0.0;
      dcsidd [i] = 0.0;}

    for(i=0; i<3; i++) dcsidsigma[i] = 0.0;



    double num, gzeta, csi, a;

    /*numerator of csi: derive as grho all components and then square the 3 parts
      [2 (grho_a[0]n_b - grho_b[0]n_a) +2 (grho_a[1]n_b - grho_b[1]n_a) + 2 (grho_a[2]n_b - grho_b[2]n_a)]/(n_a+n_b)^2   
      -> 4 (sigma_aa n_b^2 - 2 sigma_ab n_a n_b + sigma_bb n_b^2)/(n_a+n_b)^2 */

    num = sigma[0] * pow(rho[1],2) - 2.* sigma[1]*rho[0]*rho[1]+ sigma[2]*pow(rho[0],2);
    num = max(num, 1e-20);
    gzeta = sqrt(4*(num))/(dens*dens);
    gzeta = max(gzeta, MIN_GRAD);
    /*denominator of csi*/
    a = 2*pow(3.0*M_PI*M_PI*dens, 1.0/3.0);

    csi = gzeta/a;

    c_tpss_14(csi, zeta, &C, &dCdcsi, &dCdzeta);

    dzetadd[0] =  (1.0 - zeta)/dens; /*OK*/
    dzetadd[1] = -(1.0 + zeta)/dens; /*OK*/


    dcsidd [0] = 0.5*csi*(-2*sigma[1]*rho[1]+2*sigma[2]*rho[0])/num - 7./3.*csi/dens; /*OK*/
    dcsidd [1] = 0.5*csi*(-2*sigma[1]*rho[0]+2*sigma[0]*rho[1])/num - 7./3.*csi/dens; /*OK*/

    dcsidsigma[0]=  csi*pow(rho[1],2)/(2*num);   /*OK*/
    dcsidsigma[1]= -csi*rho[0]*rho[1]/num;  /*OK*/
    dcsidsigma[2]=  csi*pow(rho[0],2)/(2*num);   /*OK*/

  }

  aux = (densp[0] * max(e_PBEup, e_PBE) + densp[1] * max(e_PBEdn, e_PBE)) / dens;

  double dauxdd[2], dauxdsigma[3];
      
  if(e_PBEup > e_PBE)
    {
      //case densp[0] * e_PBEup
      dauxdd[0] = de_PBEddup[0];
      dauxdd[1] = 0.0;
      dauxdsigma[0] = de_PBEdsigmaup[0];
      dauxdsigma[1] = 0.0;
      dauxdsigma[2] = 0.0;
    }else{
    //case densp[0] * e_PBE
    dauxdd[0] = densp[0] / dens * (de_PBEdd[0] - e_PBE) + e_PBE;
    dauxdd[1] = densp[0] / dens * (de_PBEdd[1] - e_PBE);
    dauxdsigma[0] = densp[0] / dens * de_PBEdsigma[0];
    dauxdsigma[1] = densp[0] / dens * de_PBEdsigma[1];
    dauxdsigma[2] = densp[0] / dens * de_PBEdsigma[2];
  }

  if(e_PBEdn > e_PBE)
    {//case densp[1] * e_PBEdn
      dauxdd[0] += 0.0;
      dauxdd[1] += de_PBEdddn[0];
      dauxdsigma[0] += 0.0;
      dauxdsigma[1] += 0.0;
      dauxdsigma[2] += de_PBEdsigmadn[0];
    }else{//case densp[1] * e_PBE
    dauxdd[0] += densp[1] / dens * (de_PBEdd[0] - e_PBE);
    dauxdd[1] += densp[1] / dens * (de_PBEdd[1] - e_PBE) + e_PBE;
    dauxdsigma[0] += densp[1] / dens * de_PBEdsigma[0];
    dauxdsigma[1] += densp[1] / dens * de_PBEdsigma[1];
    dauxdsigma[2] += densp[1] / dens * de_PBEdsigma[2];
  }
 
  zsq=z*z;
  *e_PKZB    = (e_PBE*(1.0 + C * zsq) - (1.0 + C) * zsq * aux);
  *de_PKZBdz = dens * e_PBE * C * 2*z - dens * (1.0 + C) * 2*z * aux;  /*? think ok*/

      
  double dCdd[2];
      
  dCdd[0] = dCdzeta*dzetadd[0] + dCdcsi*dcsidd[0]; /*OK*/
  dCdd[1] = dCdzeta*dzetadd[1] + dCdcsi*dcsidd[1]; /*OK*/
      
  /* partial derivatives*/
  de_PKZBdd[0] = de_PBEdd[0] * (1.0 + C*zsq) + dens * e_PBE * dCdd[0] * zsq
    - zsq * (dens*dCdd[0] * aux + (1.0 + C) * dauxdd[0]);
  de_PKZBdd[1] = de_PBEdd[1] * (1.0 + C*zsq) + dens * e_PBE * dCdd[1] * zsq
    - zsq * (dens*dCdd[1] * aux + (1.0 + C) * dauxdd[1]);
			  
  int nder = (nspin==XC_UNPOLARIZED) ? 1 : 3;
  for(i=0; i<nder; i++){
    if(nspin==XC_UNPOLARIZED) dauxdsigma[i] /= 2.;
    double dCdsigma[i]; 
    dCdsigma[i]=  dCdcsi*dcsidsigma[i];
	
    /* partial derivatives*/
    de_PKZBdsigma[i] = de_PBEdsigma[i] * (1.0 + C * zsq) + dens * e_PBE * dCdsigma[i] * zsq
      - zsq * (dens * dCdsigma[i] * aux + (1.0 + C) * dauxdsigma[i]);

  }
} 


static void 
XC(mgga_c_tpss)(void* p, const double *rho, const double *sigma, const double *tau,
                double *energy, double *dedd, double *vsigma, double *dedtau)
{
  tpss_params *par = (tpss_params*)p;
  double dens, zeta, grad;
  double tautr, taut, tauw, z;
  double e_PKZB, de_PKZBdd[2], de_PKZBdsigma[3], de_PKZBdz;
  int i, is;
  int nspin = par->common.nspin;

  zeta = (rho[0]-rho[1])/(rho[0]+rho[1]);

  dens = rho[0];
  tautr = tau[0];
  grad  = sigma[0];

  if(nspin == XC_POLARIZED) {
    dens  += rho[1];
    tautr += tau[1];
    grad  += (2*sigma[1] + sigma[2]);
  }

  grad = max(MIN_GRAD*MIN_GRAD, grad);
  tauw = max(grad/(8.0*dens), 1.0e-12);

  taut = max(tautr, tauw);

  z = tauw/taut;
  
  double sigmatmp[3];
  sigmatmp[0] = max(MIN_GRAD*MIN_GRAD, sigma[0]);
  sigmatmp[1] = 0.0;
  sigmatmp[2] = 0.0;
  if(nspin == XC_POLARIZED) 
    {
      //sigma[1] = max(MIN_GRAD*MIN_GRAD, sigma[1]);
      sigmatmp[1] = sigma[1];
      sigmatmp[2] = max(MIN_GRAD*MIN_GRAD, sigma[2]);
    }

  /* Equation (12) */
  c_tpss_12(par->c_aux1, par->c_aux2, nspin, rho, sigmatmp, dens, zeta, z,
	    &e_PKZB, de_PKZBdd, de_PKZBdsigma, &de_PKZBdz);

  /* Equation (11) */
  {
    double z2 = z*z, z3 = z2*z;
    double dedz;
    double dzdd[2], dzdsigma[3], dzdtau;

    if(tauw >= tautr || fabs(tauw- tautr)< 1.0e-10){ 
      dzdtau = 0.0;              
      dzdd[0] = 0.0;                
      dzdd[1] = 0.0;                
      dzdsigma[0] = 0.0;
      dzdsigma[1] = 0.0;
      dzdsigma[2] = 0.0;
    }else{
      dzdtau = -z/taut;              
      dzdd[0] = - z/dens;
      dzdd[1] = 0.0;
      if (nspin == XC_POLARIZED) dzdd[1] = - z/dens;
      dzdsigma[0] = 1.0/(8*dens*taut);    
      dzdsigma[1] = 0.0;  
      dzdsigma[2] = 0.0;
      if (nspin == XC_POLARIZED) {
        dzdsigma[1] = 2.0/(8*dens*taut);    
        dzdsigma[2] = 1.0/(8*dens*taut);    
      }
    }
    
    *energy = e_PKZB * (1.0 + d*e_PKZB*z3);
    /* due to the definition of na and nb in libxc.c we need to divide by (na+nb) to recover the 
     * same energy for polarized and unpolarized calculation with the same total density */
    if(nspin == XC_UNPOLARIZED) *energy *= dens/(rho[0]+rho[1]);
	
    dedz = de_PKZBdz*(1.0 + 2.0*d*e_PKZB*z3) +  dens*e_PKZB * e_PKZB * d * 3.0*z2;  

    for(is=0; is<nspin; is++){
      dedd[is]   = de_PKZBdd[is] * (1.0 + 2.0*d*e_PKZB*z3) + dedz*dzdd[is] - e_PKZB*e_PKZB * d * z3; /*OK*/
      dedtau[is] = dedz * dzdtau; /*OK*/
    }
    int nder = (nspin==XC_UNPOLARIZED) ? 1 : 3;
    for(i=0; i<nder; i++){  
      vsigma[i] = de_PKZBdsigma[i] * (1.0 + 2.0*d*e_PKZB*z3) + dedz*dzdsigma[i];
    }
  }
}

static void tpss_init(void *p) {
  tpss_params *par = (tpss_params*)p;
  par->x_aux = (XC(func_type) *) malloc(sizeof(XC(func_type)));
  XC(func_init)(par->x_aux, XC_LDA_X, XC_UNPOLARIZED);

  par->c_aux1 = (XC(func_type) *) malloc(sizeof(XC(func_type)));
  par->c_aux2 = (XC(func_type) *) malloc(sizeof(XC(func_type)));
  XC(func_init)(par->c_aux1, XC_GGA_C_PBE, par->common.nspin);
  XC(func_init)(par->c_aux2, XC_GGA_C_PBE, XC_POLARIZED);
}

static void tpss_end(void *p) {
  tpss_params *par = (tpss_params*)p;
  XC(func_end)(par->x_aux);
  free(par->x_aux);

  XC(func_end)(par->c_aux1);
  XC(func_end)(par->c_aux2);
  free(par->c_aux1);
  free(par->c_aux2);
}

const mgga_func_info tpss_info = {
  sizeof(tpss_params),
  &tpss_init,
  &tpss_end,
  &XC(mgga_x_tpss),
  &XC(mgga_c_tpss)
};
