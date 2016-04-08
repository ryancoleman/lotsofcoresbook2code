
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <xc.h>
#include "xc_mgga.h"

/************************************************************************
 Implements Perdew, Burke & Ernzerhof Generalized Gradient Approximation
 correlation functional.

 I based this implementation on a routine from L.C. Balbas and J.M. Soler
************************************************************************/

// from old libxc util.h

#define RS(x)          (pow((3.0/(4*M_PI*x)), 1.0/3.0))
typedef struct XC(perdew_t) {
  int    nspin;
  double dens, zeta, gdmt;
  double ecunif, vcunif[2], fcunif[3];

  double  rs,  kf,  ks,  phi, t;
  double drs, dkf, dks, dphi, dt, decunif;

  double d2rs2, d2rskf, d2rsks, d2rsphi,  d2rst,  d2rsecunif;
  double        d2kf2,  d2kfks, d2kfphi,  d2kft,  d2kfecunif;
  double                 d2ks2, d2ksphi,  d2kst,  d2ksecunif;
  double                         d2phi2, d2phit, d2phiecunif;
  double                                   d2t2,   d2tecunif;
  double                                           d2ecunif2;
} XC(perdew_t);

// from old libxc util.c

/* this function converts the spin-density into total density and
	 relative magnetization */
inline void
XC(rho2dzeta)(int nspin, const double *rho, double *d, double *zeta)
{
  assert(nspin==XC_UNPOLARIZED || nspin==XC_POLARIZED);
  
  if(nspin==XC_UNPOLARIZED){
    *d    = max(MIN_DENS, rho[0]);
    *zeta = 0.0;
  }else{
    *d    = max(MIN_DENS, rho[0]+rho[1]);
    *zeta = (*d > MIN_DENS) ? (rho[0]-rho[1])/(*d) : 0.0;
  }
}

// from old libxc gga_perdew.c

static void 
XC(perdew_params)(const XC(func_type) *gga_p, const double *rho, const double *sigma, int order, XC(perdew_t) *pt)
{
  pt->nspin = gga_p->nspin;
  XC(rho2dzeta)(pt->nspin, rho, &(pt->dens), &(pt->zeta));

  const int np = 1;
  switch (order){
  case 0:
    XC(lda_exc) (gga_p, np, rho, &(pt->ecunif));
    break;
  case 1:
    XC(lda_exc_vxc)(gga_p, np, rho, &(pt->ecunif), pt->vcunif);
    break;
  case 2:
    XC(lda)(gga_p, np, rho, &(pt->ecunif), pt->vcunif, pt->fcunif, NULL);
    break;
  }

  pt->rs = RS(pt->dens);
  pt->kf = pow(3.0*M_PI*M_PI*pt->dens, 1.0/3.0);
  pt->ks = sqrt(4.0*pt->kf/M_PI);

  /* phi is bounded between 2^(-1/3) and 1 */
  pt->phi  = 0.5*(pow(1.0 + pt->zeta, 2.0/3.0) + pow(1.0 - pt->zeta, 2.0/3.0));

  /* get gdmt = |nabla n| */
  pt->gdmt = sigma[0];
  if(pt->nspin == XC_POLARIZED) pt->gdmt += 2.0*sigma[1] + sigma[2];
  if(pt->gdmt < MIN_GRAD*MIN_GRAD) pt->gdmt = MIN_GRAD*MIN_GRAD;
  pt->gdmt = sqrt(pt->gdmt);

  pt->t = pt->gdmt/(2.0 * pt->phi * pt->ks * pt->dens);

  if(order > 0)
    pt->drs = pt->dkf = pt->dks = pt->dphi = pt->dt = pt->decunif = 0.0;

  if(order > 1){
    pt->d2rs2 = pt->d2rskf = pt->d2rsks = pt->d2rsphi = pt->d2rst  = pt->d2rsecunif  = 0.0;
                pt->d2kf2  = pt->d2kfks = pt->d2kfphi = pt->d2kft  = pt->d2kfecunif  = 0.0;
		             pt->d2ks2  = pt->d2ksphi = pt->d2kst  = pt->d2ksecunif  = 0.0;
                                          pt->d2phi2  = pt->d2phit = pt->d2phiecunif = 0.0;
                                                        pt->d2t2   = pt->d2tecunif   = 0.0;
                                                                     pt->d2ecunif2   = 0.0;
  }
}

static void 
XC(perdew_potentials)(XC(perdew_t) *pt, const double *rho, double e_gga, int order,
		      double *vrho, double *vsigma,
		      double *v2rho2, double *v2rhosigma, double *v2sigma2)
{
  /* alpha = {0->rs, 1->kf, 2->ks, 3->phi, 4->t, 5->ec */
  double   dalphadd[6][2],   dFdalpha[6];
  double d2alphadd2[6][3], d2Fdalpha2[6][6];

  double dzdd[2], dpdz, d2zdd2[3], d2pdz2;
  double dtdsig, d2tdsig2;
  int is, js, ks, ns;
 
  if(order < 1) return;

  if(pt->nspin == XC_POLARIZED){
    dpdz    = 0.0;
    if(fabs(1.0 + pt->zeta) >= MIN_DENS)
      dpdz += 1.0/(3.0*pow(1.0 + pt->zeta, 1.0/3.0));
    if(fabs(1.0 - pt->zeta) >= MIN_DENS)
      dpdz -= 1.0/(3.0*pow(1.0 - pt->zeta, 1.0/3.0));

    dzdd[0] =  (1.0 - pt->zeta)/pt->dens;
    dzdd[1] = -(1.0 + pt->zeta)/pt->dens;
  }else{
    dpdz    = 0.0;
    dzdd[0] = 0.0;
  }

  dFdalpha[0] = pt->drs;
  dFdalpha[1] = pt->dkf;
  dFdalpha[2] = pt->dks;
  dFdalpha[3] = pt->dphi;
  dFdalpha[4] = pt->dt;
  dFdalpha[5] = pt->decunif;

  for(is=0; is<pt->nspin; is++){
    dalphadd[0][is] = -pt->rs/(3.0*pt->dens);
    dalphadd[1][is] =  pt->kf/(3.0*pt->dens);
    dalphadd[2][is] =  pt->ks*dalphadd[1][is]/(2.0*pt->kf);
    dalphadd[3][is] =  dpdz*dzdd[is];
    dalphadd[4][is] = -pt->t*(1.0/pt->dens + dalphadd[2][is]/pt->ks + dalphadd[3][is]/pt->phi);;
    dalphadd[5][is] = (pt->vcunif[is] - pt->ecunif)/pt->dens;
  }

  /* calculate vrho */
  if(vrho != NULL)
    for(is=0; is<pt->nspin; is++){
      if(rho[is] > MIN_DENS){
	int k;
	
	vrho[is] = e_gga;
	for(k=0; k<6; k++)
	  vrho[is] += pt->dens * dFdalpha[k]*dalphadd[k][is];
      }else{
      vrho[is] = 0.0;
      }
    }

  dtdsig  = pt->t/(2.0*pt->gdmt*pt->gdmt);

  if(vrho != NULL){ /* calculate now vsigma */
    vsigma[0] = pt->dens*pt->dt*dtdsig;
    if(pt->nspin == XC_POLARIZED){
      vsigma[1] = 2.0*vsigma[0];
      vsigma[2] =     vsigma[0];
    }
  }

  if(order < 2) return;

  /* first let us sort d2Fdalpha2 in a matrix format */
  d2Fdalpha2[0][0] = pt->d2rs2;
  d2Fdalpha2[0][1] = pt->d2rskf;
  d2Fdalpha2[0][2] = pt->d2rsks;
  d2Fdalpha2[0][3] = pt->d2rst;
  d2Fdalpha2[0][4] = pt->d2rsphi;
  d2Fdalpha2[0][5] = pt->d2rsecunif;

  d2Fdalpha2[1][0] = d2Fdalpha2[0][1];
  d2Fdalpha2[1][1] = pt->d2kf2;
  d2Fdalpha2[1][2] = pt->d2kfks;
  d2Fdalpha2[1][3] = pt->d2kft;
  d2Fdalpha2[1][4] = pt->d2kfphi;
  d2Fdalpha2[1][5] = pt->d2kfecunif;

  d2Fdalpha2[2][0] = d2Fdalpha2[0][2];
  d2Fdalpha2[2][1] = d2Fdalpha2[1][2];
  d2Fdalpha2[2][2] = pt->d2ks2;
  d2Fdalpha2[2][3] = pt->d2kst;
  d2Fdalpha2[2][4] = pt->d2ksphi;
  d2Fdalpha2[2][5] = pt->d2ksecunif;

  d2Fdalpha2[3][0] = d2Fdalpha2[0][3];
  d2Fdalpha2[3][1] = d2Fdalpha2[1][3];
  d2Fdalpha2[3][2] = d2Fdalpha2[2][3];
  d2Fdalpha2[3][3] = pt->d2phi2;
  d2Fdalpha2[3][4] = pt->d2phit;
  d2Fdalpha2[3][5] = pt->d2phiecunif;

  d2Fdalpha2[4][0] = d2Fdalpha2[0][4];
  d2Fdalpha2[4][1] = d2Fdalpha2[1][4];
  d2Fdalpha2[4][2] = d2Fdalpha2[2][4];
  d2Fdalpha2[4][3] = d2Fdalpha2[3][4];
  d2Fdalpha2[4][4] = pt->d2t2;
  d2Fdalpha2[4][5] = pt->d2tecunif;

  d2Fdalpha2[5][0] = d2Fdalpha2[0][5];
  d2Fdalpha2[5][1] = d2Fdalpha2[1][5];
  d2Fdalpha2[5][2] = d2Fdalpha2[2][5];
  d2Fdalpha2[5][3] = d2Fdalpha2[3][5];
  d2Fdalpha2[5][4] = d2Fdalpha2[4][5];
  d2Fdalpha2[5][5] = pt->d2ecunif2;

  /* now we sort d2alphadd2 */
  if(pt->nspin == XC_POLARIZED){
    d2pdz2 = 0.0;
    if(fabs(1.0 + pt->zeta) >= MIN_DENS)
      d2pdz2 += -(1.0/9.0)*pow(1.0 + pt->zeta, -4.0/3.0);
    if(fabs(1.0 - pt->zeta) >= MIN_DENS)
      d2pdz2 += -(1.0/9.0)*pow(1.0 - pt->zeta, -4.0/3.0);

    d2zdd2[0] = -2.0*dzdd[0]/pt->dens;
    d2zdd2[1] =  2.0*pt->zeta/(pt->dens*pt->dens);
    d2zdd2[2] = -2.0*dzdd[1]/pt->dens;
  }else{
    d2pdz2    = 0.0;
    d2zdd2[0] = 0.0;
  }

  ns = (pt->nspin == XC_UNPOLARIZED) ? 0 : 2;
  for(ks=0; ks<=ns; ks++){
    is = (ks == 0 || ks == 1) ? 0 : 1;
    js = (ks == 0           ) ? 0 : 1;

    d2alphadd2[0][ks] =  4.0/9.0*pt->rs/(pt->dens*pt->dens);

    d2alphadd2[1][ks] = -2.0/9.0*pt->kf/(pt->dens*pt->dens);

    d2alphadd2[2][ks] =  pt->ks/(2.0*pt->kf)*
      (d2alphadd2[1][ks] - dalphadd[1][is]*dalphadd[1][js]/(2.0*pt->kf));

    d2alphadd2[3][ks] =  d2pdz2*dzdd[is]*dzdd[js] + dpdz*d2zdd2[ks];

    d2alphadd2[4][ks] =  pt->t *
      (+2.0/(pt->dens*pt->dens)
       +2.0/(pt->ks*pt->ks)   *(dalphadd[2][is] * dalphadd[2][js])
       +2.0/(pt->phi*pt->phi) *(dalphadd[3][is] * dalphadd[3][js])
       +1.0/(pt->dens*pt->ks) *(dalphadd[2][is] + dalphadd[2][js])
       +1.0/(pt->dens*pt->phi)*(dalphadd[3][is] + dalphadd[3][js])
       +1.0/(pt->ks*pt->phi)  *(dalphadd[2][is]*dalphadd[3][js] + dalphadd[2][js]*dalphadd[3][is])
       -1.0/(pt->ks)*d2alphadd2[2][ks] -1.0/(pt->phi)*d2alphadd2[3][ks]);

    d2alphadd2[5][ks] = pt->fcunif[ks]/pt->dens -
      (pt->vcunif[is] + pt->vcunif[js] - 2.0*pt->ecunif)/(pt->dens*pt->dens);
  }

  for(ks=0; ks<=ns; ks++){
    int j, k;

    is = (ks == 0 || ks == 1) ? 0 : 1;
    js = (ks == 0           ) ? 0 : 1;

    v2rho2[ks] = 0.0;

    for(j=0; j<6; j++){
      v2rho2[ks] += dFdalpha[j]*(dalphadd[j][is] + dalphadd[j][js]);
      v2rho2[ks] += pt->dens * dFdalpha[j]*d2alphadd2[j][ks];

      for(k=0; k<6; k++)
	v2rho2[ks] +=  pt->dens * d2Fdalpha2[j][k]*dalphadd[j][is]*dalphadd[k][js];
    }
  }

  /* now we handle v2rhosigma */
  for(is=0; is<pt->nspin; is++){
    int j;
    ks = (is == 0) ? 0 : 5;

    v2rhosigma[ks] = dFdalpha[4]*dtdsig;

    for(j=0; j<6; j++)
      v2rhosigma[ks] += pt->dens * d2Fdalpha2[4][j]*dalphadd[j][is]*dtdsig;

    v2rhosigma[ks] += pt->dens * dFdalpha[4]*dalphadd[4][is]/(2.0*pt->gdmt*pt->gdmt);
  }

  if(pt->nspin == XC_POLARIZED){
    v2rhosigma[1] = 2.0*v2rhosigma[0];
    v2rhosigma[2] =     v2rhosigma[0];
    v2rhosigma[3] =     v2rhosigma[5];
    v2rhosigma[4] = 2.0*v2rhosigma[5];
  }

  /* now wwe take care of v2sigma2 */
  d2tdsig2 = -dtdsig/(2.0*pt->gdmt*pt->gdmt);
  v2sigma2[0] = pt->dens*(pt->d2t2*dtdsig*dtdsig + pt->dt*d2tdsig2);
  if(pt->nspin == XC_POLARIZED){
    v2sigma2[1] = 2.0*v2sigma2[0]; /* aa_ab */
    v2sigma2[2] =     v2sigma2[0]; /* aa_bb */
    v2sigma2[3] = 4.0*v2sigma2[0]; /* ab_ab */
    v2sigma2[4] = 2.0*v2sigma2[0]; /* ab_bb */
    v2sigma2[5] =     v2sigma2[0]; /* bb_bb */
  }
  
}

// from old libxc gga_c_pbe.c

static const double beta[4]  = {
  0.06672455060314922,  /* original PBE */
  0.046,                /* PBE sol      */
  0.089809,
  0.06672455060314922  /* PBE for revTPSS */
};

static double gamm[4];

static inline void 
pbe_eq8(int func, int order, double rs, double ecunif, double phi, 
	double *A, double *dec, double *dphi, double *drs,
	double *dec2, double *decphi, double *dphi2)
{
  double phi3, f1, df1dphi, d2f1dphi2, f2, f3, dx, d2x;

  phi3 = pow(phi, 3);
  f1   = ecunif/(gamm[func]*phi3);
  f2   = exp(-f1);
  f3   = f2 - 1.0;

  *A   = beta[func]/(gamm[func]*f3);
  if(func == 3) *A *= (1. + 0.1*rs)/(1. + 0.1778*rs);

  if(order < 1) return;

  df1dphi = -3.0*f1/phi;
  dx      = (*A)*f2/f3;

  *dec    = dx/(gamm[func]*phi3);
  *dphi   = dx*df1dphi;
  *drs    = 0.0;
  if(func == 3) *drs = beta[func]*((0.1-0.1778)/pow(1+0.1778*rs,2))/(gamm[func]*f3);

  if(func ==3) return;
  if(order < 2) return;

  d2f1dphi2 = -4.0*df1dphi/phi;
  d2x       = dx*(2.0*f2 - f3)/f3;
  *dphi2    = d2x*df1dphi*df1dphi + dx*d2f1dphi2;
  *decphi   = (d2x*df1dphi*f1 + dx*df1dphi)/ecunif;
  *dec2     = d2x/(gamm[func]*gamm[func]*phi3*phi3);
}


static void 
pbe_eq7(int func, int order, double rs, double phi, double t, double A, 
	double *H, double *dphi, double *drs, double *dt, double *dA,
	double *d2phi, double *d2phit, double *d2phiA, double *d2t2, double *d2tA, double *d2A2)
{
  double t2, phi3, f1, f2, f3;
  double df1dt, df2drs, df2dt, df1dA, df2dA;
  double d2f1dt2, d2f2dt2, d2f2dA2, d2f1dtA, d2f2dtA;

  t2   = t*t;
  phi3 = pow(phi, 3);

  f1 = t2 + A*t2*t2;
  f3 = 1.0 + A*f1;
  f2 = beta[func]*f1/(gamm[func]*f3);
  if(func == 3) f2 *= (1. + 0.1*rs)/(1. + 0.1778*rs);

  *H = gamm[func]*phi3*log(1.0 + f2);

  if(order < 1) return;

  *dphi  = 3.0*(*H)/phi;
    
  df1dt  = t*(2.0 + 4.0*A*t2);
  df2dt  = beta[func]/(gamm[func]*f3*f3) * df1dt;
  if(func == 3) df2dt*=(1. + 0.1*rs)/(1. + 0.1778*rs);
  *dt    = gamm[func]*phi3*df2dt/(1.0 + f2);
    
  df1dA  = t2*t2;
  df2dA  = beta[func]/(gamm[func]*f3*f3) * (df1dA - f1*f1);
  if(func == 3) df2dA *= (1. + 0.1*rs)/(1. + 0.1778*rs);
  *dA    = gamm[func]*phi3*df2dA/(1.0 + f2);

  df2drs = 0.0;
  *drs = 0.0;
  if(func == 3){
    df2drs = beta[func]*((0.1-0.1778)/pow(1+0.1778*rs,2))*f1/(gamm[func]*f3);
    *drs = gamm[func]*phi3*df2drs/(1.0 + f2);
  }

  if(func ==3) return;
  if(order < 2) return;

  *d2phi  = 2.0*(*dphi)/phi;
  *d2phit = 3.0*(*dt)/phi;
  *d2phiA = 3.0*(*dA)/phi;

  d2f1dt2 = 2.0 + 4.0*3.0*A*t2;
  d2f2dt2 = beta[func]/(gamm[func]*f3*f3) * (d2f1dt2 - 2.0*A/f3*df1dt*df1dt);
  *d2t2   = gamm[func]*phi3*(d2f2dt2*(1.0 + f2) - df2dt*df2dt)/((1.0 + f2)*(1.0 + f2));

  d2f1dtA = 4.0*t*t2;
  d2f2dtA = beta[func]/(gamm[func]*f3*f3) * 
    (d2f1dtA - 2.0*df1dt*(f1 + A*df1dA)/f3);
  *d2tA   = gamm[func]*phi3*(d2f2dtA*(1.0 + f2) - df2dt*df2dA)/((1.0 + f2)*(1.0 + f2));

  d2f2dA2 = beta[func]/(gamm[func]*f3*f3*f3) *(-2.0)*(2.0*f1*df1dA - f1*f1*f1 + A*df1dA*df1dA);
  *d2A2   = gamm[func]*phi3*(d2f2dA2*(1.0 + f2) - df2dA*df2dA)/((1.0 + f2)*(1.0 + f2));
}

void 
gga_c_pbe_revtpss(XC(func_type) *p, const double *rho, const double *sigma,
                  double *e, double *vrho, double *vsigma,
                  double *v2rho2, double *v2rhosigma, double *v2sigma2)
{
  gamm[0] = gamm[1] = gamm[3] = (1.0 - log(2.0))/(M_PI*M_PI);

  XC(perdew_t) pt;

  int func, order;
  double me;
  double A, dAdec, dAdphi, dAdrs, d2Adec2, d2Adecphi, d2Adphi2;
  double H, dHdphi, dHdrs, dHdt, dHdA, d2Hdphi2, d2Hdphit, d2HdphiA, d2Hdt2, d2HdtA, d2HdA2;

  d2HdphiA = 0.0;
  d2Hdphi2 = 0.0;
  d2Adphi2 = 0.0;
  d2HdA2 = 0.0;
  d2HdtA = 0.0;
  d2Hdphit = 0.0;
  d2Adecphi = 0.0;
  d2Hdt2 = 0.0;
  d2Adec2 = 0.0;
  dAdrs = 0.0;
  dAdphi = 0.0;
  dAdec = 0.0;
  dHdA = 0.0;
  dHdt = 0.0;
  dHdrs = 0.0;
  dHdphi = 0.0;

  func = 3; // for revTPSS

  order = 0;
  if(vrho   != NULL) order = 1;
  if(v2rho2 != NULL) order = 2;

  XC(perdew_params)(p, rho, sigma, order, &pt);


  pbe_eq8(func, order, pt.rs, pt.ecunif, pt.phi,
	  &A, &dAdec, &dAdphi, &dAdrs, &d2Adec2, &d2Adecphi, &d2Adphi2);

  pbe_eq7(func, order, pt.rs, pt.phi, pt.t, A, 
	  &H, &dHdphi, &dHdrs, &dHdt, &dHdA, &d2Hdphi2, &d2Hdphit, &d2HdphiA, &d2Hdt2, &d2HdtA, &d2HdA2);

  me = pt.ecunif + H;
  if(e != NULL) *e = me;

  if(order >= 1){
    pt.dphi    = dHdphi + dHdA*dAdphi;
    pt.drs     = dHdrs + dHdA*dAdrs;
    pt.dt      = dHdt;
    pt.decunif = 1.0 + dHdA*dAdec;
  }

  if(order >= 2){
    pt.d2phi2      = d2Hdphi2 + 2.0*d2HdphiA*dAdphi + dHdA*d2Adphi2 + d2HdA2*dAdphi*dAdphi;
    pt.d2phit      = d2Hdphit + d2HdtA*dAdphi;
    pt.d2phiecunif = d2HdphiA*dAdec + d2HdA2*dAdphi*dAdec + dHdA*d2Adecphi;

    pt.d2t2        = d2Hdt2;
    pt.d2tecunif   = d2HdtA*dAdec;

    pt.d2ecunif2   = d2HdA2*dAdec*dAdec + dHdA*d2Adec2;
  }

  XC(perdew_potentials)(&pt, rho, me, order, vrho, vsigma, v2rho2, v2rhosigma, v2sigma2);
}
