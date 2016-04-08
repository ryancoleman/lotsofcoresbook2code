/************************************************************************
 Implements Zhao, Truhlar
   Meta-gga M06-Local

  Correlation part
************************************************************************/

#include <math.h>
#include <stdlib.h>
#include <xc.h>
#include "xc_mgga.h"

typedef struct m06l_params {
  common_params common; // needs to be at the beginning of every functional_params
  XC(func_type) *c_aux;
  XC(func_type) *x_aux;
} m06l_params;
  
/* derivatives of x and z with respect to rho, grho and tau*/
static void 
c_m06l_zx(double x, double z, double rho, double tau, double *dxdd, double *dxdgd, double *dzdd, double *dzdtau)
{
  *dxdd = -8./3. * x * 1/rho;
  *dxdgd = 1./pow(rho,8./3.);

  *dzdd = -5./3. * 2 * tau/pow(rho, 8./3.);
  *dzdtau = 2./pow(rho, 5./3.);
}

/* Get g for Eq. (13)*/
static void
c_m06_13(double *x, double *rho, double *g_ab, double *dg_abdd, double *dg_abdgd)
{
  /*define the C_ab,i */
  static double c_ab0= 0.6042374, c_ab1= 177.6783, c_ab2= -251.3252, c_ab3=76.35173, c_ab4=-12.55699;
  double gammaCab = 0.0031 ;
  double x_ab, a; 
  double dg_abdx, dxdd_a, dxdgd_a, dzdd_a, dzdtau_a;
  double dxdd_b, dxdgd_b, dzdd_b, dzdtau_b;

  /*x = x_ba^2 = x_a^2+x_b^2*/
  x_ab = x[0] + x[1];

  a= (gammaCab*x_ab/(1+gammaCab*x_ab));

  *g_ab = c_ab0*pow(a,0)+ c_ab1*pow(a,1)+ c_ab2*pow(a,2)+c_ab3*pow(a,3)+c_ab4*pow(a,4);

  double dadx = gammaCab/pow(1+gammaCab*x_ab, 2.);
  dg_abdx = (0.0*c_ab0*pow(a,-1)+ 1.*c_ab1*pow(a,0)+ 2.*c_ab2*pow(a,1)+3.*c_ab3*pow(a,2)+4.*c_ab4*pow(a,3))*dadx;
    
  c_m06l_zx(x[0], 0.0, rho[0], 0.0, &dxdd_a, &dxdgd_a, &dzdd_a, &dzdtau_a);
  c_m06l_zx(x[1], 0.0, rho[1], 0.0, &dxdd_b, &dxdgd_b, &dzdd_b, &dzdtau_b);

  dg_abdd[0] = dg_abdx*dxdd_a; 
  dg_abdd[1] = dg_abdx*dxdd_b; 
  dg_abdgd[0] = dg_abdx*dxdgd_a; 
  dg_abdgd[1] = 0.0;
  dg_abdgd[2] = dg_abdx*dxdgd_b; 
}

/* Get g for Eq. (15)*/
static void
c_m06_15(double x, double rho, double *g_ss, double *dg_ssdd, double *dg_ssdgd)
{
  /*define the C_ss,i */
  static double c_ss0=0.5349466, c_ss1=0.5396620, c_ss2=-31.61217, c_ss3= 51.49592, c_ss4=-29.19613;
  double gammaCss = 0.06 ;
  double a; 
  double dg_ssdx, dxdd, dxdgd, dzdd, dzdtau;

  /*x = x_a^2 */

  a= (gammaCss*x/(1+gammaCss*x));

  *g_ss = c_ss0*pow(a,0)+ c_ss1*pow(a,1)+ c_ss2*pow(a,2)+c_ss3*pow(a,3)+c_ss4*pow(a,4);

  double dadx = gammaCss/pow(1+gammaCss*x, 2.);
  dg_ssdx = (0.0*c_ss0*pow(a,-1)+ 1.*c_ss1*pow(a,0)+ 2.*c_ss2*pow(a,1)+3.*c_ss3*pow(a,2)+4.*c_ss4*pow(a,3))*dadx;

  c_m06l_zx(x, 0.0, rho, 0.0, &dxdd, &dxdgd, &dzdd, &dzdtau);

  *dg_ssdd = dg_ssdx*dxdd; 
  *dg_ssdgd = dg_ssdx*dxdgd; 
  /*printf("g_ss %19.12f\n", *g_ss);*/
    
}

/* Get h_ab for Eq. (12)*/
static 
void c_m06l_hab(double *x, double *z, double *rho, double *tau, double *h_ab, double *dh_abdd, double *dh_abdgd, double *dh_abdtau)
{
  /* define the d_ab,i for Eq. (12)*/
  static double d_ab0= 0.3957626, d_ab1= -0.5614546, d_ab2= 0.01403963, d_ab3= 0.0009831442, d_ab4= -0.003577176;
  double alpha_ab = 0.00304966; 
  double hab1, dhabdd1[2], dhabdgd1[3], dhabdtau1[2];
  double x_ab, z_ab, gamma, xgamma, zgamma;
  double dgammadx, dgammadz;
  double dgammadd_a, dgammadgd_a, dgammadtau_a;
  double dgammadd_b, dgammadgd_b, dgammadtau_b;
  double dxdd_a, dxdgd_a, dzdd_a, dzdtau_a;
  double dxdd_b, dxdgd_b, dzdd_b, dzdtau_b;
  
  x_ab = x[0] + x[1];
  z_ab = z[0] + z[1];
  gamma = 1 + alpha_ab*(x_ab + z_ab);
  { /* derivatives of gamma with respect to x and z*/ 
    dgammadx = alpha_ab;        
    dgammadz = alpha_ab;
  }

  c_m06l_zx(x[0], z[0], rho[0], tau[0], &dxdd_a, &dxdgd_a, &dzdd_a, &dzdtau_a);
  c_m06l_zx(x[1], z[1], rho[1], tau[1], &dxdd_b, &dxdgd_b, &dzdd_b, &dzdtau_b);

  { /*derivatives of gamma with respect to density, gradient and kietic energy*/
    dgammadd_a   = dgammadx * dxdd_a + dgammadz * dzdd_a;
    dgammadd_b   = dgammadx * dxdd_b + dgammadz * dzdd_b;
    dgammadgd_a  = dgammadx * dxdgd_a;
    dgammadgd_b  = dgammadx * dxdgd_b;
    dgammadtau_a = dgammadz * dzdtau_a;
    dgammadtau_b = dgammadz * dzdtau_b;
  }

  xgamma = x_ab/gamma;
  zgamma = z_ab/gamma;

  /* we initialize h and collect the terms*/
  hab1    = 0.0;
  dhabdd1[0]   = dhabdd1[1]   = 0.0;
  dhabdgd1[0]  = dhabdgd1[1]  = dhabdgd1[2] = 0.0;
  dhabdtau1[0] = dhabdtau1[1] = 0.0;


  { /* first term */
    double g2=pow(gamma,2.);

    hab1         +=  d_ab0/gamma; 
    dhabdd1[0]   += -d_ab0*dgammadd_a/g2;
    dhabdd1[1]   += -d_ab0*dgammadd_b/g2;
    dhabdgd1[0]  += -d_ab0*dgammadgd_a/g2;
    dhabdgd1[1]  +=  0.0;
    dhabdgd1[2]  += -d_ab0*dgammadgd_b/g2;
    dhabdtau1[0] += -d_ab0*dgammadtau_a/g2 ;
    dhabdtau1[1] += -d_ab0*dgammadtau_b/g2 ;
  }

  { /* second term */
    double g3=pow(gamma,3.);

    hab1         += (d_ab1*xgamma + d_ab2*zgamma)/gamma;
    dhabdd1[0]   += (gamma*(d_ab1*dxdd_a+d_ab2*dzdd_a)-2*dgammadd_a*(d_ab1*x_ab+d_ab2*z_ab))/g3;
    dhabdd1[1]   += (gamma*(d_ab1*dxdd_b+d_ab2*dzdd_b)-2*dgammadd_b*(d_ab1*x_ab+d_ab2*z_ab))/g3;
    dhabdgd1[0]  += (d_ab1*dxdgd_a*gamma -2*(d_ab1*x_ab+d_ab2*z_ab)*dgammadgd_a)/g3;
    dhabdgd1[1]  += 0.0;
    dhabdgd1[2]  += (d_ab1*dxdgd_b*gamma -2*(d_ab1*x_ab+d_ab2*z_ab)*dgammadgd_b)/g3;
    dhabdtau1[0] += (d_ab2*dzdtau_a*gamma -2*(d_ab1*x_ab+d_ab2*z_ab)*dgammadtau_a)/g3;
    dhabdtau1[1] += (d_ab2*dzdtau_b*gamma -2*(d_ab1*x_ab+d_ab2*z_ab)*dgammadtau_b)/g3;
  }
  
  { /* third term */
    double g4= pow(gamma,4);

    hab1      += (d_ab3*xgamma*xgamma+d_ab4*xgamma*zgamma)/gamma;
    dhabdd1[0]   += (-3*dgammadd_a*(d_ab3*pow(x_ab,2.)+d_ab4*x_ab*z_ab)+dxdd_a*gamma*(2*d_ab3*x_ab+d_ab4*z_ab)+d_ab4*x_ab*dzdd_a*gamma)/g4;
    dhabdd1[1]   += (-3*dgammadd_b*(d_ab3*pow(x_ab,2.)+d_ab4*x_ab*z_ab)+dxdd_b*gamma*(2*d_ab3*x_ab+d_ab4*z_ab)+d_ab4*x_ab*dzdd_b*gamma)/g4;
    dhabdgd1[0]  += (-3*x_ab*(d_ab3*x_ab+d_ab4*z_ab)*dgammadgd_a+gamma*(2*d_ab3*x_ab+d_ab4*z_ab)*dxdgd_a)/g4;
    dhabdgd1[1]  += 0.0;
    dhabdgd1[2]  += (-3*x_ab*(d_ab3*x_ab+d_ab4*z_ab)*dgammadgd_b+gamma*(2*d_ab3*x_ab+d_ab4*z_ab)*dxdgd_b)/g4;
    dhabdtau1[0] += (d_ab4*x_ab*dzdtau_a*gamma-3*x_ab*(d_ab3*x_ab+d_ab4*z_ab)*dgammadtau_a)/g4;
    dhabdtau1[1] += (d_ab4*x_ab*dzdtau_b*gamma-3*x_ab*(d_ab3*x_ab+d_ab4*z_ab)*dgammadtau_b)/g4;
  }
  *h_ab = hab1;
  //derivatives
  dh_abdd[0]   = dhabdd1[0];
  dh_abdd[1]   = dhabdd1[1];
  dh_abdgd[0]  = dhabdgd1[0];
  dh_abdgd[1]  = dhabdgd1[1];
  dh_abdgd[2]  = dhabdgd1[2];
  dh_abdtau[0] = dhabdtau1[0];
  dh_abdtau[1] = dhabdtau1[1];

}

/* Get h_ss for Eq. (14)*/
static 
void c_m06l_hss(double x, double z, double rho, double tau, double *h_ss, double *dh_ssdd, double *dh_ssdgd, double *dh_ssdtau)
{
  /* define the d_ab,i for Eq. (12)*/
  static double d_ss0= 0.4650534, d_ss1= 0.1617589, d_ss2= 0.1833657, d_ss3= 0.0004692100, d_ss4= -0.004990573;
  double alpha_ss = 0.00515088; 
  double hss1, dhssdd1, dhssdgd1, dhssdtau1;
  double gamma, xgamma, zgamma;
  double dgammadx, dgammadz;
  double dgammadd, dgammadgd, dgammadtau;
  double dxdd, dxdgd, dzdd, dzdtau;
  

  gamma = 1 + alpha_ss*(x + z);
  { /* derivatives of gamma with respect to x and z*/ 
    dgammadx = alpha_ss;        
    dgammadz = alpha_ss;
  }

  c_m06l_zx(x, z, rho, tau, &dxdd, &dxdgd, &dzdd, &dzdtau);

  { /* derivatives of gamma with respect to density, gradient and kinetic energy */ 
    dgammadd   = dgammadx * dxdd + dgammadz * dzdd;
    dgammadgd  = dgammadx * dxdgd;
    dgammadtau = dgammadz * dzdtau;
  }

  xgamma = x/gamma;
  zgamma = z/gamma;

  /* we initialize h and collect the terms*/
  hss1    = 0.0;
  dhssdd1 = 0.0;
  dhssdgd1 = 0.0;
  dhssdtau1 = 0.0;


  { /* first term */
    double g2=pow(gamma,2.);

    hss1    +=  d_ss0/gamma; 
    dhssdd1   += -d_ss0*dgammadd/g2;
    dhssdgd1  += -d_ss0*dgammadgd/g2;
    dhssdtau1 += -d_ss0*dgammadtau/g2 ;
  }


  { /* second term */
    double g3=pow(gamma,3.);

    hss1      += (d_ss1*xgamma + d_ss2*zgamma)/gamma;
    dhssdd1   += (gamma*(d_ss1*dxdd+d_ss2*dzdd)-2*dgammadd*(d_ss1*x+d_ss2*z))/g3;
    dhssdgd1  += (d_ss1*dxdgd*gamma -2*(d_ss1*x+d_ss2*z)*dgammadgd)/g3;
    dhssdtau1 += (d_ss2*dzdtau*gamma -2*(d_ss1*x+d_ss2*z)*dgammadtau)/g3;
  }
  
  { /* third term */
    double g4= pow(gamma,4);

    hss1    += (d_ss3*xgamma*xgamma+d_ss4*xgamma*zgamma)/gamma;
    dhssdd1   += (-3*dgammadd*(d_ss3*pow(x,2.)+d_ss4*x*z)+dxdd*gamma*(2*d_ss3*x+d_ss4*z)+d_ss4*x*dzdd*gamma)/g4;
    dhssdgd1  += (-3*x*(d_ss3*x+d_ss4*z)*dgammadgd+gamma*(2*d_ss3*x+d_ss4*z)*dxdgd)/g4;
    dhssdtau1 += (d_ss4*x*dzdtau*gamma-3*x*(d_ss3*x+d_ss4*z)*dgammadtau)/g4;
  }
  *h_ss = hss1;
  //derivatives
  *dh_ssdd   = dhssdd1;
  *dh_ssdgd  = dhssdgd1;
  *dh_ssdtau = dhssdtau1;


}


static void 
c_m06l_para(m06l_params *p, const double *rho, const double *sigmatmp, const double *tautmp,
	    double *energy, double *dedd, double *vsigma, double *dedtau)
{
  double rho2[2], rho2s[2], x[2], z[2], zc_ss[2];
  double tau2[2], tauw[2], dens, dens1, sigma[3];
  double g_ss[2], h_ss[2], Ec_ss[2], D_ss[2]; 
  double g_ab=0.0, h_ab=0.0, Ec_ab=0.0; 
  double exunif_ss[2], vxunif_up[2], vxunif_dn[2], vxunif_ss[2];
  double exunif =0.0, exunif_ab=0.0, vxunif[2];
  //derivatives
  double dh_ssdd[2], dh_ssdgd[3], dh_ssdtau[2];
  double dg_ssdd[2], dg_ssdgd[3] ;
  double dh_abdd[2], dh_abdgd[3], dh_abdtau[2];
  double dg_abdd[2], dg_abdgd[3];
  double dEc_ssdd[2], dEc_ssdgd[3], dEc_ssdtau[2];
  double dEc_abdd[2], dEc_abdgd[3], dEc_abdtau[2];
  double dD_ssdd[2], dD_ssdgd[3], dD_ssdtau[2], dD_ssdx[2], dD_ssdz[2];
  double dxdd[2], dxdgd[2], dzdd[2], dzdtau[2];

  const double Cfermi= (3./5.)*pow(6*M_PI*M_PI,2./3.); 
  
  /* put in by cpo for const reasons */
  double sigma_[3],tau[2];
  sigma_[0] = sigmatmp[0];
  sigma_[1] = sigmatmp[1];
  sigma_[2] = sigmatmp[2];
  tau[0] = tautmp[0];
  tau[1] = tautmp[1];

  /*calculate |nabla rho|^2 */
  sigma_[0] = max(MIN_GRAD*MIN_GRAD, sigma_[0]);
  tauw[0] = max(sigma_[0]/(8.0*rho[0]), 1.0e-12);
  tau[0] = max(tauw[0], tau[0]);


  dens1 = rho[0]+rho[1];

  if(p->common.nspin== XC_UNPOLARIZED)
    {
      tau[1]  = 0.0; 

      rho2[0] = rho[0]/2.;
      rho2[1] = rho[0]/2.;	
      sigma[0] = sigma_[0]/4.;
      sigma[1] = sigma_[0]/4.;
      sigma[2] = sigma_[0]/4.;
      dens = rho[0];

      tau2[0] = tau[0]/2.;
      tau2[1] = tau[0]/2.;

    }else{
    sigma_[2] = max(MIN_GRAD*MIN_GRAD, sigma_[2]);
    tauw[1] = max(sigma_[2]/(8.0*rho[1]), 1.0e-12);
    tau[1] = max(tauw[1], tau[1]);

    rho2[0]=rho[0];
    rho2[1]=rho[1];	
    sigma[0] = sigma_[0];
    sigma[1] = sigma_[1];
    sigma[2] = sigma_[2];
    dens = rho[0]+rho[1];

    tau2[0] =tau[0];
    tau2[1] =tau[1];

  }
  //get the e_LDA(rho_a,b)
  const int np = 1;
  XC(lda_exc_vxc)(p->c_aux, np, rho2, &exunif, vxunif);
  exunif = exunif*dens;

  /*==============get the E_sigma part================*/
  /*============ spin up =============*/

  rho2s[0]=rho2[0];
  rho2s[1]=0.;	

  //get the e_LDA(rho_up,0)
  XC(lda_exc_vxc)(p->c_aux, np, rho2s, &(exunif_ss[0]), vxunif_up);
  exunif_ss[0] = exunif_ss[0] * rho2s[0];
  vxunif_ss[0] = vxunif_up[0];

  /*define variables for rho_up and zc in order to avoid x/0 -> D_ss = -inf */
  x[0] = sigma[0]/(pow(rho2s[0], 8./3.)); 
  z[0] = 2*tau2[0]/pow(rho2s[0],5./3.) - Cfermi;
  zc_ss[0] = 2*tau2[0]/pow(rho2s[0],5./3.);

  /*D_ss = 1 -x/4*(z + Cf), z+Cf = 2*tau2/pow(rho2s[0],5./3.) = zc */
  D_ss[0] = 1 - x[0]/(4. * zc_ss[0]);
  //derivatives for D_up
  dD_ssdx[0] = -1/(4 * zc_ss[0]);
  dD_ssdz[0] =  4 * x[0]/pow(4.*zc_ss[0],2.);

  c_m06l_zx(x[0], z[0], rho2s[0], tau2[0], &(dxdd[0]), &(dxdgd[0]), &(dzdd[0]), &(dzdtau[0]));
	  
  dD_ssdd[0]   = dD_ssdx[0] * dxdd[0] + dD_ssdz[0] * dzdd[0];
  dD_ssdgd[0]  = dD_ssdx[0] * dxdgd[0];
  dD_ssdtau[0] = dD_ssdz[0] * dzdtau[0];

  /*build up Eq. (14): Ec_sigmasigma*/
  c_m06_15(x[0], rho2s[0], &(g_ss[0]), &(dg_ssdd[0]), &(dg_ssdgd[0]));
  c_m06l_hss(x[0], z[0], rho2s[0], tau2[0], &(h_ss[0]), &(dh_ssdd[0]), &(dh_ssdgd[0]), &(dh_ssdtau[0]));

  Ec_ss[0] = (exunif_ss[0] * (g_ss[0]+h_ss[0]) * D_ss[0]);
  //printf("Ec_up %.9e\n", Ec_ss[0]);

  /*============== spin down =============*/

  rho2s[0]=rho2[1];
  rho2s[1]=0.;	
	  
  //get the e_LDA(0,rho_dn)
  XC(lda_exc_vxc)(p->c_aux, np, rho2s, &(exunif_ss[1]), vxunif_dn);
  exunif_ss[1] = exunif_ss[1] * rho2s[0];
  vxunif_ss[1] = vxunif_dn[0];

  /*define variables for rho_beta*/
  x[1] = sigma[2]/(pow(rho2s[0], 8./3.)); 
  z[1] = 2*tau2[1]/pow(rho2s[0],5./3.) - Cfermi;
  zc_ss[1] = 2*tau2[1]/pow(rho2s[0],5./3.);

  //printf("x1 %.9e, zc_ss%.9e\n", x[1], zc_ss[1]);
  D_ss[1] = 1 - x[1]/(4.*zc_ss[1]);
  //derivatives for D_dn
  dD_ssdx[1] = - 1/(4*zc_ss[1]);
  dD_ssdz[1] = 4*x[1]/pow(4.*zc_ss[1],2.);

  c_m06l_zx(x[1], z[1], rho2s[0], tau2[1], &(dxdd[1]), &(dxdgd[1]), &(dzdd[1]), &(dzdtau[1]));
	  
  dD_ssdd[1]   = dD_ssdx[1] * dxdd[1] + dD_ssdz[1] * dzdd[1];
  dD_ssdgd[2]  = dD_ssdx[1] * dxdgd[1];
  dD_ssdtau[1] = dD_ssdz[1] * dzdtau[1];

  c_m06_15(x[1], rho2s[0], &(g_ss[1]), &(dg_ssdd[1]), &(dg_ssdgd[2]));
  c_m06l_hss(x[1], z[1], rho2s[0], tau2[1], &(h_ss[1]), &(dh_ssdd[1]), &(dh_ssdgd[2]), &(dh_ssdtau[1]));


  //printf("exunif_ss %.9e, (g_ss[1]+h_ss[1])%.9e, D_ss %.9e\n", exunif_ss[1],(g_ss[1]+h_ss[1]),D_ss[1]);
  Ec_ss[1] = (exunif_ss[1] * (g_ss[1]+h_ss[1]) * D_ss[1]);
  //printf("Ec_dn %.9e\n", Ec_ss[1]);
	  
  // Derivatives for Ec_up and Ec_dn with respect to density and kinetic energy
  int i;
  for(i=0; i<2; i++){

    dEc_ssdd[i]   = exunif_ss[i] * dh_ssdd[i] * D_ss[i] + vxunif_ss[i] * h_ss[i] * D_ss[i] + exunif_ss[i] * h_ss[i] * dD_ssdd[i] +
      exunif_ss[i] * dg_ssdd[i] * D_ss[i] + vxunif_ss[i] * g_ss[i] * D_ss[i] + exunif_ss[i] * g_ss[i] * dD_ssdd[i];

    dEc_ssdtau[i] = exunif_ss[i] * dh_ssdtau[i] * D_ss[i] + exunif_ss[i] * h_ss[i] * dD_ssdtau[i] + exunif_ss[i] * g_ss[i] * dD_ssdtau[i];

  }
  // Derivatives for Ec_up and Ec_dn with respect to gradient
  dEc_ssdgd[0]  = exunif_ss[0] * dh_ssdgd[0] * D_ss[0] + exunif_ss[0] * h_ss[0] * dD_ssdgd[0] +
    exunif_ss[0] * dg_ssdgd[0] * D_ss[0] + exunif_ss[0] * g_ss[0] * dD_ssdgd[0];
  dEc_ssdgd[2]  = exunif_ss[1] * dh_ssdgd[2] * D_ss[1] + exunif_ss[1] * h_ss[1] * dD_ssdgd[2] + 
    exunif_ss[1] * dg_ssdgd[2] * D_ss[1] + exunif_ss[1] * g_ss[1] * dD_ssdgd[2];

	  
  /*==============get the E_ab part========================*/

  exunif_ab = exunif - exunif_ss[0] - exunif_ss[1];

  //x_ab = sigmatot[0] /(pow(rho2[0], 8./3.)) + sigmatot[2] /(pow(rho2[1], 8./3.));
  //z_ab = 2*tau2[0]/pow(rho2[0],5./3.) + 2*tau2[1]/pow(rho2[1],5./3.) - 2*Cfermi;

  /*build up Eq. (12): Ec_alphabeta*/
  c_m06_13(x, rho2, &g_ab, dg_abdd, dg_abdgd);
  c_m06l_hab(x, z, rho2, tau2, &h_ab, dh_abdd, dh_abdgd, dh_abdtau);

  Ec_ab = exunif_ab * (g_ab+h_ab);

  // Derivatives for Ec_ab with respect to density and kinetic energy
  for(i=0; i<2; i++){

    dEc_abdd[i]   = exunif_ab * (dh_abdd[i]+ dg_abdd[i]) + (vxunif[i]- vxunif_ss[i]) * (g_ab+h_ab);
    dEc_abdtau[i] = exunif_ab * dh_abdtau[i];

  }
  // Derivatives for Ec_ab with respect to gradient
  for(i=0; i<3; i++){
    dEc_abdgd[i] = exunif_ab * (dh_abdgd[i] + dg_abdgd[i]);
  }

  /*==============get the total energy E_c= E_up + E_dn + E_ab========================*/
  /*==============================and derivatives=====================================*/

  *energy = (Ec_ss[0] + Ec_ss[1] + Ec_ab)/dens1;
  //printf("Ec_ss %.9e, Ec_ss %.9e, Ec_ab %.9e\n", Ec_ss[0], Ec_ss[1], Ec_ab);

	  
  //derivative for the total correlation energy
  if(p->common.nspin== XC_UNPOLARIZED)
    {
      dedd[0]=dEc_ssdd[0] + dEc_abdd[0];
      dedd[1]=0.0;
	  
      vsigma[0]= (dEc_ssdgd[0] + dEc_abdgd[0])/2.;
      vsigma[1]= 0.0;
      vsigma[2]= 0.0;

      dedtau[0]= dEc_ssdtau[0] + dEc_abdtau[0];
      dedtau[1]= 0.0;
    }else{
    dedd[0]=dEc_ssdd[0] + dEc_abdd[0];
    dedd[1]=dEc_ssdd[1] + dEc_abdd[1];
	  
    vsigma[0]= dEc_ssdgd[0] + dEc_abdgd[0];
    vsigma[1]= 0.0;
    vsigma[2]= dEc_ssdgd[2] + dEc_abdgd[2];

    dedtau[0]= dEc_ssdtau[0] + dEc_abdtau[0];
    dedtau[1]= dEc_ssdtau[1] + dEc_abdtau[1];
  }
}

static void 
XC(mgga_c_m06l)(void *p, const double *rho, const double *sigma, const double *tau,
                double *e, double *dedd, double *vsigma, double *dedtau)
{
  c_m06l_para(p, rho, sigma, tau, e, dedd, vsigma, dedtau);
}

/* derivatives of x and z with respect to rho, grho and tau: Eq.(1) and Eq.(3)*/
static void 
x_m06l_zx(double x, double z, double rho, double tau, double *dxdd, double *dxdgd, double *dzdd, double *dzdtau)
{
  *dxdd = -8./3. * x * 1/rho;
  *dxdgd = 1./pow(rho,8./3.);

  *dzdd = -5./3. * 2* tau/pow(rho, 8./3.);
  *dzdtau = 2./pow(rho, 5./3.);
}

/* Build gamma and its derivatives with respect to rho, grho and tau: Eq. (4)*/
static void 
x_m06l_gamma(double x, double z, double rho, double tau, double *gamma, double *dgammadd, double *dgammadgd, double *dgammadtau)
{
  static double alpha = 0.00186726;   /*set alpha of Eq. (4)*/
  double dgammadx, dgammadz;
  double dxdd, dxdgd, dzdd, dzdtau;

  *gamma = 1 + alpha*(x + z);
  /*printf("gamma %19.12f\n", *gamma);*/

  { /* derivatives */ 
    dgammadx = alpha;        
    dgammadz = alpha;
  }

  x_m06l_zx(x, z, rho, tau, &dxdd, &dxdgd, &dzdd, &dzdtau);

  {
    *dgammadd = dgammadx*dxdd + dgammadz*dzdd;
    *dgammadgd = dgammadx*dxdgd;
    *dgammadtau = dgammadz*dzdtau;
  }

}

/************************************************************************
 Implements Zhao, Truhlar
   Meta-gga M06-Local

  Correlation part
************************************************************************/

/* calculate h and h derivatives with respect to rho, grho and tau: Equation (5) */
static 
void x_m06l_h(double x, double z, double rho, double tau, double *h, double *dhdd, double *dhdgd, double *dhdtau)
{
  /* parameters for h(x_sigma,z_sigma) of Eq. (5)*/
  static double d0=0.6012244, d1=0.004748822, d2=-0.008635108, d3=-0.000009308062, d4=0.00004482811;

  double h1, dhdd1, dhdgd1, dhdtau1;
  double gamma, dgammadd, dgammadgd, dgammadtau;
  double xgamma, zgamma;
  double dxdd, dxdgd, dzdd, dzdtau;
  
  x_m06l_gamma(x, z, rho, tau, &gamma, &dgammadd, &dgammadgd, &dgammadtau);

  xgamma = x/gamma;
  zgamma = z/gamma;

  /* we initialize h and its derivatives and collect the terms*/
  h1    = 0.0;
  dhdd1 = 0.0;
  dhdgd1 = 0.0;
  dhdtau1 = 0.0;


  { /* first term */
    double g2=pow(gamma,2.);

    h1      += d0/gamma; 
    dhdd1   += -d0*dgammadd/g2;
    dhdgd1  += -d0*dgammadgd/g2;
    dhdtau1 += -d0*dgammadtau/g2 ;
  }

  x_m06l_zx(x, z, rho, tau, &dxdd, &dxdgd, &dzdd, &dzdtau);
  
  { /* second term */
    double g3=pow(gamma,3.);

    h1      += (d1*xgamma + d2*zgamma)/gamma;
    dhdd1   += (gamma*(d1*dxdd+d2*dzdd)-2*dgammadd*(d1*x+d2*z))/g3;
    dhdgd1  += (d1*dxdgd*gamma -2*(d1*x+d2*z)*dgammadgd)/g3;
    dhdtau1 += (d2*dzdtau*gamma -2*(d1*x+d2*z)*dgammadtau)/g3;
  }
  
  { /* third term */
    double g4= pow(gamma,4);

    h1      += (d3*xgamma*xgamma+d4*xgamma*zgamma)/gamma;
    dhdd1   += (-3*dgammadd*(d3*pow(x,2.)+d4*x*z)+dxdd*gamma*(2*d3*x+d4*z)+d4*x*dzdd*gamma)/g4;
    dhdgd1  += (-3*x*(d3*x+d4*z)*dgammadgd+gamma*(2*d3*x+d4*z)*dxdgd)/g4;
    dhdtau1 += (d4*x*dzdtau*gamma-3*x*(d3*x+d4*z)*dgammadtau)/g4;
  }
  *h = h1;
  /*printf(" h %19.12f\n", *h);*/
  *dhdd = dhdd1;
  *dhdgd =dhdgd1;
  *dhdtau = dhdtau1;

}

/* f(w) and its derivatives with respect to rho and tau*/
static void 
x_m06l_fw(double rho, double tau, double *fw, double *dfwdd, double *dfwdtau)
{
  /*define the parameters for fw of Eq. (8) as in the reference paper*/
  static double a0= 0.3987756, a1= 0.2548219, a2= 0.3923994, a3= -2.103655, a4= -6.302147, a5= 10.97615,
    a6= 30.97273,  a7=-23.18489,  a8=-56.73480,  a9=21.60364,  a10= 34.21814, a11= -9.049762;

  double tau_lsda, t, w;
  double dtdd, dtdtau; 
  double dfwdw, dwdt, dtau_lsdadd;
  double aux =  (3./10.) * pow((6*M_PI*M_PI),2./3.); /*3->6 for nspin=2 */
	

  tau_lsda = aux * pow(rho,5./3.); 
  t = tau_lsda/tau;
  dtdtau = -t/tau; 
  w = (t - 1)/(t + 1);

  *fw = a0*pow(w,0.)+a1*pow(w,1.)+a2*pow(w,2.)+a3*pow(w,3.)+a4*pow(w,4.)+
    + a5*pow(w,5.)+a6*pow(w,6.)+a7*pow(w,7.)+a8*pow(w,8.)+a9*pow(w,9.)+a10*pow(w,10.)+a11*pow(w,11.);

  dfwdw = 0.0*a0*pow(w,-1)+1.0*a1*pow(w,0.)+2.0*a2*pow(w,1.)+3.0*a3*pow(w,2.)+4.0*a4*pow(w,3.)+
    + 5.0*a5*pow(w,4.)+6.0*a6*pow(w,5.)+7.0*a7*pow(w,6.)+8.0*a8*pow(w,7.)+9.0*a9*pow(w,8.)+
    + 10*a10*pow(w,9.)+11*a11*pow(w,10.);

  dwdt = 2/pow((t + 1),2.);

  dtau_lsdadd = aux * 5./3.* pow(rho,2./3.);
  dtdd = dtau_lsdadd/tau;

  *dfwdd =   dfwdw * dwdt * dtdd;
  *dfwdtau = dfwdw * dwdt * dtdtau;
}

static void 
x_m06l_para(m06l_params *pt, double rho, double sigma, double tau, double *energy, double *dedd, double *vsigma, double *dedtau)
{
  /*Build Eq. (6) collecting the terms Fx_PBE,  fw, e_lsda and h*/
  double grad, tauw, tau2, x, z;
  double rho2[2],sigmatot[3];
  double F_PBE, de_PBEdd[2], de_PBEdgd[3];
  double h, dhdd, dhdgd, dhdtau;
  double fw, dfwdd, dfwdtau;
  double epsx_lsda, depsx_lsdadd;
  const double Cfermi = (3./5.) * pow(6*M_PI*M_PI,2./3.);


  /* calculate |nabla rho|^2 */
  grad = sigma;
  grad = max(MIN_GRAD*MIN_GRAD, grad);
  tauw = max(grad/(8.0*rho),1.0e-12); /* tau^W = |nabla rho|^2/ 8rho */
  tau = max(tau, tauw);

  rho2[0]=rho/2.;
  rho2[1]=0.0;
  sigmatot[0] = grad/4.;
  sigmatot[1] = 0.0;
  sigmatot[2] = 0.0;
  tau2 =tau/2.;
  

  /* get the uniform gas energy and potential a MINUS was missing in the paper*/

  epsx_lsda = -(3./2.)*pow(3./(4*M_PI),1./3.)*pow(rho2[0],4./3.);

  depsx_lsdadd = -2*pow(3./(4*M_PI),1./3.)*pow(rho2[0],1./3.);

  /*get Fx for PBE*/
  const int np = 1;
  XC(gga_exc_vxc)(pt->x_aux, np, rho2, sigmatot, &F_PBE, de_PBEdd, de_PBEdgd);


  /* define x and z from Eq. (1) and Eq. (3) NOTE: we build directly x^2 */
  x = grad/(4*pow(rho2[0], 8./3.));
  
  z = 2*tau2/pow(rho2[0],5./3.) - Cfermi;  /*THERE IS A 2 IN FRONT AS IN THEOR. CHEM. ACCOUNT 120 215 (2008)*/
  
  /*get  h and fw*/
  x_m06l_h(x, z, rho2[0], tau2, &h, &dhdd, &dhdgd, &dhdtau);
  x_m06l_fw(rho2[0], tau2, &fw, &dfwdd, &dfwdtau);


  
  { /* Eq. (6)  E_x = Int F_PBE*fw + exunif*h, the factor 2 accounts for spin. */

    *energy = 2*(F_PBE*rho2[0] *fw + epsx_lsda *h);

    *dedd   = (de_PBEdd[0] *fw + F_PBE*rho2[0] * dfwdd+ depsx_lsdadd *h + epsx_lsda * dhdd);
    *dedtau = (F_PBE * dfwdtau *rho2[0] + epsx_lsda * dhdtau);

    *vsigma = (de_PBEdgd[0] *fw +  epsx_lsda*dhdgd)/2.;
  }
}


void 
XC(mgga_x_m06l)(void *p, const double *rho, const double *sigma, const double *tau,
                double *e, double *dedd, double *vsigma, double *dedtau)

{
  m06l_params *par = (m06l_params*)p;
  if(par->common.nspin == XC_UNPOLARIZED){
    double en;
    x_m06l_para(p, rho[0], sigma[0], tau[0], &en, dedd, vsigma, dedtau);
    *e = en/(rho[0]+rho[1]);

  }else{
  
  
    *e = 0.0;

    double e2na, e2nb, rhoa[2], rhob[2];

    double vsigmapart[3]; 
	  

    rhoa[0]=2*rho[0];
    rhoa[1]=0.0;
    rhob[0]=2*rho[1];
    rhob[1]=0.0;


		  
    x_m06l_para(p, rhoa[0], 4*sigma[0], 2.0*tau[0], &e2na, &(dedd[0]), &(vsigmapart[0]), &(dedtau[0]));

    x_m06l_para(p, rhob[0], 4*sigma[2], 2.0*tau[1], &e2nb, &(dedd[1]), &(vsigmapart[2]), &(dedtau[1]));
		 
    *e = (e2na + e2nb )/(2.*(rho[0]+rho[1]));
    vsigma[0] = 2*vsigmapart[0];
    vsigma[2] = 2*vsigmapart[2];
  }
}

static void m06l_init(void *p)
{
  m06l_params *par = (m06l_params*)p;
  par->c_aux = (XC(func_type) *) malloc(sizeof(XC(func_type)));
  XC(func_init)(par->c_aux, XC_LDA_C_PW, XC_POLARIZED);

  par->x_aux = (XC(func_type) *) malloc(sizeof(XC(func_type)));
  XC(func_init)(par->x_aux, XC_GGA_X_PBE, XC_POLARIZED);
}

static void m06l_end(void *p)
{
  m06l_params *par = (m06l_params*)p;
  XC(func_end)(par->c_aux);
  free(par->c_aux);

  XC(func_end)(par->x_aux);
  free(par->x_aux);
}

const mgga_func_info m06l_info = {
  sizeof(m06l_params),
  &m06l_init,
  &m06l_end,
  &XC(mgga_x_m06l),
  &XC(mgga_c_m06l),
};
