
#ifndef GPAW_XC_MGGA_H
#define GPAW_XC_MGGA_H

#define M_PI           3.14159265358979323846
#define MIN_DENS             1.0e-20
#define MIN_GRAD             1.0e-20
#define max(x,y)  ((x<y) ? y : x)

struct xc_func_type;

typedef struct gga_func_info {
  void (*init)(struct xc_func_type *par);
  void (*end)(struct xc_func_type *par);
  void (*corr)(void* par, const double *rho, const double *sigma,
               double *energy, double *dedd, double *vsigma,
               double *v2rho2, double *v2rhosigma, double *v2sigma2);
} gga_func_info;

typedef struct mgga_func_info {
  int size;
  void (*init)(void* par);
  void (*end)(void* par);
  void (*exch)(void* par, const double *rho, const double *sigma, const double *tau,
               double *energy, double *dedd, double *vsigma, double *dedtau);
  void (*corr)(void* par, const double *rho, const double *sigma, const double *tau,
               double *energy, double *dedd, double *vsigma, double *dedtau);
} mgga_func_info;
  
typedef struct common_params {
  int nspin;
  int code;
  const mgga_func_info *funcinfo;
} common_params;

#endif
