
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "xc_mgga.h"
#include "xc_gpaw.h"

extern const mgga_func_info m06l_info;
extern const mgga_func_info tpss_info;
extern const mgga_func_info revtpss_info;

static void init_common(common_params* params, int code, int nspin, const mgga_func_info *finfo) {
  params->code = code;
  params->nspin = nspin;
  params->funcinfo = finfo;
}

void init_mgga(void** params, int code, int nspin) {
  const mgga_func_info *finfo;
  if (code==20) {
    finfo = &tpss_info;
  } else if (code==21) {
    finfo = &m06l_info;
  } else if (code==22) {
    finfo = &revtpss_info;
  } else {
    // this should never happen.  forces a crash.
    assert(code>=20 && code <=22);
    finfo = NULL;
  }
  *params = malloc(finfo->size);
  init_common(*params, code, nspin, finfo);
  finfo->init(*params);
}

void end_mgga(common_params *common) {
  common->funcinfo->end(common);
  free(common);
}

void calc_mgga(void** params, int nspin, int ng,
               const double* n_g, const double* sigma_g, const double* tau_g,
               double *e_g, double *v_g, double *dedsigma_g, double *dedtau_g) {

  common_params *common = (common_params*)*params;
  // check for a changed spin (similar to a line in gpaw/libxc.py)
  if (nspin!=common->nspin) {
    int code = common->code; // save this, since we're about to destroy common
    end_mgga(common);
    init_mgga(params, code, nspin);
    common = (common_params*)*params; // init_mgga changes this
  }

  if (nspin == 1) {
    for (int g = 0; g < ng; g++) {
      // kludge n[1] because of the way TPSS was written (requires n[1]=0.0 even for unpolarized)
      double n[2];
      n[0] = n_g[g];
      n[1] = 0.0;
      if (n[0] < NMIN) n[0] = NMIN;
      // m06l is assuming that there is space for spinpolarized calculation output
      // even for non-spin-polarized.
      double etmp, vtmp[2], dedsigmatmp[3], dedtautmp[2];
      common->funcinfo->exch(*params, n, sigma_g+g, tau_g+g,
                             &etmp, vtmp, dedsigmatmp, dedtautmp);
      e_g[g] = etmp;
      v_g[g] += vtmp[0];
      dedsigma_g[g] = dedsigmatmp[0];
      dedtau_g[g] = dedtautmp[0];
      common->funcinfo->corr(*params, n, sigma_g+g, tau_g+g,
                             &etmp, vtmp, dedsigmatmp, dedtautmp);
      e_g[g] += etmp;
      e_g[g] *= n[0];
      v_g[g] += vtmp[0];
      dedsigma_g[g] += dedsigmatmp[0];
      dedtau_g[g] += dedtautmp[0];
    }
  } else {
    double etmp, ntmp[2], vtmp[2], sigmatmp[3], dedsigmatmp[3],
      tautmp[2], dedtautmp[2];
    for (int g = 0; g < ng; g++) {
      ntmp[0] = n_g[g];
      if (ntmp[0] < NMIN) ntmp[0] = NMIN;
      ntmp[1] = n_g[g+ng];
      if (ntmp[1] < NMIN) ntmp[1] = NMIN;

      sigmatmp[0] = sigma_g[g];
      sigmatmp[1] = sigma_g[g+ng];
      sigmatmp[2] = sigma_g[g+ng+ng];

      tautmp[0] = tau_g[g];
      tautmp[1] = tau_g[g+ng];

      // kludge: mgga_x_tpss requires dedsigma[1] set to 0, since it doesn't calculate it.
      dedsigmatmp[1]=0.0;

      common->funcinfo->exch(*params, ntmp, sigmatmp, tautmp,
                             &etmp, vtmp, dedsigmatmp, dedtautmp);
      e_g[g] = etmp;
      v_g[g] += vtmp[0];
      v_g[g+ng] += vtmp[1];
      dedsigma_g[g] = dedsigmatmp[0];
      dedsigma_g[g+ng] = dedsigmatmp[1];
      dedsigma_g[g+ng+ng] = dedsigmatmp[2];
      dedtau_g[g] = dedtautmp[0];
      dedtau_g[g+ng] = dedtautmp[1];
      common->funcinfo->corr(*params, ntmp, sigmatmp, tautmp,
                             &etmp, vtmp, dedsigmatmp, dedtautmp);
      e_g[g] += etmp;
      e_g[g] *= ntmp[0]+ntmp[1];
      v_g[g] += vtmp[0];
      v_g[g+ng] += vtmp[1];
      dedsigma_g[g] += dedsigmatmp[0];
      dedsigma_g[g+ng] += dedsigmatmp[1];
      dedsigma_g[g+ng+ng] += dedsigmatmp[2];
      dedtau_g[g] += dedtautmp[0];
      dedtau_g[g+ng] += dedtautmp[1];
    }
  }
}
