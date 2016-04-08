/*  Copyright (C) 2003-2007  CAMP
 *  Copyright (C) 2007-2008  CAMd
 *  Please see the accompanying LICENSE file for further information. */

#include <Python.h>
#define PY_ARRAY_UNIQUE_SYMBOL GPAW_ARRAY_API
#define NO_IMPORT_ARRAY
#include <numpy/arrayobject.h>
#include <assert.h>
#include <xc.h>
#include "xc_gpaw.h"
#include "../extensions.h"

typedef struct
{
  PyObject_HEAD
  /* exchange-correlation energy second derivatives */
  void (*get_fxc)(XC(func_type) *func, double point[7], double der[5][5]);
  XC(func_type) xc_functional;
  XC(func_type) x_functional;
  XC(func_type) c_functional;
  XC(func_type) *functional[2]; /* store either x&c, or just xc */
  int nspin; /* must be common to x and c, so declared redundantly here */
} lxcXCFunctionalObject;

void XC(lda_fxc_fd)(const XC(func_type) *p, const double *rho, double *fxc);

/* a general call for an LDA functional - finite difference */
void get_fxc_fd_lda(XC(func_type) *func, double point[7], double der[5][5])
{
  double v2rho2[3], v2rhosigma[6], v2sigma2[6];

  for(int i=0; i<3; i++) v2rho2[i] = 0.0;
  for(int i=0; i<6; i++){
    v2rhosigma[i] = 0.0;
    v2sigma2[i]    = 0.0;
  }
  XC(lda_fxc_fd)(func, point, v2rho2);
  der[0][0] = v2rho2[0];
  der[0][1] = der[1][0] = v2rho2[1];
  der[1][1] = v2rho2[2];
  der[0][2] = der[2][0] = v2rhosigma[0];
  der[0][3] = der[3][0] = v2rhosigma[1];
  der[0][4] = der[4][0] = v2rhosigma[2];
  der[1][2] = der[2][1] = v2rhosigma[3];
  der[1][3] = der[3][1] = v2rhosigma[4];
  der[1][4] = der[4][1] = v2rhosigma[5];
  der[2][2] = v2sigma2[0];
  der[2][3] = der[3][2] = v2sigma2[1];
  der[2][4] = der[4][2] = v2sigma2[2];
  der[3][3] = v2sigma2[3];
  der[3][4] = der[4][3] = v2sigma2[4];
  der[4][4] = v2sigma2[5];
}

// finite difference calculation of second functional derivative
// stolen from libxc/testsuite/xc-consistency.c

double get_point(XC(func_type) *func, double point[7], double *e, double der[5], int which)
{
  const int np = 1;
  switch(func->info->family)
    {
    case XC_FAMILY_LDA:
      XC(lda_exc_vxc)(func, np, &(point[0]), e, &(der[0]));
      break;
    case XC_FAMILY_GGA:
    case XC_FAMILY_HYB_GGA:
      XC(gga_exc_vxc)(func, np, &(point[0]), &(point[2]),
                      e, &(der[0]), &(der[2]));
      break;
    }

  if(which == 0)
    return (*e)*(point[0] + point[1]);
  else
    return der[which-1];
}

void first_derivative(XC(func_type) *func, double point[7], double der[5], int which,
                      int nspin)
{
  int i;

  for(i=0; i<5; i++){
    const double delta = 5e-10;

    double dd, p[5], v[5];
    int j;

    if(nspin==1 && (i!=0 && i!=2)){
      der[i] = 0.0;
      continue;
    }

    dd = point[i]*delta;
    if(dd < delta) dd = delta;

    for(j=0; j<5; j++) p[j] = point[j];

    if(point[i]>=3.0*dd){ /* centered difference */
      double e, em1, em2, ep1, ep2;

      p[i] = point[i] + dd;
      ep1 = get_point(func, p, &e, v, which);

      p[i] = point[i] + 2*dd;
      ep2 = get_point(func, p, &e, v, which);

      p[i] = point[i] - dd;  /* backward point */
      em1 = get_point(func, p, &e, v, which);

      p[i] = point[i] - 2*dd;  /* backward point */
      em2 = get_point(func, p, &e, v, which);

      der[i]  = 1.0/2.0*(ep1 - em1);
      der[i] += 1.0/12.0*(em2 - 2*em1 + 2*ep1 - ep2);

      der[i] /= dd;

    }else{  /* we use a 5 point forward difference */
      double e, e1, e2, e3, e4, e5;

      p[i] = point[i];
      e1 = get_point(func, p, &e, v, which);

      p[i] = point[i] + dd;
      e2 = get_point(func, p, &e, v, which);

      p[i] = point[i] + 2.0*dd;
      e3 = get_point(func, p, &e, v, which);

      p[i] = point[i] + 3.0*dd;
      e4 = get_point(func, p, &e, v, which);

      p[i] = point[i] + 4.0*dd;
      e5 = get_point(func, p, &e, v, which);

      der[i]  =          (-e1 + e2);
      der[i] -=  1.0/2.0*( e1 - 2*e2 + e3);
      der[i] +=  1.0/3.0*(-e1 + 3*e2 - 3*e3 + e4);
      der[i] -=  1.0/4.0*( e1 - 4*e2 + 6*e3 - 4*e4 + e5);

      der[i] /= dd;
    }
  }
}

void first_derivative_spinpaired(XC(func_type) *func, double point[7], double der[5],
                                 int which)
{
  first_derivative(func, point, der, which, XC_UNPOLARIZED);
}

void first_derivative_spinpolarized(XC(func_type) *func, double point[7], double der[5],
                                    int which)
{
  first_derivative(func, point, der, which, XC_POLARIZED);
}

void second_derivatives_spinpaired(XC(func_type) *func, double point[7], double der[5][5])
{
  int i;

  for(i=0; i<5; i++){
    first_derivative_spinpaired(func, point, der[i], i+1);
  }
}

void second_derivatives_spinpolarized(XC(func_type) *func, double point[7], double der[5][5])
{
  int i;

  for(i=0; i<5; i++){
    first_derivative_spinpolarized(func, point, der[i], i+1);
  }
}

/* a general call for a functional - finite difference */
void get_fxc_fd_spinpaired(XC(func_type) *func, double point[7], double der[5][5])
{
  second_derivatives_spinpaired(func, point, der);
}

/* a general call for a functional - finite difference */
void get_fxc_fd_spinpolarized(XC(func_type) *func, double point[7], double der[5][5])
{
  second_derivatives_spinpolarized(func, point, der);
}

static void lxcXCFunctional_dealloc(lxcXCFunctionalObject *self)
{
  for (int i=0; i<2; i++)
    if (self->functional[i] != NULL) xc_func_end(self->functional[i]);

  PyObject_DEL(self);
}

static PyObject*
lxcXCFunctional_is_gga(lxcXCFunctionalObject *self, PyObject *args)
{
  int success = 0; /* assume functional is not GGA */
  // check family of most-complex functional
  if (self->functional[0]->info->family == XC_FAMILY_GGA ||
      self->functional[0]->info->family == XC_FAMILY_HYB_GGA) success = XC_FAMILY_GGA;
  return Py_BuildValue("i", success);
}

static PyObject*
lxcXCFunctional_is_mgga(lxcXCFunctionalObject *self, PyObject *args)
{
  int success = 0; /* assume functional is not MGGA */
  // check family of most-complex functional
  if (self->functional[0]->info->family == XC_FAMILY_MGGA) success = XC_FAMILY_MGGA;
  return Py_BuildValue("i", success);
}

static PyObject*
lxcXCFunctional_CalculateFXC_FD_SpinPaired(lxcXCFunctionalObject *self, PyObject *args)
{
  PyArrayObject* n_array;              /* rho */
  PyArrayObject* v2rho2_array;         /* d2E/drho2 */
  PyArrayObject* a2_array = 0;         /* |nabla rho|^2*/
  PyArrayObject* v2rhosigma_array = 0; /* d2E/drhod|nabla rho|^2 */
  PyArrayObject* v2sigma2_array = 0;   /* d2E/drhod|nabla rho|^2 */
  if (!PyArg_ParseTuple(args, "OO|OOO", &n_array, &v2rho2_array, /* object | optional objects*/
                        &a2_array, &v2rhosigma_array, &v2sigma2_array))
    return NULL;

  /* find nspin */
  int nspin = self->nspin;

  assert(nspin == XC_UNPOLARIZED); /* we are spinpaired */
  assert (self->functional[0]->info->family != XC_FAMILY_MGGA);

  int ng = PyArray_DIMS(n_array)[0]; /* number of grid points */

  const double* n_g = DOUBLEP(n_array); /* density on the grid */
  double* v2rho2_g = DOUBLEP(v2rho2_array); /* v on the grid */

  const double* a2_g = 0; /* a2 on the grid */
  double* v2rhosigma_g = 0;    /* d2Ednda2 on the grid */
  double* v2sigma2_g = 0;    /* d2Eda2da2 on the grid */

  if ((self->functional[0]->info->family == XC_FAMILY_GGA) ||
      (self->functional[0]->info->family == XC_FAMILY_HYB_GGA))
    {
      a2_g = DOUBLEP(a2_array);
      v2rhosigma_g = DOUBLEP(v2rhosigma_array);
      v2sigma2_g = DOUBLEP(v2sigma2_array);
    }

  self->get_fxc = get_fxc_fd_spinpaired;
  /* ################################################################ */
  for (int g = 0; g < ng; g++)
    {
      double n = n_g[g];
      if (n < NMIN)
        n = NMIN;
      double a2 = 0.0; /* initialize for lda */
      if ((self->functional[0]->info->family == XC_FAMILY_GGA) ||
          (self->functional[0]->info->family == XC_FAMILY_HYB_GGA))
        {
          a2 = a2_g[g];
        }

      double point[7]; /* generalized point */
      // from http://www.tddft.org/programs/octopus/wiki/index.php/Libxc:manual
      // rhoa rhob sigmaaa sigmaab sigmabb taua taub
      // \sigma[0] = \nabla n_\uparrow \cdot \nabla n_\uparrow \qquad
      // \sigma[1] = \nabla n_\uparrow \cdot \nabla n_\downarrow \qquad
      // \sigma[2] = \nabla n_\downarrow \cdot \nabla n_\downarrow \qquad

      double derivative[5][5]; /* generalized derivative */
      double v2rho2[3];
      double v2rhosigma[6];
      double v2sigma2[6];
      // one that uses this: please add description of spin derivative order notation
      // (see c/libxc/src/gga_perdew.c) MDTMP

      for(int i=0; i<3; i++) v2rho2[i] = 0.0;
      for(int i=0; i<6; i++){
        v2rhosigma[i] = 0.0;
        v2sigma2[i]    = 0.0;
      }

      for(int j=0; j<7; j++)
        {
          point[j] = 0.0;
        }

      for(int i=0; i<5; i++)
        {
          for(int j=0; j<5; j++)
            {
              derivative[i][j] = 0.0;
            }
        }

      point[0] = n;   /* -> rho */
      point[2] = a2;  /* -> sigma */

      for (int i=0; i<2; i++) {
        XC(func_type) *func = self->functional[i];
        if (func == NULL) continue;
        self->get_fxc(func, point, derivative);
        v2rho2[0] = derivative[0][0];
        v2rho2[1] = derivative[0][1]; // XC_POLARIZED
        v2rho2[2] = derivative[1][1]; // XC_POLARIZED
        v2rhosigma[0] = derivative[0][2];
        v2rhosigma[1] = derivative[0][3]; // XC_POLARIZED
        v2rhosigma[2] = derivative[0][4]; // XC_POLARIZED
        v2rhosigma[3] = derivative[1][2]; // XC_POLARIZED
        v2rhosigma[4] = derivative[1][3]; // XC_POLARIZED
        v2rhosigma[5] = derivative[1][4]; // XC_POLARIZED
        v2sigma2[0] = derivative[2][2]; /* aa_aa */
        v2sigma2[1] = derivative[2][3]; // XC_POLARIZED /* aa_ab */
        v2sigma2[2] = derivative[2][4]; // XC_POLARIZED /* aa_bb */
        v2sigma2[3] = derivative[3][3]; // XC_POLARIZED /* ab_ab */
        v2sigma2[4] = derivative[3][4]; // XC_POLARIZED /* ab_bb */
        v2sigma2[5] = derivative[4][4]; // XC_POLARIZED /* bb_bb */
        switch(func->info->family)
          {
          case XC_FAMILY_HYB_GGA:
          case XC_FAMILY_GGA:
            v2rhosigma_g[g] += v2rhosigma[0];
            v2sigma2_g[g] += v2sigma2[0];
            // don't break here since we need LDA values as well
          case XC_FAMILY_LDA:
            v2rho2_g[g] += v2rho2[0];
          }
      }

    }
  Py_RETURN_NONE;
}

// Below are changes made by cpo@slac.stanford.edu for libxc 1.2.0
// which allows passing of arrays of points to libxc routines.

// The fundamental design idea (to try to minimize code-duplication) is that
// all libxc routines have input/output arrays that get processed in
// common ways with three special exceptions: n_sg, e_g, dedn_sg.  The
// struct "xcptrlist" is used to keep track of these pointers.

// Two libxc features prevent us from using a straightforward
// interface:
// 1) libxc calls memset(0) on output arrays, preventing us
//    from adding x/c contributions "in place" without scratch arrays
// 2) for spin-polarized calculations libxc wants spin indices to be
//    dense in memory, whereas GPAW probably loops over grid indices
//    more often, so we want to keep those dense in memory.
// I asked Miguel Marques to remove the memset, and to add a "stride"
// argument to libxc routines to address the above.  He says he will
// consider it in the future.  In the meantime we have to "block"
// over gridpoints using some scratch memory.

// What is supported:
// - combined xc-functional mode
// - separate x,c functionals.
// - separate x,c can have differing complexities (e.g. one GGA, one LDA)
// - "exc_vxc" style routines for LDA/GGA/MGGA both unpolarized/polarized
// - "fxc" style routines for LDA/GGA both unpolarized/polarized
// To support a libxc routine other than exc_vxc/fxc one needs to
// copy a "Calculate" routine and change the pointer list setup, and
// associated libxc function calls.

// number of gridpoints we will "block" over when doing xc calculation
#define BLOCKSIZE 1024

// this is the maximum number of BLOCKSIZE arrays that will be put
// into scratch (depends on the "spinsize" values for the various
// arrays.  currently determined by fxc, which has input spinsizes
// of 2+3 and output spinsizes of 3+6+6 (totalling 20).
#define MAXARRAYS 20
#define LIBXCSCRATCHSIZE (BLOCKSIZE*MAXARRAYS)

static double *scratch=NULL;

// we don't use lapl, but libxc needs space for them.
static double *scratch_lapl=NULL;
static double *scratch_vlapl=NULL;

// special cases for array behaviors:
// flag to indicate we need to add to existing values for dedn_sg
#define DEDN_SG 1
// flag to indicate we need to apply NMIN cutoff to n_sg
#define N_SG 2
// flag to indicate we need to multiply by density for e_g
#define E_G 4

typedef struct xcptr {
  double *p;
  int special;
  int spinsize;
} xcptr;

#define MAXPTR 10

typedef struct xcptrlist {
  int num;
  xcptr p[MAXPTR];
} xcptrlist;

typedef struct xcinfo {
  int nspin;
  bool spinpolarized;
  int ng;
} xcinfo;

// these 3 functions make the spin index closest in memory ("gather") or the
// farthest apart in memory ("scatter").  "scatteradd" adds to previous results.

static void gather(const double* src, double* dst, int np, int stride, int nspins) {
  const double *dstend = dst+np*nspins;
  const double *srcend = src+nspins*stride;
  do {
    const double *s = src;
    do {
      *dst++ = *s; s+=stride;
    } while (s<srcend);
    src++; srcend++;
  } while (dst<dstend);
}

static void scatter(const double* src, double* dst, int np, int stride, int nspins) {
  const double *srcend = src+np*nspins;
  const double *dstend = dst+nspins*stride;
  do {
    double *d = dst;
    do {
      *d = *src++; d+=stride;
    } while (d<dstend);
    dst++; dstend++;
  } while (src<srcend);
}

static void scatteradd(const double* src, double* dst, int np, int stride, int nspins) {
  const double *srcend = src+np*nspins;
  const double *dstend = dst+nspins*stride;
  do {
    double *d = dst;
    do {
      *d += *src++; d+=stride;
    } while (d<dstend);
    dst++; dstend++;
  } while (src<srcend);
}

// set up the pointers into the scratch area, leaving space for each of the arrays

static void setupblockptrs(const xcinfo *info,
                           const xcptrlist *inlist, const xcptrlist *outlist,
                           double **inblocklist, double **outblocklist,
                           int blocksize) {
  // set up the block pointers we are going to use in the "scratch" space
  double *next = scratch;
  for (int i=0; i<inlist->num; i++) {
    inblocklist[i] = next;
    next+=blocksize*inlist->p[i].spinsize;
  }
  for (int i=0; i<outlist->num; i++) {
    outblocklist[i] = next;
    next+=blocksize*outlist->p[i].spinsize;
  }
  // check that we fit in the scratch space
  // if we don't, then we need to increase MAXARRAY
  assert((next - scratch) <= LIBXCSCRATCHSIZE);
}

// copy a piece of the full data into the block for processing by libxc

static void data2block(const xcinfo *info,
                       const xcptrlist *inlist, double *inblocklist[],
                       int blocksize) {

  // copy data into the block, taking into account special cases
  for (int i=0; i<inlist->num; i++) {
    double *ptr = inlist->p[i].p; double* block = inblocklist[i];
    if (info->spinpolarized) {
      gather(ptr,block,blocksize,info->ng,inlist->p[i].spinsize);
      if (inlist->p[i].special&N_SG)
        for (int i=0; i<blocksize*2; i++) block[i] = (block[i]<NMIN) ? NMIN : block[i];
    } else {
      // don't copy sigma and tau for non-spin-polarized.
      // use input arrays instead to save time. have to
      // copy n_g however, because of the NMIN patch.
      if (inlist->p[i].special&N_SG) for (int i=0; i<blocksize; i++) block[i] = (ptr[i]<NMIN) ? NMIN : ptr[i];
    }
  }
}

// copy the data from the block back into its final resting place

static void block2data(const xcinfo *info, double *outblocklist[], const xcptrlist *outlist,
                       const double *n_sg, int blocksize) {
  for (int i=0; i<outlist->num; i++) {
    double *ptr = outlist->p[i].p; double* block = outblocklist[i];
    if (outlist->p[i].special&E_G) {
      if (info->spinpolarized) {
        for (int i=0; i<blocksize; i++)
          ptr[i]=(n_sg[i*2]+n_sg[i*2+1])*block[i];
      } else {
        for (int i=0; i<blocksize; i++) ptr[i]=n_sg[i]*block[i];
      }
    } else if (outlist->p[i].special&DEDN_SG) {
      if (info->spinpolarized) {
        scatteradd(block,ptr,blocksize,info->ng,outlist->p[i].spinsize); // need to add to pre-existing values
      } else {
        for (int i=0; i<blocksize; i++) ptr[i]+=block[i]; // need to add to pre-existing values
      }
    } else {
      if (info->spinpolarized) {
        scatter(block,ptr,blocksize,info->ng,outlist->p[i].spinsize);
      } else {
        memcpy(ptr,block,blocksize*sizeof(double));
      }
    }
  }
}

// copy the data from the block back into its final resting place, but add to previous results

static void block2dataadd(const xcinfo *info, double *outblocklist[], const xcptrlist *outlist,
                          const double *n_sg, int blocksize, int noutcopy) {
  for (int i=0; i<noutcopy; i++) {
    double *ptr = outlist->p[i].p; double* block = outblocklist[i];
    if (outlist->p[i].special&E_G) {
      if (info->spinpolarized) {
        for (int i=0; i<blocksize; i++)
          ptr[i]+=(n_sg[i*2]+n_sg[i*2+1])*block[i];
      } else {
        for (int i=0; i<blocksize; i++) ptr[i]+=n_sg[i]*block[i];
      }
    } else {
      if (info->spinpolarized) {
        scatteradd(block,ptr,blocksize,info->ng,outlist->p[i].spinsize);
      } else {
        for (int i=0; i<blocksize; i++) ptr[i]+=block[i];
      }
    }
  }
}

static PyObject*
lxcXCFunctional_Calculate(lxcXCFunctionalObject *self, PyObject *args)
{
  PyArrayObject* py_n_sg = NULL;
  PyArrayObject* py_sigma_xg = NULL;
  PyArrayObject* py_e_g = NULL;
  PyArrayObject* py_dedn_sg = NULL;
  PyArrayObject* py_dedsigma_xg = NULL;
  PyArrayObject* py_tau_sg = NULL;
  PyArrayObject* py_dedtau_sg = NULL;

  if (!PyArg_ParseTuple(args, "OOO|OOOO", &py_e_g,
                        &py_n_sg, &py_dedn_sg,
                        &py_sigma_xg, &py_dedsigma_xg,
                        &py_tau_sg, &py_dedtau_sg))
    return NULL;

  xcinfo info;
  info.nspin = self->nspin;
  info.spinpolarized = (info.nspin==2);
  info.ng = PyArray_DIMS(py_e_g)[0];

  xcptrlist inlist,outlist;
  inlist.num=0;
  outlist.num=0;

  int blocksize = BLOCKSIZE;
  int remaining = info.ng;

  // setup pointers using most complex functional
  switch(self->functional[0]->info->family)
    {
    case XC_FAMILY_MGGA:
      inlist.p[2].p = DOUBLEP(py_tau_sg);
      inlist.p[2].special = 0;
      inlist.p[2].spinsize = 2;
      inlist.num++;
      outlist.p[3].p = DOUBLEP(py_dedtau_sg);
      outlist.p[3].special = 0;
      outlist.p[3].spinsize = 2;
      outlist.num++;
      // don't break here since MGGA also needs GGA ptrs
    case XC_FAMILY_HYB_GGA:
    case XC_FAMILY_GGA:
      inlist.p[1].p = DOUBLEP(py_sigma_xg);
      inlist.p[1].special = 0;
      inlist.p[1].spinsize = 3;
      inlist.num++;
      outlist.p[2].p = DOUBLEP(py_dedsigma_xg);
      outlist.p[2].special = 0;
      outlist.p[2].spinsize = 3;
      outlist.num++;
      // don't break here since GGA also needs LDA ptrs
    case XC_FAMILY_LDA:
      inlist.p[0].p = DOUBLEP(py_n_sg);
      inlist.p[0].special = N_SG;
      inlist.p[0].spinsize = 2;
      inlist.num += 1;
      outlist.p[0].p = DOUBLEP(py_e_g);
      outlist.p[0].special = E_G;
      outlist.p[0].spinsize = 1;
      outlist.p[1].p = DOUBLEP(py_dedn_sg);
      outlist.p[1].special = DEDN_SG;
      outlist.p[1].spinsize = 2;
      outlist.num += 2;
    }

  assert(inlist.num < MAXPTR);
  assert(outlist.num < MAXPTR);

  double *inblock[MAXPTR];
  double *outblock[MAXPTR];
  setupblockptrs(&info, &inlist, &outlist, &inblock[0], &outblock[0], blocksize);

  do {
    blocksize = blocksize<remaining ? blocksize : remaining;
    data2block(&info, &inlist, inblock, blocksize);
    double *n_sg = inblock[0];
    double *sigma_xg, *tau_sg;
    if (info.spinpolarized) {
      sigma_xg = inblock[1];
      tau_sg = inblock[2];
    } else {
      sigma_xg = inlist.p[1].p;
      tau_sg = inlist.p[2].p;
    }
    double *e_g = outblock[0];
    double *dedn_sg = outblock[1];
    double *dedsigma_xg = outblock[2];
    double *dedtau_sg = outblock[3];
    for (int i=0; i<2; i++) {
      if (self->functional[i] == NULL) continue;
      XC(func_type) *func = self->functional[i];
      int noutcopy=0;
      switch(func->info->family)
        {
        case XC_FAMILY_LDA:
          xc_lda_exc_vxc(func, blocksize, n_sg, e_g, dedn_sg);
          noutcopy = 2; // potentially decrease the size for block2dataadd if second functional less complex.
          break;
        case XC_FAMILY_HYB_GGA:
        case XC_FAMILY_GGA:
          xc_gga_exc_vxc(func, blocksize,
                         n_sg, sigma_xg, e_g,
                         dedn_sg, dedsigma_xg);
          noutcopy = 3; // potentially decrease the size for block2dataadd if second functional less complex.
          break;
        case XC_FAMILY_MGGA:
          xc_mgga_exc_vxc(func, blocksize, n_sg, sigma_xg, scratch_lapl, tau_sg,
                          e_g, dedn_sg, dedsigma_xg, scratch_vlapl, dedtau_sg);
          noutcopy = 4; // potentially decrease the size for block2dataadd if second functional less complex.
          break;
        }
      // if we have more than 1 functional, add results
      // canonical example: adding "x" results to "c"
      if (i==0)
        block2data(&info, &outblock[0], &outlist, n_sg, blocksize);
      else
        block2dataadd(&info, &outblock[0], &outlist, n_sg, blocksize, noutcopy);
    }

    for (int i=0; i<inlist.num; i++) inlist.p[i].p+=blocksize;
    for (int i=0; i<outlist.num; i++) outlist.p[i].p+=blocksize;

    remaining -= blocksize;
  } while (remaining>0);

  Py_RETURN_NONE;
}

static PyObject*
lxcXCFunctional_CalculateFXC(lxcXCFunctionalObject *self, PyObject *args)
{
  PyArrayObject* py_n_sg=NULL;
  PyArrayObject* py_v2rho2_xg=NULL;
  PyArrayObject* py_sigma_xg=NULL;
  PyArrayObject* py_v2rhosigma_yg=NULL;
  PyArrayObject* py_v2sigma2_yg=NULL;
  if (!PyArg_ParseTuple(args, "OO|OOO", &py_n_sg, &py_v2rho2_xg,
                        &py_sigma_xg, &py_v2rhosigma_yg, &py_v2sigma2_yg))
    return NULL;

  xcinfo info;
  info.nspin = self->nspin;
  info.spinpolarized = (info.nspin==2);
  info.ng = (info.spinpolarized) ? PyArray_DIMS(py_n_sg)[0]/2 : PyArray_DIMS(py_n_sg)[0];

  xcptrlist inlist,outlist;
  inlist.num=0;
  outlist.num=0;

  int blocksize = BLOCKSIZE;
  int remaining = info.ng;

  // setup pointers using most complex functional
  switch(self->functional[0]->info->family)
    {
    case XC_FAMILY_MGGA:
      // not supported
      assert(self->functional[0]->info->family != XC_FAMILY_MGGA);
      // don't break here since MGGA also needs GGA ptrs
    case XC_FAMILY_HYB_GGA:
    case XC_FAMILY_GGA:
      inlist.p[1].p = DOUBLEP(py_sigma_xg);
      inlist.p[1].special = 0;
      inlist.p[1].spinsize = 3;
      inlist.num++;
      outlist.p[1].p = DOUBLEP(py_v2rhosigma_yg);
      outlist.p[1].special = 0;
      outlist.p[1].spinsize = 6;
      outlist.p[2].p = DOUBLEP(py_v2sigma2_yg);
      outlist.p[2].special = 0;
      outlist.p[2].spinsize = 6;
      outlist.num+=2;
      // don't break here since GGA also needs LDA ptrs
    case XC_FAMILY_LDA:
      inlist.p[0].p = DOUBLEP(py_n_sg);
      inlist.p[0].special = N_SG;
      inlist.p[0].spinsize = 2;
      inlist.num += 1;
      outlist.p[0].p = DOUBLEP(py_v2rho2_xg);
      outlist.p[0].special = 0;
      outlist.p[0].spinsize = 3;
      outlist.num++;
    }

  assert(inlist.num < MAXPTR);
  assert(outlist.num < MAXPTR);

  double *inblock[MAXPTR];
  double *outblock[MAXPTR];
  setupblockptrs(&info, &inlist, &outlist, &inblock[0], &outblock[0], blocksize);

  do {
    blocksize = blocksize<remaining ? blocksize : remaining;
    data2block(&info, &inlist, inblock, blocksize);
    double *n_sg = inblock[0];
    double *sigma_xg;
    if (info.spinpolarized) {
      sigma_xg = inblock[1];
    } else {
      sigma_xg = inlist.p[1].p;
    }
    double *v2rho2 = outblock[0];
    double *v2rhosigma = outblock[1];
    double *v2sigma2 = outblock[2];
    for (int i=0; i<2; i++) {
      if (self->functional[i] == NULL) continue;
      XC(func_type) *func = self->functional[i];
      int noutcopy=0;
      switch(func->info->family)
        {
        case XC_FAMILY_LDA:
          xc_lda_fxc(func, blocksize, n_sg, v2rho2);
          noutcopy = 1; // potentially decrease the size for block2dataadd if second functional less complex.
          break;
        case XC_FAMILY_HYB_GGA:
        case XC_FAMILY_GGA:
          xc_gga_fxc(func, blocksize, n_sg, sigma_xg,
                     v2rho2, v2rhosigma, v2sigma2);
          noutcopy = 3; // potentially decrease the size for block2dataadd if second functional less complex.
          break;
        case XC_FAMILY_MGGA:
          // not supported by GPAW yet, so crash
          assert (func->info->family!=XC_FAMILY_MGGA);
          break;
        }
      // if we have more than 1 functional, add results
      // canonical example: adding "x" results to "c"
      if (i==0)
        block2data(&info, &outblock[0], &outlist, n_sg, blocksize);
      else
        block2dataadd(&info, &outblock[0], &outlist, n_sg, blocksize, noutcopy);
    }

    for (int i=0; i<inlist.num; i++) inlist.p[i].p+=blocksize;
    for (int i=0; i<outlist.num; i++) outlist.p[i].p+=blocksize;

    remaining -= blocksize;
  } while (remaining>0);

  Py_RETURN_NONE;
}

static PyMethodDef lxcXCFunctional_Methods[] = {
  {"is_gga",
   (PyCFunction)lxcXCFunctional_is_gga, METH_VARARGS, 0},
  {"is_mgga",
   (PyCFunction)lxcXCFunctional_is_mgga, METH_VARARGS, 0},
  {"calculate_fxc_fd_spinpaired",
   (PyCFunction)lxcXCFunctional_CalculateFXC_FD_SpinPaired, METH_VARARGS, 0},
  {"calculate",
   (PyCFunction)lxcXCFunctional_Calculate, METH_VARARGS, 0},
  {"calculate_fxc_spinpaired",
   (PyCFunction)lxcXCFunctional_CalculateFXC, METH_VARARGS, 0},
  {NULL, NULL, 0, NULL}
};

static PyObject* lxcXCFunctional_getattr(PyObject *obj, char *name)
{
  return Py_FindMethod(lxcXCFunctional_Methods, obj, name);
}

PyTypeObject lxcXCFunctionalType = {
  PyObject_HEAD_INIT(NULL)
  0,
  "lxcXCFunctional",
  sizeof(lxcXCFunctionalObject),
  0,
  (destructor)lxcXCFunctional_dealloc,
  0,
  lxcXCFunctional_getattr
};


PyObject * NewlxcXCFunctionalObject(PyObject *obj, PyObject *args)
{
  int xc, x, c; /* functionals identifier number */
  int nspin; /* XC_UNPOLARIZED or XC_POLARIZED  */

  if (!scratch) {
    scratch = (double*)malloc(LIBXCSCRATCHSIZE*sizeof(double));
    const int laplsize = BLOCKSIZE*sizeof(double)*2;
    scratch_lapl = (double*)malloc(laplsize);
    memset(scratch_lapl,0,laplsize);
    scratch_vlapl = (double*)malloc(laplsize);
  }

  if (!PyArg_ParseTuple(args, "iiii", &xc, &x, &c, &nspin)) {
    return NULL;
  }

  /* checking if the numbers xc x c are valid is done at python level */

  lxcXCFunctionalObject *self = PyObject_NEW(lxcXCFunctionalObject,
					     &lxcXCFunctionalType);

  if (self == NULL){
    return NULL;
  }

  assert(nspin==XC_UNPOLARIZED || nspin==XC_POLARIZED);
  self->nspin = nspin; /* must be common to x and c, so declared redundantly */

  int number,family,familyx,familyc;

  if (xc != -1) {
    xc_family_from_id(xc,&family,&number);
    assert (family != XC_FAMILY_UNKNOWN);
    XC(func_init)(&self->xc_functional, xc, nspin);
    self->functional[0]=&self->xc_functional;
    self->functional[1]=NULL;
  } else {
    assert (x!=-1 || c!=-1);

    if (x!=-1) {
      xc_family_from_id(x,&familyx,&number);
      assert (familyx != XC_FAMILY_UNKNOWN);
      XC(func_init)(&self->x_functional, x, nspin);
    }
    if (c!=-1) {
      xc_family_from_id(c,&familyc,&number);
      assert (familyc != XC_FAMILY_UNKNOWN);
      XC(func_init)(&self->c_functional, c, nspin);
    }

    if (x!=-1 && c!=-1) {
      /* put most complex functional first */
      /* important for later loops over functionals */
      if (familyx == XC_FAMILY_MGGA) {
        self->functional[0]=&self->x_functional;
        self->functional[1]=&self->c_functional;
      } else if (familyc == XC_FAMILY_MGGA) {
        self->functional[0]=&self->c_functional;
        self->functional[1]=&self->x_functional;
      } else if (familyx == XC_FAMILY_GGA || familyx == XC_FAMILY_HYB_GGA) {
        self->functional[0]=&self->x_functional;
        self->functional[1]=&self->c_functional;
      } else {
        // either c is GGA, or both are LDA (so don't care)
      self->functional[0]=&self->c_functional;
      self->functional[1]=&self->x_functional;
      }
    } else if (x!=-1) {
      self->functional[0]=&self->x_functional;
      self->functional[1]=NULL;
    } else if (c!=-1) {
      self->functional[0]=&self->c_functional;
      self->functional[1]=NULL;
    }
  }

  return (PyObject*)self;
}

PyObject * lxcXCFuncNum(PyObject *obj, PyObject *args)
{
  char *funcname;
  if (!PyArg_ParseTuple(args, "s", &funcname)) {
    return NULL;
  }

  int num = XC(functional_get_number)(funcname);
  if (num != -1)
    return Py_BuildValue("i",num);
  else
    Py_RETURN_NONE;
}
