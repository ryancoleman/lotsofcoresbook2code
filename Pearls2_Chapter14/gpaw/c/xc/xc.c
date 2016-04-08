/*  Copyright (C) 2003-2007  CAMP
 *  Copyright (C) 2007-2009  CAMd
 *  Please see the accompanying LICENSE file for further information. */

#include <Python.h>
#define PY_ARRAY_UNIQUE_SYMBOL GPAW_ARRAY_API
#define NO_IMPORT_ARRAY
#include <numpy/arrayobject.h>
#include "xc_gpaw.h"
#include "../extensions.h"

//
//          __  2
// a2    = |\/n|
//
//         dE
// dedrs = ---
//         dr
//           s
//
//            dE
// deda2 = ---------
//            __  2
//         d(|\/n| )
//

void init_mgga(void** params, int code, int nspin);
void calc_mgga(void** params, int nspin, int ng,
               const double* n_g, const double* sigma_g, const double* tau_g,
               double *e_g, double *v_g, double *dedsigma_g, double *dedtau_g);

double pbe_exchange(const xc_parameters* par,
                    double n, double rs, double a2,
                    double* dedrs, double* deda2);
double pbe_correlation(double n, double rs, double zeta, double a2,
                       bool gga, bool spinpol,
                       double* dedrs, double* dedzeta, double* deda2);
double pw91_exchange(const xc_parameters* par,
                     double n, double rs, double a2,
                     double* dedrs, double* deda2);
double pw91_correlation(double n, double rs, double zeta, double a2,
                        bool gga, bool spinpol,
                        double* dedrs, double* dedzeta, double* deda2);
double rpbe_exchange(const xc_parameters* par,
                     double n, double rs, double a2,
                     double* dedrs, double* deda2);
double beefvdw_exchange(const xc_parameters* par,
                        double n, double rs, double a2,
                        double* dedrs, double* deda2);

//
typedef struct
{
  PyObject_HEAD
  double (*exchange)(const xc_parameters* par,
                     double n, double rs, double a2,
                     double* dedrs, double* deda2);
  double (*correlation)(double n, double rs, double zeta, double a2,
                        bool gga, bool spinpol,
                        double* dedrs, double* dedzeta, double* deda2);
  xc_parameters par;
  // below added by cpo for mgga functionals outside of libxc (TPSS, M06L, etc.)
  void* mgga;
} XCFunctionalObject;

static void XCFunctional_dealloc(XCFunctionalObject *self)
{
  PyObject_DEL(self);
}

static PyObject*
XCFunctional_calculate(XCFunctionalObject *self, PyObject *args)
{
  PyArrayObject* e_array;
  PyArrayObject* n_array;
  PyArrayObject* v_array;
  PyArrayObject* sigma_array = 0;
  PyArrayObject* dedsigma_array = 0;
  PyArrayObject* tau_array = 0;
  PyArrayObject* dedtau_array = 0;

  if (!PyArg_ParseTuple(args, "OOO|OOOO", &e_array, &n_array, &v_array,
                        &sigma_array, &dedsigma_array, &tau_array, &dedtau_array))
    return NULL;

  int ng = 1;
  for (int d = 0; d < PyArray_NDIM(e_array); d++)
    ng *= PyArray_DIM(e_array, d);

  xc_parameters* par = &self->par;

  double* e_g = DOUBLEP(e_array);
  const double* n_g = DOUBLEP(n_array);
  double* v_g = DOUBLEP(v_array);

  const double* sigma_g = 0;
  double* dedsigma_g = 0;
  if (par->gga)
    {
      sigma_g = DOUBLEP(sigma_array);
      dedsigma_g = DOUBLEP(dedsigma_array);
    }

  const double* tau_g = 0;
  double* dedtau_g = 0;
  if (self->mgga)
    {
      tau_g = DOUBLEP(tau_array);
      dedtau_g = DOUBLEP(dedtau_array);
    }

  if (self->mgga) {
    int nspin = PyArray_DIM(n_array, 0) == 1 ? 1 : 2;
    calc_mgga(&self->mgga, nspin, ng, n_g, sigma_g, tau_g, e_g, v_g, dedsigma_g, dedtau_g);
    Py_RETURN_NONE;
  }

  if (PyArray_DIM(n_array, 0) == 1)
    for (int g = 0; g < ng; g++)
      {
        double n = n_g[g];
        if (n < NMIN)
          n = NMIN;
        double rs = pow(C0I / n, THIRD);
        double dexdrs;
        double dexda2;
        double ex;
        double decdrs;
        double decda2;
        double ec;
        if (par->gga)
          {
            double a2 = sigma_g[g];
            ex = self->exchange(par, n, rs, a2, &dexdrs, &dexda2);
            ec = self->correlation(n, rs, 0.0, a2, 1, 0, &decdrs, 0, &decda2);
            dedsigma_g[g] = n * (dexda2 + decda2);
          }
        else
          {
            ex = self->exchange(par, n, rs, 0.0, &dexdrs, 0);
            ec = self->correlation(n, rs, 0.0, 0.0, 0, 0, &decdrs, 0, 0);
          }
        e_g[g] = n * (ex + ec);
        v_g[g] += ex + ec - rs * (dexdrs + decdrs) / 3.0;
      }
  else
    {
      const double* na_g = n_g;
      double* va_g = v_g;
      const double* nb_g = na_g + ng;
      double* vb_g = va_g + ng;

      const double* sigma0_g = 0;
      const double* sigma1_g = 0;
      const double* sigma2_g = 0;
      double* dedsigma0_g = 0;
      double* dedsigma1_g = 0;
      double* dedsigma2_g = 0;

      const xc_parameters* par = &self->par;
      if (par->gga)
        {
          sigma0_g = sigma_g;
          sigma1_g = sigma0_g + ng;
          sigma2_g = sigma1_g + ng;
          dedsigma0_g = dedsigma_g;
          dedsigma1_g = dedsigma0_g + ng;
          dedsigma2_g = dedsigma1_g + ng;
        }

      for (int g = 0; g < ng; g++)
        {
          double na = 2.0 * na_g[g];
          if (na < NMIN)
            na = NMIN;
          double rsa = pow(C0I / na, THIRD);
          double nb = 2.0 * nb_g[g];
          if (nb < NMIN)
            nb = NMIN;
          double rsb = pow(C0I / nb, THIRD);
          double n = 0.5 * (na + nb);
          double rs = pow(C0I / n, THIRD);
          double zeta = 0.5 * (na - nb) / n;
          double dexadrs;
          double dexada2;
          double exa;
          double dexbdrs;
          double dexbda2;
          double exb;
          double decdrs;
          double decdzeta;
          double decda2;
          double ec;
          if (par->gga)
            {
              exa = self->exchange(par, na, rsa, 4.0 * sigma0_g[g],
                                   &dexadrs, &dexada2);
              exb = self->exchange(par, nb, rsb, 4.0 * sigma2_g[g],
                                   &dexbdrs, &dexbda2);
              double a2 = sigma0_g[g] + 2 * sigma1_g[g] + sigma2_g[g];
              ec = self->correlation(n, rs, zeta, a2, 1, 1,
                                     &decdrs, &decdzeta, &decda2);
              dedsigma0_g[g] = 2 * na * dexada2 + n * decda2;
              dedsigma1_g[g] = 2 * n * decda2;
              dedsigma2_g[g] = 2 * nb * dexbda2 + n * decda2;
            }
          else
            {
              exa = self->exchange(par, na, rsa, 0.0, &dexadrs, 0);
              exb = self->exchange(par, nb, rsb, 0.0, &dexbdrs, 0);
              ec = self->correlation(n, rs, zeta, 0.0, 0, 1,
                                     &decdrs, &decdzeta, 0);
            }
          e_g[g] = 0.5 * (na * exa + nb * exb) + n * ec;
          va_g[g] += (exa + ec -
                      (rsa * dexadrs + rs * decdrs) / 3.0 -
                      (zeta - 1.0) * decdzeta);
          vb_g[g] += (exb + ec -
                      (rsb * dexbdrs + rs * decdrs) / 3.0 -
                      (zeta + 1.0) * decdzeta);
        }
    }
  Py_RETURN_NONE;
}

static PyMethodDef XCFunctional_Methods[] = {
    {"calculate",
     (PyCFunction)XCFunctional_calculate, METH_VARARGS, 0},
    {NULL, NULL, 0, NULL}
};

static PyObject* XCFunctional_getattr(PyObject *obj, char *name)
{
    return Py_FindMethod(XCFunctional_Methods, obj, name);
}

PyTypeObject XCFunctionalType = {
  PyObject_HEAD_INIT(NULL)
  0,
  "XCFunctional",
  sizeof(XCFunctionalObject),
  0,
  (destructor)XCFunctional_dealloc,
  0,
  XCFunctional_getattr
};

PyObject * NewXCFunctionalObject(PyObject *obj, PyObject *args)
{
  int code;
  PyArrayObject* parameters = 0;
  if (!PyArg_ParseTuple(args, "i|O", &code, &parameters))
    return NULL;

  XCFunctionalObject *self = PyObject_NEW(XCFunctionalObject,
                                          &XCFunctionalType);
  if (self == NULL)
    return NULL;

  self->mgga = NULL;

  self->par.gga = 1;

  self->correlation = pbe_correlation;
  self->exchange = pbe_exchange;

  if (code == -1) {
    // LDA
    self->par.gga = 0;
  }
  else if (code == 0) {
    // PBE
    self->par.kappa = 0.804;
  }
  else if (code == 1) {
    // revPBE
    self->par.kappa = 1.245;
  }
  else if (code == 2) {
    // RPBE
    self->exchange = rpbe_exchange;
  }
  else if (code == 14) {
    // PW91
    self->exchange = pw91_exchange;
  }
  else if (code == 20 || code == 21 || code == 22) {
    // MGGA
    const int nspin = 1; // a guess, perhaps corrected later in calc_mgga
    init_mgga(&self->mgga,code,nspin);
  }
  else {
    assert (code == 17);
    // BEEF-vdW
    self->exchange = beefvdw_exchange;
    int n = PyArray_DIM(parameters, 0);
    assert(n <= 110);
    double* p = (double*)PyArray_BYTES(parameters);
    for (int i = 0; i < n; i++)
      self->par.parameters[i] = p[i];
    self->par.nparameters = n / 2;
  }
  return (PyObject*)self;
}
