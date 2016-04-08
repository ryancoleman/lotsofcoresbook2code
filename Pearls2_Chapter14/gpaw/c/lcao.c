/*  Copyright (C) 2003-2007  CAMP
 *  Copyright (C) 2007-2009  CAMd
 *  Please see the accompanying LICENSE file for further information. */

#include <complex.h>
#include "localized_functions.h"
#include "extensions.h"
#include "bmgs/bmgs.h"

#ifdef GPAW_NO_UNDERSCORE_BLAS
#  define dgemv_ dgemv
#  define dgemm_ dgemm
#endif

int dgemv_(char *trans, int *m, int * n,
	   double *alpha, double *a, int *lda,
	   double *x, int *incx, double *beta,
	   double *y, int *incy);
int dgemm_(char *transa, char *transb, int *m, int * n,
	   int *k, const double *alpha, double *a, int *lda,
	   double *b, int *ldb, double *beta,
	   double *c, int *ldc);


//                    +-----------n
//  +----m   +----m   | +----c+m  |
//  |    |   |    |   | |    |    |
//  |  b | = |  v | * | |  a |    |
//  |    |   |    |   | |    |    |
//  0----+   0----+   | c----+    |
//                    |           |
//                    0-----------+
void cut(const double* a, const int n[3], const int c[3],
	 const double* v,
	 double* b, const int m[3])
{
  a += c[2] + (c[1] + c[0] * n[1]) * n[2];
  for (int i0 = 0; i0 < m[0]; i0++)
    {
      for (int i1 = 0; i1 < m[1]; i1++)
        {
	  for (int i2 = 0; i2 < m[2]; i2++)
	    b[i2] = v[i2] * a[i2];
          a += n[2];
          b += m[2];
          v += m[2];
        }
      a += n[2] * (n[1] - m[1]);
    }
}


PyObject * overlap(PyObject* self, PyObject *args)
{
  PyObject* lfs_b_obj;
  PyArrayObject* m_b_obj;
  PyArrayObject* phase_bk_obj;
  PyArrayObject* vt_sG_obj;
  PyArrayObject* Vt_skmm_obj;
  if (!PyArg_ParseTuple(args, "OOOOO", &lfs_b_obj, &m_b_obj, &phase_bk_obj,
			&vt_sG_obj, &Vt_skmm_obj))
    return NULL;

  int nk = PyArray_DIMS(phase_bk_obj)[1];
  int nm = PyArray_DIMS(Vt_skmm_obj)[2];
  int nspins = PyArray_DIMS(vt_sG_obj)[0];

  const long *m_b = LONGP(m_b_obj);
  const double complex *phase_bk = COMPLEXP(phase_bk_obj);
  const double *vt_sG = DOUBLEP(vt_sG_obj);
  double *Vt_smm = 0;
  double complex *Vt_skmm = 0;

  if (nk == 0)
    Vt_smm = DOUBLEP(Vt_skmm_obj);
  else
    Vt_skmm = COMPLEXP(Vt_skmm_obj);

  int nb = PyList_Size(lfs_b_obj);

  int nmem = 0;
  double* a1 = 0;
  for (int b1 = 0; b1 < nb; b1++)
    {
      const LocalizedFunctionsObject* lf1 =
	(const LocalizedFunctionsObject*)PyList_GetItem(lfs_b_obj, b1);
      int m1 = m_b[b1];
      int nao1 = lf1->nf;
      double* f1 = lf1->f;
      double* vt1 = GPAW_MALLOC(double, lf1->ng0 * nspins);
      for (int s = 0; s < nspins; s++)
	bmgs_cut(vt_sG + s * lf1->ng, lf1->size, lf1->start,
		 vt1 + s * lf1->ng0, lf1->size0);
      for (int b2 = b1; b2 < nb; b2++)
	{
	  const LocalizedFunctionsObject* lf2 =
	    (const LocalizedFunctionsObject*)PyList_GetItem(lfs_b_obj, b2);
	  int beg[3];
	  int end[3];
	  int size[3];
	  int beg1[3];
	  int beg2[3];
	  bool overlap = true;
	  for (int c = 0; c < 3; c++)
	    {
	      beg[c] = MAX(lf1->start[c], lf2->start[c]);
	      end[c] = MIN(lf1->start[c] + lf1->size0[c],
			   lf2->start[c] + lf2->size0[c]);
	      size[c] = end[c] - beg[c];
	      if (size[c] <= 0)
		{
		  overlap = false;
		  continue;
		}
	      beg1[c] = beg[c] - lf1->start[c];
	      beg2[c] = beg[c] - lf2->start[c];
	    }
	  int nao2 = lf2->nf;
	  if (overlap)
	    {
	      int ng = size[0] * size[1] * size[2];
	      int n = ng * (nao1 + nao2) + nao1 * nao2;
	      if (n > nmem)
		{
		  if (nmem != 0)
		    free(a1);
		  nmem = n;
		  a1 = GPAW_MALLOC(double, nmem);
		}
	      double* a2 = a1 + ng * nao1;
	      double* H = a2 + ng * nao2;
	      double* f2 = lf2->f;
	      double* vt2 = lf2->w;
	      double dv = lf1->dv;
	      int m2 = m_b[b2];
	      if (b2 > b1)
		for (int i = 0; i < nao2; i++)
		  bmgs_cut(f2 + i * lf2->ng0, lf2->size0, beg2,
			   a2 + i * ng, size);
	      else
		  a2 = f2;
	      for (int s = 0; s < nspins; s++)
		{
		  if (b2 > b1)
		    {
		      bmgs_cut(vt1 + s * lf1->ng0, lf1->size0, beg1, vt2, size);
		      for (int i = 0; i < nao1; i++)
			cut(f1 + i * lf1->ng0, lf1->size0, beg1, vt2,
			    a1 + i * ng, size);
		    }
		  else
		    {
		      for (int i1 = 0; i1 < nao1; i1++)
			for (int g = 0; g < ng; g++)
			  a1[i1 * ng + g] = (vt1[g + s * lf1->ng0] *
					     f1[i1 * ng + g]);
		    }
		  double zero = 0.0;
		  dgemm_("t", "n", &nao2, &nao1, &ng, &dv,
			 a2, &ng, a1, &ng, &zero, H, &nao2);
		  if (nk == 0)
		    {
		      double* Vt_mm = (Vt_smm + s * nm * nm + m1 + m2 * nm);
		      if (b2 == b1)
			for (int i1 = 0; i1 < nao1; i1++)
			  for (int i2 = i1; i2 < nao2; i2++)
			    Vt_mm[i1 + i2 * nm] += H[i2 + i1 * nao2];
		      else if (m1 == m2)
			for (int i1 = 0; i1 < nao1; i1++)
			  for (int i2 = i1; i2 < nao2; i2++)
			    Vt_mm[i1 + i2 * nm] += (H[i2 + i1 * nao2] +
						    H[i1 + i2 * nao2]);
		      else
			for (int ii = 0, i1 = 0; i1 < nao1; i1++)
			  for (int i2 = 0; i2 < nao2; i2++, ii++)
			    Vt_mm[i1 + i2 * nm] += H[ii];
		    }
		  else
		    for (int k = 0; k < nk; k++)
		      {
			double complex* Vt_mm = (Vt_skmm +
						 (s * nk + k) * nm * nm +
						 m1 + m2 * nm);
			if (b2 == b1)
			  for (int i1 = 0; i1 < nao1; i1++)
			    for (int i2 = i1; i2 < nao2; i2++)
			      Vt_mm[i1 + i2 * nm] += H[i2 + i1 * nao2];
			else
			  {
			    double complex phase = \
			      (phase_bk[b1 * nk + k] *
			       conj(phase_bk[b2 * nk + k]));
			    if (m1 == m2)
			      for (int i1 = 0; i1 < nao1; i1++)
				for (int i2 = i1; i2 < nao2; i2++)
				  Vt_mm[i1 + i2 * nm] += \
				    (phase * H[i2 + i1 * nao2] +
				     conj(phase) * H[i1 + i2 * nao2]);
			    else
			      for (int ii = 0, i1 = 0; i1 < nao1; i1++)
				for (int i2 = 0; i2 < nao2; i2++, ii++)
				  Vt_mm[i1 + i2 * nm] += phase * H[ii];
			  }
		      }
		}
	    }
	}
      free(vt1);
    }
  if (nmem != 0)
    free(a1);
  Py_RETURN_NONE;
}
