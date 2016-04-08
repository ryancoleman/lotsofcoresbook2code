/*  Copyright (C) 2003-2007  CAMP
 *  Copyright (C) 2007-2009  CAMd
 *  Copyright (C) 2005-2007  CSC - IT Center for Science Ltd.
 *  Please see the accompanying LICENSE file for further information. */

#include <Python.h>
#define PY_ARRAY_UNIQUE_SYMBOL GPAW_ARRAY_API
#define NO_IMPORT_ARRAY
#include <numpy/arrayobject.h>
#include "extensions.h"
#ifdef GPAW_NO_UNDERSCORE_LAPACK
#  define dlamch_ dlamch
#  define dsyev_  dsyev
#  define zheev_  zheev
#  define dsyevr_ dsyevr
#  define zheevr_ zheevr 
#  define dsygv_  dsygv
#  define dsygvx_ dsygvx
#  define dhegv_  dhegv
#  define zhegv_  zhegv
#  define zhegvx_ zhegvx
#  define dgeev_  dgeev
#  define dpotrf_ dpotrf
#  define dpotri_ dpotri
#  define zpotrf_ zpotrf
#  define zpotri_ zpotri
#  define dtrtri_ dtrtri
#  define ztrtri_ ztrtri
#  define dsytrf_ dsytrf
#  define zsytrf_ zsytrf
#  define dgetrf_ dgetrf
#  define zgetrf_ zgetrf
#  define dsytri_ dsytri
#  define zsytri_ zsytri
#  define dgetri_ dgetri
#  define zgetri_ zgetri
#  define zgbsv_  zgbsv
#  define zgttrf_ zgttrf
#  define zgttrs_ zgttrs
#  define ilaenv_ ilaenv
#endif

double dlamch_(char* cmach);

void dsyev_(char *jobz, char *uplo, int *n, 
            double *a, int *lda, double *w, double *work, int *lwork,
            int *info);
void zheev_(char *jobz, char *uplo, int *n,
            void *a, int *lda, double *w, void *work,
            int *lwork, double *rwork, int *lrwork, int *info);
void dsyevr_(char *jobz, char *range, char *uplo, int *n, 
             double *a, int *lda, 
             double *vl, double *vu, int *il, int*iu, double *abstol,
             int *m, double *w, double *z, int *ldz, int *isuppz, 
             double *work, int *lwork, int *iwork, int *liwork,
             int *info);
void zheevr_(char *jobz, char *range, char *uplo, int *n, 
             void *a, int *lda, 
             double *vl, double *vu, int *il, int *iu, double *abstol,
             int *m, double *w, void *z, int *ldz, int *isuppz,
             void *work, int *lwork, double *rwork, int *lrwork, 
             int *iwork, int *liwork,
             int *info);           
void dsygv_(int *itype, char *jobz, char *uplo, int *
           n, double *a, int *lda, double *b, int *ldb,
           double *w, double *work, int *lwork, int *info);
void dsygvx_(int *itype, char *jobz, char *range, char *uplo,
             int *n, void *a, int *lda, void *b, int *ldb,
             double *vl, double *vu, int *il, int *iu, double *abstol,
             int *m, double *w, void *z, int *ldz, void *work,
             int *lwork, int *iwork, int *ifail,
             int *info);
void zhegv_(int *itype, char *jobz, char *uplo, int *
           n, void *a, int *lda, void *b, int *ldb,
           double *w, void *work, int *lwork,
           double *rwork,
           int *lrwork, int *info);
void zhegvx_(int *itype, char *jobz, char *range, char *uplo, 
             int *n, void *a, int *lda, void *b, int *ldb,
             double *vl, double *vu, int *il, int *iu, double *abstol,
             int *m, double *w, void *z, int *ldz, void *work, 
             int *lwork, double *rwork, int *iwork, int *ifail, 
             int *info);
void dpotrf_(char *uplo, int *n, double *a, int *
            lda, int *info);
void dpotri_(char *uplo, int *n, double *a, int *
            lda, int *info);
void zpotrf_(char *uplo, int *n, void *a,
            int *lda, int *info);
void zpotri_(char *uplo, int *n, void *a,
            int *lda, int *info);
void dgeev_(char *jovl, char *jobvr, int *n, double *a, int *lda,
           double *wr, double *wl,
           double *vl, int *ldvl, double *vr, int *ldvr,
           double *work, int *lwork, int *info);

void dtrtri_(char *uplo,char *diag, int *n, void *a,
            int *lda, int *info );
void ztrtri_(char *uplo,char *diag, int *n, void *a,
             int *lda, int *info );

void dsytrf_(char *uplo, int *n, double *a, int *lda, int *ipiv,
            double *work, int *lwork, int *info);
void zsytrf_(char *uplo, int *n, void *a, int *lda, int *ipiv,
            void *work, int *lwork, int *info);

void dgetrf_(int *n, int *m, double *a, int *lda, int *ipiv, int *info);
void zgetrf_(int *n, int *m, void *a, int *lda, int *ipiv, int *info);

void dsytri_(char *uplo, int *n, double *a, int *lda, int *ipiv,
             double *work, int *info);
void zsytri_(char *uplo, int *n, void *a, int *lda, int *ipiv,
             void *work, int *info);

void dgetri_(int *n, double *a, int *lda, int *ipiv,
             double *work, int *lwork, int *info);
void zgetri_(int *n, void *a, int *lda, int *ipiv,
             void *work, int *lwork, int *info);
void zgbsv_(int*n, int* kl, int* ku, int* nrhs, void* ab, int*ldab,
            int*ipiv, void* b, int*ldb, int*info);
void zgttrf_(int* n, void* dl, void* d, void* du,
            void* du2, int* ipiv, int* info);
void zgttrs_(char* tran, int* n, int* nrhs, void* dl,
               void* d, void* du, void* du2,
               int* ipiv, void* b, int* ldb, int* info);
int ilaenv_(int* ispec, char* name, char* opts, int* n1, 
            int* n2, int* n3, int* n4);

PyObject* diagonalize(PyObject *self, PyObject *args)
{
  PyArrayObject* a;
  PyArrayObject* w;
  if (!PyArg_ParseTuple(args, "OO", &a, &w))
    return NULL;
  int n = PyArray_DIMS(a)[0];
  int lda = n;
  int info = 0;
  if (PyArray_DESCR(a)->type_num == NPY_DOUBLE)
    {
      int lwork = 3 * n + 1;
      double* work = GPAW_MALLOC(double, lwork);
      dsyev_("V", "U", &n, DOUBLEP(a), &lda,
             DOUBLEP(w), work, &lwork, &info);
      free(work);
    }
  else
    {
      int lwork = 2 * n + 1;
      int lrwork = 3 * n + 1;
      void* work = GPAW_MALLOC(double_complex, lwork);
      double* rwork = GPAW_MALLOC(double, lrwork);
      zheev_("V", "U", &n, (void*)COMPLEXP(a), &lda,
             DOUBLEP(w),
             work, &lwork, rwork, &lrwork, &info);
      free(work);
      free(rwork);
    }
  return Py_BuildValue("i", info);
}

PyObject* diagonalize_mr3(PyObject *self, PyObject *args)
{
  PyArrayObject* a;
  PyArrayObject* w;
  PyArrayObject* z;
  if (!PyArg_ParseTuple(args, "OOO", &a, &w, &z))
    return NULL;
  char jobz = 'V';
  char range = 'A';
  char uplo = 'U';
  int n = PyArray_DIMS(a)[0];
  int lda = MAX(1, n);
  double vl, vu;
  int il, iu;
  double abstol = dlamch_("Safe minimum");
  int m = n; /* assume we find all eigenvalues */ 
  int ldz = lda;
  int info = 0;
  int* isuppz = GPAW_MALLOC(int, 2*m);
  if (PyArray_DESCR(a)->type_num == NPY_DOUBLE)
    {
      /* Minimum workspace plus a little extra */
      int lwork = 26 * n + 1;
      int liwork = 10 * n + 1;
      double* work = GPAW_MALLOC(double, lwork);
      int* iwork = GPAW_MALLOC(int, liwork);
      dsyevr_(&jobz, &range, &uplo, &n, 
              DOUBLEP(a), &lda,
              &vl, &vu, &il, &iu, &abstol, 
              &m, DOUBLEP(w), DOUBLEP(z), &ldz, isuppz, 
              work, &lwork, iwork, &liwork,
              &info);
      free(work);
      free(iwork);
    }
  else
    {
      /* Minimum workspace plus a little extra */
      int lwork = 2 * n + 1;
      int lrwork = 24 * n + 1;
      int liwork = 10 * n + 1;
      void* work = GPAW_MALLOC(double_complex, lwork);
      double* rwork = GPAW_MALLOC(double, lrwork);
      int* iwork = GPAW_MALLOC(int, liwork);
      zheevr_(&jobz, &range, &uplo, &n, 
              (void*)COMPLEXP(a), &lda,
              &vl, &vu, &il, &iu, &abstol, 
              &m,  DOUBLEP(w), (void*)COMPLEXP(z), &ldz, isuppz, 
              work, &lwork, rwork, &lrwork,
              iwork, &liwork, 
              &info);
      free(work);
      free(rwork);
      free(iwork);
    }
  free(isuppz);
  // If this fails, fewer eigenvalues than request were computed
  assert (m == n);
  return Py_BuildValue("i", info);
}

PyObject* general_diagonalize(PyObject *self, PyObject *args)
{
  PyArrayObject* a;
  PyArrayObject* w;
  PyArrayObject* b;
  PyArrayObject* z;
  int iu = -1;
  if (!PyArg_ParseTuple(args, "OOO|Oi", &a, &w, &b, &z, &iu))
    return NULL;
  int itype = 1;
  char jobz = 'V';
  char range = 'I';
  char uplo = 'U';
  int n = PyArray_DIMS(a)[0];
  int lda = MAX(1, n);
  int ldb = lda;
  double vl, vu;
  int il = 1;
  double abstol = dlamch_("Safe minimum");
  int m; 
  int ldz = lda;
  int info = 0;
  int ispec = 1;
  int dummy = -1;
  int NB = ilaenv_(&ispec, "dsytrd", &uplo, &n, &dummy, &dummy, &dummy);

  if (PyArray_DESCR(a)->type_num == NPY_DOUBLE)
    {
      if (iu == -1)
        {
          int lwork = MAX((NB + 2) * n, 3 * n + 1);
          double* work = GPAW_MALLOC(double, lwork);
          dsygv_(&itype, &jobz, &uplo, &n, DOUBLEP(a), &lda,
                 DOUBLEP(b), &ldb, DOUBLEP(w),
                 work, &lwork, &info);
          free(work);
        }
      else
        {
          int lwork = MAX((NB + 3) * n, 8 * n);
          int liwork = 5 * n;
          double* work = GPAW_MALLOC(double, lwork);
          int* iwork = GPAW_MALLOC(int, liwork);
          int* ifail = GPAW_MALLOC(int, n);
          dsygvx_(&itype, &jobz, &range, &uplo, &n, 
                  DOUBLEP(a), &lda, DOUBLEP(b), &ldb,
                  &vl, &vu, &il, &iu, &abstol, 
                  &m, DOUBLEP(w), DOUBLEP(z), &ldz, 
                  work, &lwork, iwork, ifail,
                  &info);
          free(iwork);
          free(work);
          free(ifail);
          assert (m == iu);
        }
    }
  else
    {
      if (iu == -1)
        {
          int lwork = MAX((NB + 1) * n, 2 * n + 1);
          int lrwork = MAX(1, 3 * n + 1);
          void* work = GPAW_MALLOC(double_complex, lwork);
          double* rwork = GPAW_MALLOC(double, lrwork);
          zhegv_(&itype, &jobz, &uplo, &n, (void*)COMPLEXP(a), &lda,
                 (void*)COMPLEXP(b), &lda,
                 DOUBLEP(w),
                 work, &lwork, rwork, &lrwork, &info);
          free(work);
          free(rwork);
        }
      else
        {
          int lwork = MAX((NB + 1) * n, 2 * n);
          int lrwork = 7 * n;
          int liwork = 5 * n;
          void* work = GPAW_MALLOC(double_complex, lwork);
          double* rwork = GPAW_MALLOC(double, lrwork);
          int* iwork = GPAW_MALLOC(int, liwork);
          int* ifail = GPAW_MALLOC(int, n);
          zhegvx_(&itype, &jobz, &range, &uplo, &n,
                  (void*)COMPLEXP(a), &lda, (void*)COMPLEXP(b), &ldb,
                  &vl, &vu, &il, &iu, &abstol,
                  &m,  DOUBLEP(w), (void*)COMPLEXP(z), &ldz,
                  work, &lwork, rwork, iwork, ifail, &info);
          free(work);
          free(rwork);
          free(iwork);
          free(ifail);
          assert (m == iu);
        }
    }
  return Py_BuildValue("i", info);
}

PyObject* inverse_cholesky(PyObject *self, PyObject *args)
{
  PyArrayObject* a;
  if (!PyArg_ParseTuple(args, "O", &a))
    return NULL;
  int n = PyArray_DIMS(a)[0];
  int lda = MAX(1, n);
  int info = 0;

  if (PyArray_DESCR(a)->type_num == NPY_DOUBLE)
    {
      dpotrf_("U", &n, (void*)DOUBLEP(a), &lda, &info);
      if (info == 0)
        {
          dtrtri_("U", "N", &n, (void*)DOUBLEP(a), &lda, &info);
          if (info == 0)
            {
              /* Make sure that the other diagonal is zero */
              double* ap = DOUBLEP(a);
              ap++;
              for (int i = 0; i < n - 1; i++)
                {
                  memset(ap, 0, (n-1-i) * sizeof(double));
                  ap += n + 1;
                }
            }
        }
    }
  else
    {
      zpotrf_("U", &n, (void*)COMPLEXP(a), &lda, &info);
      if (info == 0)
        {
          ztrtri_("U", "N", &n, (void*)DOUBLEP(a), &lda, &info);
          if (info == 0)
            {
              /* Make sure that lower diagonal is zero */
              double_complex* ap = COMPLEXP(a);
              ap++;
              for (int i = 0; i < n - 1; i++)
                {
                  memset(ap, 0, (n-1-i) * sizeof(double_complex));
                  ap += n + 1;
                }
            }
        }
    }
  return Py_BuildValue("i", info);
}

void swap(double *a, double *b) {
  double tmp=*b;
  *b = *a;
  *a = tmp;
}
void transpose(double *A, int n) {
  int i, j;
  int in=0;
  for(i=0;i<n-1;i++) {
    for(j=i+1;j<n;j++)
      swap(A+in+j,A+j*n+i);
    in+=n;
  }
}
void print(double *A, int n) {
  int i,j;
  for(i=0;i<n;i++) {
    if(i) printf(" (");
    else printf("((");
    for(j=0;j<n;j++) {
      printf(" %g",A[n*i+j]);
    }
    if(i==n-1) printf("))\n");
    else printf(")\n");
  }
}
PyObject* right_eigenvectors(PyObject *self, PyObject *args)
/* Return eigenvalues and right eigenvectors of a
 * nonsymmetric eigenvalue problem
 */
{
  PyArrayObject* A;
  PyArrayObject* v; /* eigenvectors */
  PyArrayObject* w; /* eigenvalues */
  if (!PyArg_ParseTuple(args, "OOO", &A, &w, &v))
    return NULL;
  int n = PyArray_DIMS(A)[0];
  int lda = n;
  int info = 0;
  if (PyArray_DESCR(A)->type_num == NPY_DOUBLE)
    {
      int lwork = -1;
      double* work = GPAW_MALLOC(double, 1);
      double* wr = GPAW_MALLOC(double, n);
      double* wi = GPAW_MALLOC(double, n);
      int ldvl = 1;
      int ldvr = n;
      double* vl = 0;
      int i;
      /* get size of work needed */
      dgeev_("No eigenvectors left", "Vectors right",
             &n, DOUBLEP(A), &lda, wr, wi,
             vl, &ldvl, DOUBLEP(v), &ldvr, work, &lwork, &info);
      lwork = (int) work[0];
      free(work);
      work = GPAW_MALLOC(double, lwork);

      transpose(DOUBLEP(A),n); /* transform to Fortran form */
      dgeev_("No eigenvectors left", "Vectors right",
             &n, DOUBLEP(A), &lda, wr, wi,
             vl, &ldvl, DOUBLEP(v), &ldvr, work, &lwork, &info);

      for(i=0;i<n;i++) {
        if(wi[i]!=0.)
          printf("<diagonalize_nonsymmetric> dgeev i=%d,wi[i]=%g\n",
                 i,wi[i]);
        DOUBLEP(w)[i]=wr[i];
      }
      free(wr);
      free(wi);
      free(work);
    }
  return Py_BuildValue("i", info);
}

PyObject* inverse_general(PyObject *self, PyObject *args)
{
  PyArrayObject* a;
  if (!PyArg_ParseTuple(args, "O", &a))
    return NULL;
  int n = PyArray_DIMS(a)[0];
  int m = n;
  int lda = n;
  int lwork = n;
  int* ipiv = GPAW_MALLOC(int, n);
  int info = 0;
  if (PyArray_DESCR(a)->type_num == NPY_DOUBLE)
    {
      double* work = GPAW_MALLOC(double, lwork);
      dgetrf_(&n, &m, DOUBLEP(a), &lda, ipiv, &info);
      dgetri_(&n, DOUBLEP(a), &lda, ipiv, work, &lwork, &info);
      free(work);
    }
  else
    {
      void *work = GPAW_MALLOC(double_complex, lwork);
      zgetrf_(&n, &m, (void*)COMPLEXP(a), &lda, ipiv, &info);
      zgetri_(&n, (void*)COMPLEXP(a), &lda, ipiv, work, &lwork, &info);
      free(work);
    }
  free(ipiv);
  return Py_BuildValue("i", info);
}

PyObject* inverse_symmetric(PyObject *self, PyObject *args)
{
  PyArrayObject* a;
  if (!PyArg_ParseTuple(args, "O", &a))
    return NULL;
  int n = PyArray_DIMS(a)[0];
  int lda = n;
  int lwork =n;
  int* ipiv = GPAW_MALLOC(int, n);
  int info = 0;
  if (PyArray_DESCR(a)->type_num == NPY_DOUBLE)
    {
      double* work = GPAW_MALLOC(double, lwork);
      dsytrf_("U", &n, DOUBLEP(a), &lda, ipiv, work, &lwork, &info);
      dsytri_("U", &n, DOUBLEP(a), &lda, ipiv, work, &info);
      free(work);
    }
  else
    {
      void *work = GPAW_MALLOC(double_complex, lwork);
      zsytrf_("U", &n, (void*)COMPLEXP(a), &lda, ipiv, work, &lwork, &info);
      zsytri_("U", &n, (void*)COMPLEXP(a), &lda, ipiv, work, &info);
      free(work);
    }
  free(ipiv);
  return Py_BuildValue("i", info);
}

PyObject* linear_solve_band(PyObject *self, PyObject *args)
{
 PyArrayObject* a;
 PyArrayObject* b;
 int  kl, ku, info=0, *ipiv;
 if(!PyArg_ParseTuple(args,"OOii",&a, &b,&kl,&ku))
   return NULL;
int n=PyArray_DIMS(a)[0];
int ldab=PyArray_DIMS(a)[1];
int ldb=PyArray_DIMS(b)[0];
int nrhs=PyArray_DIMS(b)[1];
   ipiv = GPAW_MALLOC(int, n);
   zgbsv_(&n, &kl,&ku, &nrhs, (void*)COMPLEXP(a), &ldab, ipiv, (void*)COMPLEXP(b), &ldb, &info);
   free(ipiv);
 return Py_BuildValue("i",info);
}

PyObject* linear_solve_tridiag(PyObject *self, PyObject *args)
{
 PyArrayObject* A;
 PyArrayObject* du;
 PyArrayObject* du2;
 PyArrayObject* dl;
 PyArrayObject* phi;
 int dim=0, one=1, info=0;
 if(!PyArg_ParseTuple(args,"iOOOOO", &dim, &A, &du, &dl, &du2, &phi))
   return NULL;
 int ldb = dim;
 int *ipiv = GPAW_MALLOC(int, dim);
 zgttrf_(&dim, (void*)COMPLEXP(dl), (void*)COMPLEXP(A), (void*)COMPLEXP(du), (void*)COMPLEXP(du2), ipiv, &info);
 zgttrs_("N", &dim, &one, (void*)COMPLEXP(dl), (void*)COMPLEXP(A), (void*)COMPLEXP(du),
                                   (void*)COMPLEXP(du2), ipiv, (void*)COMPLEXP(phi), &ldb, &info);
 free(ipiv);
 return Py_BuildValue("i",info);
}
