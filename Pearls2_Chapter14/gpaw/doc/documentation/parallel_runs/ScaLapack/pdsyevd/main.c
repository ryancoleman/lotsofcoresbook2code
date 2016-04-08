#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "mpi.h"

#define MIN(x, y) ((x) < (y) ? (x) : (y))
#define MAX(x, y) ((x) > (y) ? (x) : (y))

// BLACS
#ifdef GPAW_NO_UNDERSCORE_CBLACS
#define   Cblacs_barrier_  Cblacs_barrier
#define   Cblacs_exit_     Cblacs_exit
#define   Cblacs_get_      Cblacs_get
#define   Cblacs_gridexit_ Cblacs_gridexit
#define   Cblacs_gridinfo_ Cblacs_gridinfo
#define   Cblacs_gridinit_ Cblacs_gridinit
#define   Cblacs_pinfo_    Cblacs_pinfo
#define   Cblacs_pnum_     Cblacs_pnum
#define   Cblacs_setup_    Cblacs_setup
#endif

#ifdef GPAW_NO_UNDERSCORE_BLACS
#define   dgebr2d_  dgebr2d
#define   dgebs2d_  dgebs2d
#define   zgebr2d_  zgebr2d
#define   zgebs2d_  zgebs2d
#endif

void Cblacs_barrier_(int ConTxt, char *scope);

void Cblacs_exit_(int NotDone);

void Cblacs_get_(int ConTxt, int what, int* val);

void Cblacs_gridexit_(int ConTxt);

void Cblacs_gridinfo_(int ConTxt, int *nprow, int *npcol,
                      int *myrow, int *mycol);

void Cblacs_gridinit_(int *ConTxt, char* order, int nprow, int npcol);

void Cblacs_pinfo_(int *mypnum, int *nprocs);

int Cblacs_pnum_(int ConTxt, int prow, int pcol);

void Cblacs_setup_(int *mypnum, int *nprocs);

void dgebr2d_(int *ConTxt, char* scope, char* top, int *m, int *n,
              double *A, int *lda, int *rsrc, int *csrc);

void dgebs2d_(int *ConTxt, char* scope, char* top, int *m, int *n,
              double *A, int *lda);

void zgebr2d_(int *ConTxt, char* scope, char* top, int *m, int *n,
              double *A, int *lda, int *rsrc, int *csrc);

void zgebs2d_(int *ConTxt, char* scope, char* top, int *m, int *n,
              double *A, int *lda);
// End of BLACS

// ScaLapack

#ifdef GPAW_NO_UNDERSCORE_SCALAPACK
#define   descinit_  descinit
#define   numroc_  numroc
#define   pdelset_  pdelset
#define   pzelset_  pzelset
#define   pdgemr2d_  pdgemr2d
#define   pdlamch_  pdlamch

#define   pdpotrf_  pdpotrf
#define   pdpotri_  pdpotri
#define   pzpotrf_  pzpotrf
#define   pzpotri_  pzpotri
#define   pdtrtri_  pdtrtri
#define   pztrtri_  pztrtri

#define   pdsyevd_  pdsyevd
#define   pdsyev_  pdsyev
#define   pdsyevx_  pdsyevx
#define   pdsygvx_  pdsygvx
#define   pzheev_  pzheev
#define   sl_init_  sl_init
#endif

#ifdef GPAW_NO_UNDERSCORE_CSCALAPACK
#define   Cpdgemr2d_  Cpdgemr2d
#endif

#ifdef GPAW_SECOND_UNDERSCORE_SL_INIT
#define   sl_init_  sl_init__
#endif

void descinit_(int* desc, int* m, int* n, int* mb, int* nb, int* irsrc,
               int* icsrc, int* ictxt,
               int* lld, int* info);

int numroc_(int* n, int* nb, int* iproc, int* isrcproc, int* nprocs);

void pdelset_(double* a,int* ia,int* ja,int* desca,double* alpha);

//void pdelset2_(double* alpha,double* a,int* ia,int* ja,int* desca,double* beta);

void pzelset_(void* a,int* ia,int* ja,int* desca,void* alpha);

int pdgemr2d_(int* m, int*n, double* a, int* ia, int* ja, int* desca,
              double* b, int* ib, int* jb, int* descb, int* ctxt);

double pdlamch_(int* ictxt, char* cmach);

void Cpdgemr2d_(int m, int n, double *A, int IA, int JA, int *descA,
                double *B, int IB, int JB, int *descB, int gcontext);

// cholesky
void pdpotrf_(char *uplo, int *n, double* a, int *ia, int* ja, int* desca,
              int *info);
void pdpotri_(char *uplo, int *n, double* a, int *ia, int* ja, int* desca,
              int *info);

void pzpotrf_(char *uplo, int *n, void* a, int *ia, int* ja, int* desca,
              int *info);
void pzpotri_(char *uplo, int *n, void* a, int *ia, int* ja, int* desca,
              int *info);


void pdtrtri_(char *uplo, char *diag, int *n, double* a, int *ia, int* ja,
              int* desca, int *info);
void pztrtri_(char *uplo, char *diag, int *n, void* a, int *ia, int* ja,
              int* desca, int *info);


void pdsyevd_(char *jobz, char *uplo, int *n, double* a, int *ia, int* ja,
              int* desca, double *w,
              double* z, int *iz, int* jz, int* descz, double *work,
              int *lwork, int *iwork, int *liwork, int *info);

void pdsyev_(char *jobz, char *uplo, int *n, double* a, int *ia, int* ja,
             int* desca, double *w,
             double* z, int *iz, int* jz, int* descz, double *work,
             int *lwork, int *info);

void pdsyevx_(char *jobz, char *range, char *uplo, int *n, double* a,
              int *ia, int* ja, int* desca, double* vl,
              double* vu, int* il, int* iu, double* abstol, int* m, int* nz,
              double* w, double* orfac, double* z, int *iz,
              int* jz, int* descz, double *work, int *lwork, int *iwork,
              int *liwork, int *ifail, int *iclustr, double* gap, int *info);

void pdsygvx_(int *ibtype, char *jobz, char *range, char *uplo, int *n,
              double* a, int *ia, int* ja,
              int* desca, double* b, int *ib, int* jb, int* descb,
              double* vl, double* vu, int* il, int* iu,
              double* abstol, int* m, int* nz, double* w, double* orfac,
              double* z, int *iz, int* jz, int* descz,
              double *work, int *lwork, int *iwork, int *liwork, int *ifail,
              int *iclustr, double* gap, int *info);

void pzheev_(char *jobz, char *uplo, int *n, double* a, int *ia, int* ja,
             int* desca, double *w, double* z, int *iz, int* jz,
             int* descz, double *work, int *lwork, double *rwork,
             int *lrwork, int *info);

void sl_init_(int* ictxt, int* nprow, int* npcol);
// End of ScaLapack

void buildCustomMatrix(double* m, int row, int col)
{
     int r,c;

     for(r = 0; r < row; r++)
     {
          for(c = 0 ; c < col; c++) {
               m[r*row+c] = -(r+c+1);
               if (r != c) {m[r*row+c] = m[r*row+c] / (row+col+1);};
               //// if (r != c) {m[r*row+c] = 0.0;};
          }
     }
}

void printMatrix(double* m, int row, int col, int myrow, int mycol)
{
     int r,c;

     printf("p(%d, %d)\n", myrow, mycol);
     for(r = 0; r < row; r++)
     {
          for(c = 0 ; c < col; c++)
               printf("%f ", m[r*row+c]);
          printf("\n");
     }
}

int main(int argc, char **argv) {

     static int minusone = -1;
     static int zero = 0;
     static int one = 1;

     // root process
     int root = zero;

     // the size of the blocks the distributed matrix is split into
     // (applies to both rows and columns)
     int mb = 1;
     int nb = mb;
     // the number of rows in the process grid
     // over which the matrix is distributed
     int nprow = 1;
     // the number of columns in the process grid
     // over which the matrix is distributed
     int npcol = 1;
     // dimensions of the matrix to diagonalize
     int n = 2;
     int m = n;
     int lda = n;
     int ldb = n;

     int ibtype = 1;
     int info = 0;
     int iam = -1;

     int ConTxt;

     // 2x2 blacs grid
     nprow = 1;
     npcol = 1;

     MPI_Init(&argc, &argv);

     // initialize the grid
     sl_init_(&ConTxt, &nprow, &npcol);

     // get information back about the grid
     int myrow = -1;
     int mycol = -1;
     Cblacs_gridinfo_(ConTxt, &nprow, &npcol, &myrow, &mycol);

     int pnum;

     char TOP = ' '; // ?
     char scope = 'A'; // All grid

     int rsrc = 0;
     int csrc = 0;

     // matrix to be diagonalized
     double* a = (double*)malloc(n*m * sizeof(double));

     // eigenvalues
     double* eigvals = (double*)malloc(n * sizeof(double));

     int row = n;
     int col = m;

     if (myrow == 0 && mycol == 0) {
          buildCustomMatrix(a, row, col);
     }

     if (myrow == 0 && mycol == 0) {
          printMatrix(a, row, col, myrow, mycol);
     }

     if (myrow != -1 && mycol != -1) {

          // Returns the system process number of the process in the process grid
          pnum = Cblacs_pnum_(ConTxt, myrow, mycol);

          // build the descriptor
          int desc0[9];
          descinit_(desc0, &m, &n, &m, &n, &rsrc, &csrc, &ConTxt, &m, &info);

          // distribute the full matrix to the process grid
          if (pnum == root)
          {
               //printf("C diagonalize ScaLapack\n");
               // http://icl.cs.utk.edu/lapack-forum/viewtopic.php?p=87&
               // distribute the matrix
               // Uncomment for PDSYEV:
               //printf("run %d\n", 1);
               dgebs2d_(&ConTxt,&scope,&TOP,&m,&n,a,&lda);
               //printf("run %d\n", 2);
               // Uncomment for PZHEEV:
               //zgebs2d_(&ConTxt,&scope,&TOP,&m,&n,&a,&lda);
          }
          else
          {
               // receive the matrix
               // Uncomment for PDSYEV:
               //printf("run %d\n", 3);
               dgebr2d_(&ConTxt,&scope,&TOP,&m,&n,a,&lda,&rsrc,&csrc);
               //printf("run %d\n", 4);
               // Uncomment for PZHEEV:
               //zgebr2d_(&ConTxt,&scope,&TOP,&m,&n,&a,&lda,&rsrc,&csrc);
          }

          int desc[9];

          // get the size of the distributed matrix
          int locM = numroc_(&m, &mb, &myrow, &rsrc, &nprow);
          int locN = numroc_(&n, &nb, &mycol, &csrc, &npcol);
          // allocate the distributed matrices
          double* mata = (double*)malloc(locM*locN * sizeof(double));
          // allocate the distributed matrix of eigenvectors
          double* z = (double*)malloc(locM*locN * sizeof(double));

          int lld = MAX(1,locM);

          // build the descriptor
          descinit_(desc, &m, &n, &mb, &nb, &rsrc, &csrc, &ConTxt, &lld, &info);
          // build the distributed matrices
          for(int i1 = one; i1 < m+one; ++i1)
               for(int i2 = one; i2 < n+one; ++i2)
               {
                    // Uncomment for PDSYEV:
                    // http://icl.cs.utk.edu/lapack-forum/viewtopic.php?t=321&sid=f2ac0a3c06c66d74a2fbd65c222ffdb0
                    pdelset_(mata,&i1,&i2,desc,&a[(i1-one)*n+i2-one]);
                    // Uncomment for PZHEEV:
                    //pzelset_(mata,&i1,&i2,desc,&a(i1-one,i2-one));
               }

          char jobz = 'V'; // eigenvectors also
          char range = 'A'; // all eiganvalues
          char uplo = 'U'; // work with upper

          double vl, vu;
          int il, iu;

          char cmach = 'U';

          double abstol = 2.0 * pdlamch_(&ConTxt, &cmach);

          int eigvalm, nz;

          double orfac = -1.0;
          //double orfac = 0.001;

          int* ifail;
          ifail = (int*)malloc(m * sizeof(int));

          int* iclustr;
          iclustr =  (int*)malloc(2*nprow*npcol * sizeof(int));

          double* gap;
          gap =  (double*)malloc(nprow*npcol * sizeof(double));

          double* work;
          work = (double*)malloc(3 * sizeof(double));
          int querylwork = -1;
          int* iwork;
          iwork = (int*)malloc(1 * sizeof(int));
          int queryliwork = 1;

          //pdsyevx_(&jobz, &range, &uplo, &n, mata, &one, &one, desc, &vl,
          //         &vu, &il, &iu, &abstol, &eigvalm, &nz, eigvals, &orfac, z, &one,
          //         &one, desc, work, &querylwork, iwork, &queryliwork, ifail, iclustr, gap, &info);
          pdsyevd_(&jobz, &uplo, &n, mata, &one, &one, desc, eigvals,
                   z, &one, &one, desc,
                   work, &querylwork, iwork, &queryliwork, &info);
          //pdsyev_(&jobz, &uplo, &m, mata, &one, &one, desc, eigvals,
          //        z, &one, &one, desc, work, &querylwork, &info);

          int lwork = (int)work[0];
          //printf("lwork %d\n", lwork);
          free(work);
          int liwork = (int)iwork[0];
          //printf("liwork %d\n", liwork);
          free(iwork);

          work = (double*)malloc(lwork * sizeof(double));
          iwork = (int*)malloc(liwork * sizeof(int));

          //pdsyevx_(&jobz, &range, &uplo, &n, mata, &one, &one, desc, &vl,
          //         &vu, &il, &iu, &abstol, &eigvalm, &nz, eigvals, &orfac, z, &one,
          //         &one, desc, work, &lwork, iwork, &liwork, ifail, iclustr, gap, &info);
          pdsyevd_(&jobz, &uplo, &n, mata, &one, &one, desc, eigvals,
                   z, &one, &one, desc,
                   work, &lwork, iwork, &liwork, &info);
          //pdsyev_(&jobz, &uplo, &m, mata, &one, &one, desc, eigvals,
          //        z, &one, &one, desc, work, &lwork, &info);

          free(work);
          free(iwork);

          Cpdgemr2d_(m, n, z, one, one, desc, a, one, one, desc0, ConTxt);
          // a constains eigenvectors now

          free(gap);
          free(iclustr);
          free(ifail);

          free(z);
          free(mata);
          // clean up the grid
          Cblacs_gridexit_(ConTxt);
     }
     free(a);
     free(eigvals);
     Cblacs_exit_(zero);
     //MPI_Finalize();
     exit(0);
}
