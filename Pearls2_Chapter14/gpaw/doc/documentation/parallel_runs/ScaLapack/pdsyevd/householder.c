#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#ifdef BGP
#include <mpix.h>
#endif

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
#define   Csys2blacs_handle_ Csys2blacs_handle
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

int Csys2blacs_handle_(MPI_Comm SysCtxt);

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
#define   pdlaset_ pdlaset
#define   pzlaset_ pzlaset
#define   pdgemr2d_  pdgemr2d
#define   pdlamch_  pdlamch

#define   pdpotrf_  pdpotrf
#define   pdpotri_  pdpotri
#define   pzpotrf_  pzpotrf
#define   pzpotri_  pzpotri
#define   pdtrtri_  pdtrtri
#define   pztrtri_  pztrtri

#define   pdsytrd_  pdsytrd

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

void pzelset_(void* a,int* ia,int* ja,int* desca,void* alpha);

void pdlaset_(char* uplo, int* m, int* n, double* alpha, double* beta,
	      double* a, int* ia, int* ja, int* desca);

void pzlaset_(char* uplo, int* m, int* n, void* alpha, void* beta,
	      void* a, int* ia, int* ja, int* desca);

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

void pdsytrd_(char *uplo, int *n, double* a, int *ia, int* ja,
              int* desca, double *d, double *e, double *tau,
              double *work, int *lwork, int *info);

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

// int MPIX_Get_last_algorithm(MPI_Comm comm, int *lastalgorithm);

int main(int argc, char *argv[]) {

     // Some constants
     int minusone = -1;
     int zero = 0;
     int one = 1;
     double dzero = 0.0; 

     // ConText
     int ConTxt = minusone;

     // order
     char order = 'R';
     char scope = 'A';

     // root process
     int root = zero;

     // BLACS/SCALAPACK parameters
     // the size of the blocks the distributed matrix is split into
     // (applies to both rows and columns)
     int mb = 32;
     int nb = mb; // PDSYEVxxx constraint

     // the number of rows and columns in the processor grid
     // only square processor grids due to C vs. Fortran ordering
     int nprow = 2;
     int npcol = nprow; // only square processor grids, 

     // starting row and column in grid, do not change
     int rsrc = zero; 
     int csrc = zero;

     // dimensions of the matrix to diagonalize
     int m = 1000;
     int n = m; // only square matrices

     int info = zero;

     // Rest of code will only work for:
     // nprow = npcol
     // mb = nb;
     // m = n;
     // rsrc = crsc;

     // Paramteres for Trivial Matrix
     double alpha = 0.1; // off-diagonal
     double beta = 75.0; // diagonal
     
     // For timing:
     double tdiag0, tdiag, ttotal0, ttotal;

     // BLACS Communicator
     MPI_Comm blacs_comm;
     int nprocs;
     int iam;
     int myrow, mycol;
     int result;

     MPI_Init(&argc, &argv);
     MPI_Barrier(MPI_COMM_WORLD);
     ttotal0 = MPI_Wtime();
     MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
     MPI_Comm_rank(MPI_COMM_WORLD, &iam);

     if (argc > one) {
       nprow = strtod(argv[1],NULL);
       m = strtod(argv[2],NULL);
       npcol = nprow;
       n = m;
     }
    
     // We can do this on any subcommunicator.
#ifdef CartComm
     int dim[2];
     int pbc[2];
     dim[0] = nprow;
     dim[1] = npcol;
     pbc[0] = 0;
     pbc[1] = 0;
     MPI_Cart_create(MPI_COMM_WORLD, 2, dim, pbc, 1, &blacs_comm);
#else
     blacs_comm = MPI_COMM_WORLD;
#endif


     if (iam == root) {
       printf("world size %d \n",nprocs);
       printf("n %d \n", n);
       printf("nprow %d \n", nprow);
       printf("npcol %d \n", npcol);
#ifdef BGP
       MPIX_Get_property(blacs_comm, MPIDO_RECT_COMM, &result);
       if (result) printf("this is a rectangular communicator\n");
#endif
     }

     // initialize the grid
     // The lines below are equivalent to the one call to:
     if (blacs_comm != MPI_COMM_NULL) {
       ConTxt = Csys2blacs_handle_(blacs_comm);
       Cblacs_gridinit_(&ConTxt, &order, nprow, npcol);

       // get information back about the grid
       Cblacs_gridinfo_(ConTxt, &nprow, &npcol, &myrow, &mycol);
     }

     if (ConTxt != minusone) {

          int desc[9];

          // get the size of the distributed matrix
          int locM = numroc_(&m, &mb, &myrow, &rsrc, &nprow);
          int locN = numroc_(&n, &nb, &mycol, &csrc, &npcol);

	  // printf ("locM = %d \n", locM);
	  // printf ("locN = %d \n", locN);

          int lld = MAX(one,locM);

          // build the descriptor
          descinit_(desc, &m, &n, &mb, &nb, &rsrc, &csrc, &ConTxt, &lld, &info);
          // Allocate arrays
	  double* mata = malloc(locM*locN * sizeof(double));
	  double* matd = malloc(locM*locN * sizeof(double));
	  double* mate = malloc(locM*locN * sizeof(double));
	  double* tau  = malloc(locM*locN * sizeof(double));

	  double* work;
          work = malloc(3 * sizeof(double));
          int querylwork = minusone;
          
          // Build a trivial distributed matrix: Diagonal matrix
	  char uplo = 'L'; // work with upper
	  pdlaset_(&uplo, &m, &n, &alpha, &beta, mata, &one, &one, desc);

	  // First there is a workspace query
          pdsytrd_(&uplo, &n, mata, &one, &one, desc, matd, mate, tau,
		   work, &querylwork, &info);


          int lwork = (int)work[0];
          //printf("lwork %d\n", lwork);
          free(work);

          work = (double*)malloc(lwork * sizeof(double));
	  // This is actually performs the tridiagonalization
	  tdiag0 = MPI_Wtime();
          Cblacs_barrier(ConTxt, &scope);
          pdsytrd_(&uplo, &n, mata, &one, &one, desc, matd, mate, tau,
		   work, &lwork, &info);
          Cblacs_barrier(ConTxt, &scope);
          tdiag = MPI_Wtime() - tdiag0;

          free(work);
          free(mata);
	  free(matd);
	  free(mate);
	  free(tau);

          // Destroy BLACS grid
          Cblacs_gridexit_(ConTxt);

	  // Print time
	  if (myrow == zero && mycol == zero) {
	    if (info != zero) {
	      printf("info = %d \n", info);
	    }
	    
	    printf("Time (s) diag: %f\n", tdiag);
	  }

     }

     MPI_Barrier(MPI_COMM_WORLD);
     ttotal = MPI_Wtime() - ttotal0;
     if (iam == 0)
          printf("Time (s) total: %f\n", ttotal);
     MPI_Finalize();
}
