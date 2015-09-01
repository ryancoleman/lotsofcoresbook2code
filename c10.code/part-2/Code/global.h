#ifndef H_GLOBAL
#define H_GLOBAL

// Declare user defined structures

#define MAXLEN 200   // max length of char arrays
#define MAXLINES 10000000 // max lines in input files

// JB 
#ifdef __INTEL_OFFLOAD
#define _OFFLOADABLE __declspec(target(mic))
#else
#define _OFFLOADABLE 
#endif

#define ALLOC alloc_if(1)
#define FREE free_if(1)
#define RETAIN free_if(0)
#define REUSE alloc_if(0)

//#include <mkl.h>
#include "local_barrier.h"

// useful macros
#define min(x,y) (((x) < (y)) ? (x) : (y))

// SJP: Funny that this isn't provided.
#include <stdlib.h>
_OFFLOADABLE
inline omp_get_max_teams()
{
    return 59;
}

/* JB struct for thread private data in the gamma_3d algorithm */
typedef struct
{
    int i0;
    int j0;
    int k0;
    int i1;
    int j1;
    int k1;
    int xsize;
    double* mvec;
    double* xvec;
    double *yvec[4];    // a hack to enable re-use between offload and offload_team implementation.
    double* xdiff;
    double* ixdiff;
    double* intgrlvec;  // array of x in calculate_gamma_3d 
    double* plijkz;     // array of plijk*z in calculate_gamma_3d
    barrier_t* bar;     // ptr to struct for core barrier data
} thread_g3d_data_t;

// Declare external variables

extern char alpha_data_file[MAXLEN];
extern char alpha_size_file[MAXLEN];
extern char bessel_data_file[MAXLEN];
extern char bessel_size_file[MAXLEN];
extern char transfer_data_file[MAXLEN];
extern char transfer_size_file[MAXLEN];
extern char source_data_file[MAXLEN];
extern char source_size_file[MAXLEN];
extern char limber_data_file[MAXLEN];
extern char limber_size_file[MAXLEN];
extern char cls_data_file[MAXLEN];
extern char decompose_data_file[MAXLEN];
extern char decompose_size_file[MAXLEN];
extern char eigen_data_file[MAXLEN];
extern char eigen_tri_data_file[MAXLEN];
extern char ortho_data_file[MAXLEN];
extern char orthol_data_file[MAXLEN];
extern char ortho_tri_data_file[MAXLEN];
extern char orthol_tri_data_file[MAXLEN];
extern char lambda_data_file[MAXLEN];
extern char lambdal_data_file[MAXLEN];
extern char lambda_tri_data_file[MAXLEN];
extern char lambdal_tri_data_file[MAXLEN];
extern char gamma_data_file[MAXLEN];
extern char gamma_tri_data_file[MAXLEN];
extern char proj_data_file[MAXLEN];
extern char proj_size_file[MAXLEN];
extern char proj_tri_data_file[MAXLEN];
extern char proj_tri_size_file[MAXLEN];
extern char modes_data_file[MAXLEN];
extern char modes_tri_data_file[MAXLEN];
extern char flat_data_file[MAXLEN];
extern char flat_size_file[MAXLEN];
extern char l_data_file[MAXLEN];
extern char l_size_file[MAXLEN];
extern char bispectrum_file[MAXLEN];
extern char trispectrum_file[MAXLEN];
extern char interpolated_l_data_file[MAXLEN];
extern char interpolated_l_size_file[MAXLEN];
extern char interpolated_bispectrum_file[MAXLEN];
extern char beam_noise_file[MAXLEN];
extern char transfer_wgt_file[MAXLEN];
extern char lens_data_file[MAXLEN];
extern char restart_file[MAXLEN];
extern char bload_file[MAXLEN];
extern char edist_file[MAXLEN];
extern char edist_tri_file[MAXLEN];
extern char mdist_file[MAXLEN];
extern char mdist_tri_file[MAXLEN];

extern double deltaphi;
extern int alpha_max;

// Flags
//extern int use_l_file;
extern char lini_file[MAXLEN];
extern int eflag_order_prim;
extern int eflag_order_late;
extern int eflag_lmax;
extern int clflag_uselens;

// named

extern int initilise(char *inifile);
extern int basis_functions_legendre(double **array, int size, int pmax, double max, double *values);
extern int basis_functions_fourier(double **array, int size, int pmax, double max, double *values);
extern int basis_functions_sinlog(double *array_flat, int size, int pmax, double min, double max, double *values);
extern int basis_functions_bi(double *array_flat, int size, int pmax, double max, double *values);
extern _OFFLOADABLE double calculate_geometric(int l1, int l2, int l3);
extern _OFFLOADABLE double calculate_xint(int n, int xsize, double* x, double* y, int l1, int l2, int l3, double* xdiff, double* ixdiff);
extern void precompute_gamma_3d();
extern _OFFLOADABLE void decompose_gamma_3d_mpi(thread_g3d_data_t* thread_data, int rank, int numranks, double workshare);
_OFFLOADABLE void calculate_gamma_3d_nested(thread_g3d_data_t* thread_data, int maxblocksize);
_OFFLOADABLE void gamma_3d_offload(double * restrict gamma_flat, int gamma_flat_size, int rank, int numranks, double workshare);
_OFFLOADABLE void gamma_3d_offload_teams(double * restrict gamma_flat, int gamma_flat_size, int rank, int numranks, double workshare);
_OFFLOADABLE void gamma_3d_offload_nested(double * restrict gamma_flat, int gamma_flat_size, int rank, int numranks, double workshare);
_OFFLOADABLE void gamma_3d_host(double * restrict gamma_flat, int gamma_flat_size, int rank, int numranks, double workshare);

extern void print_affinity(int team_id, int thread_id);

// fixed data

extern int init_lmax(int l);
extern _OFFLOADABLE int get_lmax();

extern int set_terms_prim();
extern _OFFLOADABLE int get_terms_prim();
extern _OFFLOADABLE int get_pmax_prim();
extern int create_order_prim();
extern _OFFLOADABLE void find_perm_prim(int r,int* p1,int* p2,int* p3);
extern int set_terms_late();
extern int get_terms_late();
extern _OFFLOADABLE int get_pmax_late();
extern int create_order_late();
extern _OFFLOADABLE void find_perm_late(int r,int* p1,int* p2,int* p3);

extern int create_cl(int l_size);
extern int load_cl(int l_size ,int* l_values);
extern double* get_cl_array();
extern double get_cl(int i);

extern int load_BN(int l_size ,int* l_values);
extern int load_TL(int l_size ,int* l_values);
extern int load_PBN(int l_size ,int* l_values);
extern int load_lens(int l_size ,int* l_values);
extern int create_beam(int l_size);
extern double* get_beam_array();
extern double get_beam(int i);
extern int create_noise(int l_size);
extern double* get_noise_array();
extern double get_noise(int i);
extern int create_t_wgt(int l_size);
extern int create_lens(int l_size);
extern double get_lens_tt(int i);
extern double get_lens_tp(int i);

extern int create_basis_late(int ksize, double kmax, double *k);
extern _OFFLOADABLE double get_basis_late(int i, int j);
extern _OFFLOADABLE double* get_basis_late_array();

extern int create_beta(int l_size, int x_size);
extern _OFFLOADABLE double* get_beta_array();
extern int read_beta();

extern _OFFLOADABLE int get_b_lsize();
extern _OFFLOADABLE int get_b_xsize();
extern int get_b_lvec(int *vec);
extern _OFFLOADABLE int get_b_xvec(double *vec);

extern double get_orthol(int i,int j);
extern int read_orthol();

extern double get_lambdal(int i,int j);
extern int read_lambdal();

extern int create_gamma();
extern int update_gamma(double *results);
extern int output_gamma();
extern _OFFLOADABLE double permsix(int i, int j, int k);

 // arrays
extern void ivector_read(int *n, char filename[MAXLEN], int grid[(*n)]);
extern int *create_ivector(int length);
extern double *create_vector(int length);
extern double **create_array(int dim_x, int dim_y);
extern void destroy_array(double **array);
extern void array_write(int *n, char filename[MAXLEN], double grid[(*n)]);
extern void array_read(int *n, char filename[MAXLEN], double grid[(*n)]);

// timing
extern double csecond();

// io
extern int load_txt_dbl(char* filename, int columns, double* values, int* size);
extern int load_two(char* filename, double* values, int* size);
extern int load_three(char* filename, double* values, int* size);
extern int load_three_int(char* filename, int* values, int* size);
extern double get_tau0();

#endif /* H_GLOBAL INCLUDED */
