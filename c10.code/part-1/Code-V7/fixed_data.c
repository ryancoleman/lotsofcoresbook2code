#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_spline.h>
#include "global.h"

_OFFLOADABLE static double *basis_late_flat;
_OFFLOADABLE static double *beta_flat;
static double **orthol;
static double **lambdal;
static double **gamma_;
static double *cl;
static double *lens_tt;
static double *lens_tp;
static double *beam;
static double *noise;
static double *t_wgt;

_OFFLOADABLE static int lmax;
_OFFLOADABLE static int pmax_prim;
_OFFLOADABLE static int pmax_late;
_OFFLOADABLE static int terms_prim;
static int terms_late;
_OFFLOADABLE static int *order_prim_flat;
_OFFLOADABLE static int *order_late_flat;
_OFFLOADABLE static int b_lsize;
_OFFLOADABLE static int b_xsize;
static int *b_lvec;
_OFFLOADABLE static double *b_xvec;

static long get_term_flat_value(int *term, int pmax);
static void sort_terms_array(int *flat_array, int lo, int hi, int pmax);

int init_lmax(int l)
{
	lmax = l;

	// JB offload lmax value
	int offload_target;
    MPI_Comm_rank(MPI_COMM_WORLD, &offload_target);
	#pragma offload_transfer target(mic:offload_target) in(lmax)
	return 0;
}

_OFFLOADABLE int get_lmax()
{
	return lmax;
}

/**
 * Jb This takes a {a1,a2,a3} triple and flattens it into one value
 * for the comparison operator in the sort terms function.
 * pmax is the largest size the a1,a2,a3 can take
 */
static long get_term_flat_value(int *term, int pmax)
{
    long pmax1 = (long) pmax+1;
    long val = pmax1 * (pmax1 * term[0] + term[1]) + term[2];
    return val;
}

/*
 * JB given an array of {a1,a2,a3} triples, sort the array in order of
 * acscending flat value of the triples. Implemented with quicksort.
 */
static void sort_terms_array(int *flat_array, int lo, int hi, int pmax)
{
    // alias flat array as a n*3 dimensional array
    int (* restrict array)[3] = (int (*restrict)[3]) flat_array;
    int pivot, j, i;

    if (lo < hi)
    {
        pivot = lo;
        i = lo;
        j = hi;

        // partition
        while (i<j)
        {
            while (get_term_flat_value(&(array[i][0]),pmax) \
                    < get_term_flat_value(&(array[pivot][0]), pmax) )
            {
                i++;
                if (i==hi) break;
            }
            while (get_term_flat_value(&(array[j][0]), pmax) \
                    > get_term_flat_value(&(array[pivot][0]), pmax) )
            {
                j--;
                if (j==lo) break;
            }
            if (i >= j) break;

            // swap 
            int temp[] = {array[i][0], array[i][1], array[i][2]};
            array[i][0] = array[j][0];
            array[i][1] = array[j][1];
            array[i][2] = array[j][2];
            array[j][0] = temp[0];
            array[j][1] = temp[1];
            array[j][2] = temp[2];
        }
        // put pivot item at a[j]
        int temp[] = {array[pivot][0],array[pivot][1],array[pivot][2]};
        array[pivot][0] = array[j][0];
        array[pivot][1] = array[j][1];
        array[pivot][2] = array[j][2];
        array[j][0] = temp[0];
        array[j][1] = temp[1];
        array[j][2] = temp[2];

        sort_terms_array(flat_array, lo, j-1, pmax);
        sort_terms_array(flat_array, j+1, hi, pmax);
    }
    return;
}

int set_terms_prim(){
	
	int r,i,j,k,n,s;
	int* order_raw = malloc(sizeof(int)*3*MAXLINES);
	int* size = malloc(sizeof(int));
	
	terms_prim = alpha_max;
	pmax_prim = 0;
	create_order_prim();
	int (* restrict order_prim)[3] = (int (*restrict)[3]) order_prim_flat;
			
	switch(eflag_order_prim){
		case 1:
			i=j=k=s=0;
			for (n=0;n<terms_prim;n++){
				order_prim[n][0] = i;
				order_prim[n][1] = j;
				order_prim[n][2] = k;
				if(k-j<2&&j-i<2&&k-i<2){
					s++;
					i=j=0;
					k=s;
				} else if(k-j>1){
					k--;
					j++;
				} else {
					i++;
					j=i;
					k=s-2*i;
				}
				if (order_prim[n][2]>pmax_prim)pmax_prim = order_prim[n][2];
			}
		break;
		case 2:
			i=j=k=s=0;
			for (n=0;n<terms_prim;n++){
				order_prim[n][0] = i;
				order_prim[n][1] = j;
				order_prim[n][2] = k;
				if(i<j){
					j--;
					k++;
				} else if(i>0){
					i--;
					j = (s-i)/2;
					k = s-j-i;
				} else {
					s++;
					i = s/3;
					j = (s-i)/2;
					k = s-i-j;
				}
				if (order_prim[n][2]>pmax_prim)pmax_prim = order_prim[n][2];
			}
		break;
		case 3:
			load_three_int(edist_file, order_raw, size);
			if(terms_prim>*size){
				printf("eigenfunction %d exceeds that in distance file\n",terms_prim);
				exit;
			} else {
				for (n=0;n<terms_prim;n++){
					order_prim[n][0] = order_raw[3*n];
					order_prim[n][1] = order_raw[3*n+1];
					order_prim[n][2] = order_raw[3*n+2];
					if (order_prim[n][2]>pmax_prim)pmax_prim = order_prim[n][2];
				}
			}
		break;
		case 4:
			terms_prim = alpha_max;
			pmax_prim = 0;
			//create_order_prim(); JB bug
			load_three_int(edist_file, order_raw, size);
			if(terms_prim>*size){
				printf("eigenfunction %d exceeds that in distance file\n",terms_prim);
				exit;
			} else {
				for (n=0;n<terms_prim-1;n++){
					order_prim[n][0] = order_raw[3*n];
					order_prim[n][1] = order_raw[3*n+1];
					order_prim[n][2] = order_raw[3*n+2];
					if (order_prim[n][2]>pmax_prim)pmax_prim = order_prim[n][2];
				}
			}
			
			order_prim[0][0] = 0;
			order_prim[0][1] = 0;
			order_prim[0][2] = 0;
			
			for (n=terms_prim-1;n>1;n--){
				order_prim[n][0] = order_prim[n-1][0];
				order_prim[n][1] = order_prim[n-1][1];
				order_prim[n][2] = order_prim[n-1][2];
			}
			
			order_prim[1][0] = pmax_prim+1;
			order_prim[1][1] = pmax_prim+1;
			order_prim[1][2] = pmax_prim+2;
		break;
		case 5:
			load_three_int(edist_file, order_raw, size);
			if(terms_prim>*size){
				printf("eigenfunction %d exceeds that in distance file\n",terms_prim);
				exit;
			} else {
				for (n=0;n<terms_prim;n++){
					order_prim[n][0] = order_raw[3*n];
					order_prim[n][1] = order_raw[3*n+1];
					order_prim[n][2] = order_raw[3*n+2];
					if (order_prim[n][2]>pmax_prim)pmax_prim = order_prim[n][2];
				}
			}
		break;
		default:
			printf("invalid mode specified for eigen function ordering\n");
			exit;
		break;
	}

	pmax_prim = 0;
	for (n=0;n<terms_prim;n++){
		if(order_prim[n][2]>pmax_prim) pmax_prim = order_prim[n][2];
	}

    #ifdef SORT_TERMS
	// jb sort the prim terms
	sort_terms_array(order_prim_flat, 0, terms_prim-1, pmax_prim);
	#endif

    int offload_target;
    MPI_Comm_rank(MPI_COMM_WORLD, &offload_target);
	// jb offload pmax_prim and terms_prim
	#pragma offload_transfer target(mic:offload_target) in(pmax_prim) \
	                                     in(order_prim_flat[0:terms_prim*3] : ALLOC RETAIN) \
	                                     in(terms_prim)
    
	return 0;
}

_OFFLOADABLE int get_terms_prim(){
	return terms_prim;
}

_OFFLOADABLE int get_pmax_prim(){
	return pmax_prim;
}

int create_order_prim(){
	order_prim_flat = (int *)malloc(terms_prim*3*sizeof(int));
	
	return 0;
}

_OFFLOADABLE void find_perm_prim(int r,int* p1,int* p2,int* p3){

	int (* restrict order_prim)[3] = (int (*restrict)[3]) order_prim_flat;
	
	*p1 = order_prim[r][0];
	*p2 = order_prim[r][1];
	*p3 = order_prim[r][2];
	return;
}

int set_terms_late(){
	
	int r,i,j,k,n,s;
	int* order_raw = malloc(sizeof(int)*3*MAXLINES);
	int* size = malloc(sizeof(int));
	
	terms_late = alpha_max;
	pmax_late = 0;
	create_order_late();
	int (* restrict order_late)[3] = (int (*restrict)[3]) order_late_flat;
	
	switch(eflag_order_late){
		case 1:
			i=j=k=s=0;
			for (n=0;n<terms_late;n++){
				order_late[n][0] = i;
				order_late[n][1] = j;
				order_late[n][2] = k;
				if(k-j<2&&j-i<2&&k-i<2){
					s++;
					i=j=0;
					k=s;
				} else if(k-j>1){
					k--;
					j++;
				} else {
					i++;
					j=i;
					k=s-2*i;
				
				}
				if (order_late[n][2]>pmax_late)pmax_late = order_late[n][2];
			}
		break;
		case 2:
			i=j=k=s=0;
			for (n=0;n<terms_late;n++){
				order_late[n][0] = i;
				order_late[n][1] = j;
				order_late[n][2] = k;
				if(i<j){
					j--;
					k++;
				} else if(i>0){
					i--;
					j = (s-i)/2;
					k = s-j-i;
				} else {
					s++;
					i = s/3;
					j = (s-i)/2;
					k = s-i-j;
				}
				if (order_late[n][2]>pmax_late)pmax_late = order_late[n][2];
			}
		break;
		case 3:
			terms_late = alpha_max;
			pmax_late = 0;
			//create_order_late(); jb bug
			load_three_int(mdist_file, order_raw, size);
			if(terms_late>*size){
				printf("eigenfunction %d exceeds that in distance file\n",terms_late);
				exit;
			} else {
				for (n=0;n<terms_late;n++){
					order_late[n][0] = order_raw[3*n];
					order_late[n][1] = order_raw[3*n+1];
					order_late[n][2] = order_raw[3*n+2];
					if (order_late[n][2]>pmax_late)pmax_late = order_late[n][2];
				}
			}
		break;
		case 4:
			terms_late = alpha_max;
			pmax_late = 0;
			//create_order_late(); jb bug
			load_three_int(mdist_file, order_raw, size);
			if(terms_late>*size){
				printf("eigenfunction %d exceeds that in distance file\n",terms_late);
				exit;
			} else {
				for (n=0;n<terms_late-1;n++){
					order_late[n][0] = order_raw[3*n];
					order_late[n][1] = order_raw[3*n+1];
					order_late[n][2] = order_raw[3*n+2];
					if (order_late[n][2]>pmax_late)pmax_late = order_late[n][2];
				}
			}
			
			order_late[0][0] = 0;
			order_late[0][1] = 0;
			order_late[0][2] = 0;
			
			for (n=terms_late-1;n>1;n--){
				order_late[n][0] = order_late[n-1][0];
				order_late[n][1] = order_late[n-1][1];
				order_late[n][2] = order_late[n-1][2];
			}
			
			order_late[1][0] = pmax_late+1;
			order_late[1][1] = pmax_late+1;
			order_late[1][2] = pmax_late+2;
		break;
		case 5:
			terms_late = alpha_max;
			pmax_late = 0;
			//create_order_late(); jb bug
			load_three_int(mdist_file, order_raw, size);
			if(terms_late>*size){
				printf("eigenfunction %d exceeds that in distance file\n",terms_late);
				exit;
			} else {
				for (n=0;n<terms_late-3;n++){
					order_late[n][0] = order_raw[3*n];
					order_late[n][1] = order_raw[3*n+1];
					order_late[n][2] = order_raw[3*n+2];
					if (order_late[n][2]>pmax_late)pmax_late = order_late[n][2];
				}
			}
			
			order_late[0][0] = 0;
			order_late[0][1] = 0;
			order_late[0][2] = 0;
			
			for (n=terms_late-1;n>3;n--){
				order_late[n][0] = order_late[n-3][0];
				order_late[n][1] = order_late[n-3][1];
				order_late[n][2] = order_late[n-3][2];
			}
			
// 			order_late[0][0] = 1;
// 			order_late[0][1] = 1;
// 			order_late[0][2] = 1;
			
			order_late[1][0] = pmax_late+1;
			order_late[1][1] = pmax_late+1;
			order_late[1][2] = pmax_late+1;
			
			order_late[2][0] = pmax_late+2;
			order_late[2][1] = pmax_late+2;
			order_late[2][2] = pmax_late+3;
			
			order_late[3][0] = pmax_late+4;
			order_late[3][1] = pmax_late+4;
			order_late[3][2] = pmax_late+5;
		break;
		case 6:
			terms_late = alpha_max;
			pmax_late = 0;
			// create_order_late(); jb bug
			load_three_int(mdist_file, order_raw, size);
			if(terms_late>*size){
				printf("eigenfunction %d exceeds that in distance file\n",terms_late);
				exit;
			} else {
				for (n=0;n<terms_late-4;n++){
					order_late[n][0] = order_raw[3*n];
					order_late[n][1] = order_raw[3*n+1];
					order_late[n][2] = order_raw[3*n+2];
					if (order_late[n][2]>pmax_late)pmax_late = order_late[n][2];
				}
			}
			
			order_late[0][0] = 0;
			order_late[0][1] = 0;
			order_late[0][2] = 0;
			
			for (n=terms_late-1;n>4;n--){
				order_late[n][0] = order_late[n-4][0];
				order_late[n][1] = order_late[n-4][1];
				order_late[n][2] = order_late[n-4][2];
			}
			
// 			order_late[0][0] = 1;
// 			order_late[0][1] = 1;
// 			order_late[0][2] = 1;
			
			order_late[1][0] = pmax_late+1;
			order_late[1][1] = pmax_late+1;
			order_late[1][2] = pmax_late+2;
			
			order_late[2][0] = pmax_late+3;
			order_late[2][1] = pmax_late+5;
			order_late[2][2] = pmax_late+8;
			
			order_late[3][0] = pmax_late+3;
			order_late[3][1] = pmax_late+6;
			order_late[3][2] = pmax_late+7;
			
			order_late[4][0] = pmax_late+4;
			order_late[4][1] = pmax_late+5;
			order_late[4][2] = pmax_late+6;
		break;
		case 7:
			terms_late = alpha_max;
			pmax_late = 0;
			//create_order_late(); jb bug
			load_three_int(mdist_file, order_raw, size);
			if(terms_late>*size){
				printf("eigenfunction %d exceeds that in distance file\n",terms_late);
				exit;
			} else {
				for (n=0;n<terms_late-1;n++){
					order_late[n][0] = order_raw[3*n];
					order_late[n][1] = order_raw[3*n+1];
					order_late[n][2] = order_raw[3*n+2];
					if (order_late[n][2]>pmax_late)pmax_late = order_late[n][2];
				}
			}
			
			order_late[0][0] = 0;
			order_late[0][1] = 0;
			order_late[0][2] = 0;
			
			for (n=terms_late-1;n>1;n--){
				order_late[n][0] = order_late[n-1][0];
				order_late[n][1] = order_late[n-1][1];
				order_late[n][2] = order_late[n-1][2];
			}
			
			order_late[1][0] = pmax_late+1;
			order_late[1][1] = pmax_late+1;
			order_late[1][2] = pmax_late+2;
		break;
		case 8:
			terms_late = alpha_max;
			pmax_late = 0;
			//create_order_late(); jb bug
			load_three_int(mdist_file, order_raw, size);
			if(terms_late>*size){
				printf("eigenfunction %d exceeds that in distance file\n",terms_late);
				exit;
			} else {
				for (n=0;n<terms_late;n++){
					order_late[n][0] = order_raw[3*n];
					order_late[n][1] = order_raw[3*n+1];
					order_late[n][2] = order_raw[3*n+2];
					if (order_late[n][2]>pmax_late)pmax_late = order_late[n][2];
				}
			}
		break;
		case 9:
			terms_late = alpha_max;
			pmax_late = 0;
			//create_order_late(); jb bug
			load_three_int(mdist_file, order_raw, size);
			if(terms_late>*size){
				printf("eigenfunction %d exceeds that in distance file\n",terms_late);
				exit;
			} else {
				for (n=0;n<terms_late-1;n++){
					order_late[n][0] = order_raw[3*n];
					order_late[n][1] = order_raw[3*n+1];
					order_late[n][2] = order_raw[3*n+2];
					if (order_late[n][2]>pmax_late)pmax_late = order_late[n][2];
				}
			}
			
			order_late[0][0] = 0;
			order_late[0][1] = 0;
			order_late[0][2] = 0;
			
			for (n=terms_late-1;n>1;n--){
				order_late[n][0] = order_late[n-1][0];
				order_late[n][1] = order_late[n-1][1];
				order_late[n][2] = order_late[n-1][2];
			}
			
			order_late[1][0] = pmax_late+1;
			order_late[1][1] = pmax_late+1;
			order_late[1][2] = pmax_late+2;
		break;
		default:
			printf("invalid mode specified for eigen function ordering\n");
			exit;
		break;
	}

	pmax_late = 0;
	for (n=0;n<terms_late;n++){
		if(order_late[n][2]>pmax_late) pmax_late = order_late[n][2];
	}

    #ifdef SORT_TERMS
    // jb sort the terms late
	sort_terms_array(order_late_flat, 0, terms_late-1, pmax_late);
    #endif

	// jb offload order_late_flat
	int offload_target;
    MPI_Comm_rank(MPI_COMM_WORLD, &offload_target);
	#pragma offload_transfer target(mic:offload_target) in(pmax_late) in(order_late_flat[0:terms_late*3] : ALLOC RETAIN)
	return 0;
}

int get_terms_late(){
	return terms_late;
}

_OFFLOADABLE int get_pmax_late(){
	return pmax_late;
}

int create_order_late(){
	//order_late = (int **)create_iarray(terms_late,3);
	order_late_flat = (int *) malloc(terms_late*3*sizeof(int));
	return 0;
}

_OFFLOADABLE void find_perm_late(int r,int* p1,int* p2,int* p3){
    	
	//*p1 = order_late[r][0];
	//*p2 = order_late[r][1];
	//*p3 = order_late[r][2];
	int (* restrict order_late)[3] = (int (*restrict)[3]) order_late_flat;
	*p1 = order_late[r][0];
	*p2 = order_late[r][1];
	*p3 = order_late[r][2];
	return;
}

int create_cl(int l_size){

	cl = (double *)create_vector(l_size);
	return 0;
}

int load_cl(int l_size ,int* l_values){
	
	double* cls_data = malloc( sizeof(double)*MAXLINES*2);
	int* cl_len = malloc( sizeof(int));
	load_two(cls_data_file, cls_data, cl_len);
	
	int cl_size = *cl_len;
	
	double cl_raw[cl_size];
	double c_raw[cl_size];
	double pt,clpt;
	
	int i,j;
	
	j=0;
	cl_raw[0] = 0;
	c_raw[0] = 0;
	for (i=0; i<cl_size; i++){
		cl_raw[i] = cls_data[j++];
		c_raw[i] = cls_data[j++]*cl_raw[i]*(cl_raw[i]+1);
// 		printf("%d\t%e\t%e\n", i, cl_raw[i], c_raw[i]);
	}
	
	if (l_values[l_size-1]>cl_raw[cl_size-1]){
		printf("Cls do not contain enough l's, max data %d max cl: %d\n", l_values[l_size-1], (int)cl_raw[cl_size-1]);
		return 1;
		exit;
	}
	
	gsl_spline* spcl =  gsl_spline_alloc (gsl_interp_cspline, cl_size);
	gsl_interp_accel* acccl = gsl_interp_accel_alloc();
	
	gsl_spline_init(spcl,cl_raw,c_raw,cl_size);
	
	for (i=0; i<l_size; i++){
		pt = (double)l_values[i];
		cl[i] = 0;
		if(pt!=0){
			clpt = gsl_spline_eval(spcl,pt,acccl);
			cl[i] = clpt/(pt*(pt+1));
		}
	}
	
	gsl_spline_free(spcl);
	gsl_interp_accel_free(acccl);
	
	return 0;
}

double* get_cl_array() {
	if (clflag_uselens==1) return lens_tt;
	else return cl;
}
double get_cl(int i){
	double value;
	if(clflag_uselens==1){
		value = lens_tt[i];
	}else{
		value = cl[i];
	}
// 	printf("%d\t%e\t%e\n",i,cl[i],value);
	return value;
}

int load_BN(int l_size ,int* l_values){
	
	double* BN_data = malloc( sizeof(double)*MAXLINES*3);
	int* BN_len = malloc( sizeof(int));
	load_three(beam_noise_file, BN_data, BN_len);
	
	int BN_size = *BN_len;
	
	double l_raw[BN_size];
	double beam_raw[BN_size];
	double noise_raw[BN_size];
	
	int i,j;
	double pt;
	
	j=0;
	for (i=0; i<BN_size; i++){
		l_raw[i] = BN_data[j++];
		beam_raw[i] = BN_data[j++];
		noise_raw[i] = BN_data[j++];
// 		printf("%d\t%e\t%e\t%e\n",i,l_raw[i],beam_raw[i],noise_raw[i]);
// 		noise_raw[i] = 0.0;
	}
	if (l_values[l_size-1]>l_raw[BN_size-1]){
		printf("BN do not contain enough l's, max data %d max BN: %d\n", l_values[l_size-1], (int)l_raw[BN_size]);
		return 1;
		exit;
	}
	
	
	
	gsl_spline* spB =  gsl_spline_alloc (gsl_interp_cspline, BN_size);
	gsl_spline* spN =  gsl_spline_alloc (gsl_interp_cspline, BN_size);
	gsl_interp_accel* accB = gsl_interp_accel_alloc();
	gsl_interp_accel* accN = gsl_interp_accel_alloc();
	
	gsl_spline_init(spB,l_raw,beam_raw,BN_size);
	gsl_spline_init(spN,l_raw,noise_raw,BN_size);
	
	for (i=0; i<l_size; i++){
		pt = (double)l_values[i];
		beam[i] = 1.0;
		if(pt!=0){
			beam[i] = gsl_spline_eval(spB,pt,accB);
		}
	}
	
	for (i=0; i<l_size; i++){
		pt = (double)l_values[i];
		noise[i] = 0;
		if(pt!=0){
			noise[i] = gsl_spline_eval(spN,pt,accN);
		}
	}
	
	gsl_spline_free(spB);
	gsl_spline_free(spN);
	gsl_interp_accel_free(accB);
	gsl_interp_accel_free(accN);
	
	return 0;
}

int load_TL(int l_size ,int* l_values){
	
	double* T_data = malloc( sizeof(double)*MAXLINES*2);
	int* T_len = malloc( sizeof(int));
	load_txt_dbl(transfer_wgt_file, 2, T_data, T_len);
	
	int T_size = *T_len;
	
	double l_raw[T_size];
	double Tl_raw[T_size];
	
	int i,j;
	double pt;
	
	j=0;
	for (i=0; i<T_size; i++){
		l_raw[i] = T_data[j++];
		Tl_raw[i] = T_data[j++];
// 		printf("%d\t%d\t%d\t%e\n", i,j,(int)l_raw[i],Tl_raw[i]);
	}
	if (l_values[l_size-1]>l_raw[T_size-1]){
		printf("T do not contain enough l's, max data %d max T: %d\n", l_values[l_size-1], (int)l_raw[T_size]);
		return 1;
		exit;
	}
	
	gsl_spline* spT =  gsl_spline_alloc (gsl_interp_cspline, T_size);
	gsl_interp_accel* accT = gsl_interp_accel_alloc();
	
	gsl_spline_init(spT,l_raw,Tl_raw,T_size);
	
	for (i=0; i<l_size; i++){
		pt = (double)l_values[i];
		t_wgt[i] = 0.0;
		if(pt!=0){
			t_wgt[i] = gsl_spline_eval(spT,pt,accT);
		}
// 		printf("%d\t%d\t%e\n", i,l_values[i],t_wgt[i]);
	}
	
	gsl_spline_free(spT);
	gsl_interp_accel_free(accT);
	
	return 0;
}

int create_beam(int l_size){
	beam = (double *)create_vector(l_size);
	return 0;
}

double* get_beam_array() { return beam; }
double get_beam(int i){
	double value = beam[i];
	return value;
}

int create_noise(int l_size){
	noise = (double *)create_vector(l_size);
	return 0;
}

double* get_noise_array() { return noise; }
double get_noise(int i){
	double value = noise[i];
	return value;
}

int create_t_wgt(int l_size){
	t_wgt = (double *)create_vector(l_size);
	return 0;
}

int load_lens(int l_size ,int* l_values){
	
	double* lens_data = malloc( sizeof(double)*MAXLINES*3);
	int* lens_len = malloc( sizeof(int));
	
	load_txt_dbl(lens_data_file, 3, lens_data, lens_len);
	
	int lens_size = *lens_len;
	int lmax = (int)get_lmax();
// 	printf("lmax %d\t%d\n",lmax,lens_size);
	double l_raw[lens_size];
	double cltt_raw[lens_size];
	double cltp_raw[lens_size];
	
	int i,j;
	double pt;
	
	j=0;
	for (i=0; i<lmax+1; i++){
		l_raw[i] = lens_data[j++];
		cltt_raw[i] = lens_data[j++];
		cltp_raw[i] = lens_data[j++];
	}
	
	if (l_values[l_size-1]>l_raw[lmax]){
		printf("lens do not contain enough l's, max data %d max lens: %d\n", l_values[l_size-1], (int)l_raw[lmax]);
		return 1;
		exit;
	}
	
	
	gsl_spline* sptt =  gsl_spline_alloc (gsl_interp_cspline, lmax+1);
	gsl_spline* sptp =  gsl_spline_alloc (gsl_interp_cspline, lmax+1);
	gsl_interp_accel* acctt = gsl_interp_accel_alloc();
	gsl_interp_accel* acctp = gsl_interp_accel_alloc();
	
	gsl_spline_init(sptt,l_raw,cltt_raw,lmax+1);
	gsl_spline_init(sptp,l_raw,cltp_raw,lmax+1);
	
	for (i=0; i<l_size; i++){
		pt = (double)l_values[i];
		lens_tt[i] = 0.0;
		lens_tp[i] = 0.0;
		if(pt!=0){
			lens_tt[i] = gsl_spline_eval(sptt,pt,acctt);
			lens_tp[i] = gsl_spline_eval(sptp,pt,acctp);
		}
	}
	
	gsl_spline_free(sptt);
	gsl_spline_free(sptp);
	gsl_interp_accel_free(acctt);
	gsl_interp_accel_free(acctp);
	
	return 0;
}

int create_lens(int l_size){
	lens_tt = (double *)create_vector(l_size);
	lens_tp = (double *)create_vector(l_size);
	return 0;
}

double get_lens_tt(int i){
	double value = lens_tt[i];
	return value;
}

double get_lens_tp(int i){
	double value = lens_tp[i];
	return value;
}
	
int create_basis_late(int size, double max, double *vec){

	int i,j,n;
	double c,x,y,z,w;
	
	//basis_late = (double **)create_array(size,pmax_late+1);
	basis_late_flat = (double*) malloc(size * (pmax_late+1) * sizeof *basis_late_flat);
	double (*restrict basis_late)[pmax_late+1] = (double (*restrict)[pmax_late+1]) basis_late_flat;
	
	if(eflag_order_late<4){
		basis_functions_bi(basis_late_flat, size, pmax_late, max, vec);
	}else if(eflag_order_late==4){
		
		double **basis_temp = (double **)create_array(size,pmax_late-1);
	 	basis_functions_fourier(basis_temp, size, pmax_late-2, max, vec);
		
		for (i=0;i<size;i++){
			for (j=0;j<pmax_late-1;j++){
				basis_late[i][j] = basis_temp[i][j];
			}
		}
		destroy_array(basis_temp);
		
		basis_late[0][pmax_late-1] = 0.0;
		basis_late[0][pmax_late] = 0.0;
		
		for (i=1;i<size;i++){
			x = vec[i]*(vec[i]+1e0)/(max*(max+1e0));
			y = (2e0*vec[i]+1e0)/(2e0*max+1e0);
			z = sqrt(x);
			w = pow(y,1e0/6e0);
			basis_late[i][pmax_late-1] = w/z;
			basis_late[i][pmax_late] = w*z;
		}
		
	}else if(eflag_order_late==5){
		double **basis_temp = (double **)create_array(size,pmax_late-4);
	 	basis_functions_fourier(basis_temp, size, pmax_late-5, max, vec);
		
		for (i=0;i<size;i++){
			for (j=0;j<pmax_late-4;j++){
				basis_late[i][j] = basis_temp[i][j];
			}
		}
		destroy_array(basis_temp);
		
		basis_late[0][pmax_late-4] = t_wgt[0];
		basis_late[0][pmax_late-3] = 0.0;
		basis_late[0][pmax_late-2] = 0.0;
		basis_late[0][pmax_late-1] = 0.0;
		basis_late[0][pmax_late] = 0.0;
		
		for (i=1;i<size;i++){
			x = vec[i]*(vec[i]+1.0)/(max*(max+1e0));
			y = (2e0*vec[i]+1.0)/(2e0*max+1.0);
			z = sqrt(x);
			w = pow(y,1.0/6e0);
			basis_late[i][pmax_late-4] = t_wgt[i];
			basis_late[i][pmax_late-3] = w/z;
			basis_late[i][pmax_late-2] = w*z;
			basis_late[i][pmax_late-1] = t_wgt[i]*w/z;
			basis_late[i][pmax_late] = t_wgt[i]*w*z;
		}
	}else if(eflag_order_late==6){
		double **basis_temp = (double **)create_array(size,pmax_late-7);
	 	basis_functions_fourier(basis_temp, size, pmax_late-8, max, vec);
		
		for (i=0;i<size;i++){
			for (j=0;j<pmax_late-7;j++){
				basis_late[i][j] = basis_temp[i][j];
			}
		}
		destroy_array(basis_temp);
		
		basis_late[0][pmax_late-7] = 0.0;
		basis_late[0][pmax_late-6] = 0.0;
		basis_late[0][pmax_late-5] = 0.0;
		basis_late[0][pmax_late-4] = 0.0;
		basis_late[0][pmax_late-3] = 0.0;
		basis_late[0][pmax_late-2] = 0.0;
		basis_late[0][pmax_late-1] = 0.0;
		basis_late[0][pmax_late] = 0.0;
		
		for (i=1;i<size;i++){
			x = vec[i]*(vec[i]+1.0)/(max*(max+1e0));
			y = (2e0*vec[i]+1.0)/(2e0*max+1.0);
			z = sqrt(x);
			w = pow(y,1.0/6e0);
			c = sqrt(get_cl(i)+get_noise(i)/(get_beam(i)*get_beam(i)));
			if(c!=0.0)c = w/c;
			basis_late[i][pmax_late-7] = w/z;
			basis_late[i][pmax_late-6] = w*z;
			basis_late[i][pmax_late-5] = c;
			basis_late[i][pmax_late-4] = c*x;
			basis_late[i][pmax_late-3] = c*get_lens_tt(i);
			basis_late[i][pmax_late-2] = c*get_lens_tp(i);
			basis_late[i][pmax_late-1] = c*x*get_lens_tt(i);
			basis_late[i][pmax_late] = c*x*get_lens_tp(i);
		}
	}else if(eflag_order_late==7){
		
		double **basis_temp = (double **)create_array(size,pmax_late-1);
	 	basis_functions_bi(&basis_temp[0][0], size, pmax_late-2, max, vec);
		
		for (i=0;i<size;i++){
			for (j=0;j<pmax_late-1;j++){
				basis_late[i][j] = basis_temp[i][j];
			}
		}
		destroy_array(basis_temp);
		
		basis_late[0][pmax_late-1] = 0.0;
		basis_late[0][pmax_late] = 0.0;
		
		for (i=1;i<size;i++){
			x = vec[i]*(vec[i]+1e0)/(max*(max+1e0));
			y = (2e0*vec[i]+1e0)/(2e0*max+1e0);
			z = sqrt(x);
			w = pow(y,1e0/6e0);
			basis_late[i][pmax_late-1] = w/z;
			basis_late[i][pmax_late] = w*z;
		}
	}else if(eflag_order_late==8){
		double min = 2e0;
		basis_functions_sinlog(basis_late_flat, size, pmax_late, min, max, vec);
	}else if(eflag_order_late==9){
		
		double **basis_temp = (double **)create_array(size,pmax_late-1);
	 	basis_functions_legendre(basis_temp, size, pmax_late-2, max, vec);
		
		for (i=0;i<size;i++){
			for (j=0;j<pmax_late-1;j++){
				basis_late[i][j] = basis_temp[i][j];
			}
		}
		destroy_array(basis_temp);
		
		basis_late[0][pmax_late-1] = 0.0;
		basis_late[0][pmax_late] = 0.0;
		
		for (i=1;i<size;i++){
			x = vec[i]*(vec[i]+1e0)/(max*(max+1e0));
			y = (2e0*vec[i]+1e0)/(2e0*max+1e0);
			z = sqrt(x);
			w = pow(y,1e0/6e0);
			basis_late[i][pmax_late-1] = w/z;
			basis_late[i][pmax_late] = w*z;
		}
	}

	// JB offload basis_late_flat
	int basis_late_flat_size = size * (pmax_late+1);
	int offload_target;
    MPI_Comm_rank(MPI_COMM_WORLD, &offload_target);
	#pragma offload_transfer target(mic:offload_target) in(basis_late_flat[0:basis_late_flat_size] : ALLOC RETAIN)
	return 0;
}

_OFFLOADABLE double* get_basis_late_array()
{
    return basis_late_flat;
}

double get_basis_late(int i, int p)
{
	double (*restrict basis_late)[pmax_late+1] = (double (*restrict)[pmax_late+1]) basis_late_flat;
	return basis_late[i][p];
}

int create_beta(int l_size, int x_size){
    //int xsize_pad = (b_xsize + 7) & ~7;
	//beta_flat = (double*) _mm_malloc(l_size*(pmax_prim+1)*xsize_pad * sizeof(double), 64);
	beta_flat = (double*) malloc(l_size*(pmax_prim+1)*x_size * sizeof(double));
	
	return 0;
}

_OFFLOADABLE double* get_beta_array() {
    return beta_flat;
}

double get_beta(int l, int n, int i)
{
	double (*restrict beta)[pmax_prim+1][b_xsize] = (double (*restrict)[pmax_prim+1][b_xsize]) beta_flat;
	return beta[l][n][i];
}

int read_beta(){

	int number = 2;
	int *sizes = create_ivector(number);
	ivector_read(&number, proj_size_file, &sizes[0]);
	b_lsize = sizes[0]-1;
	b_xsize = sizes[1]-1;
	b_lvec = create_ivector(b_lsize);
	b_xvec = create_vector(b_xsize);
	int data_size = (b_lsize+1) * (b_xsize+1);
	double **data = create_array(b_lsize+1, b_xsize+1);
	char filename[100];

    //int xsize_pad = (b_xsize + 7) & ~7;
	create_beta(b_lsize, b_xsize);
	int n,l,i;
	double (*restrict beta)[pmax_prim+1][b_xsize] = (double (*restrict)[pmax_prim+1][b_xsize]) beta_flat;
	
	if(eflag_order_prim!=4)
	{
		for(n=0;n<pmax_prim+1;n++)
		{
			char suffix[3] = "";
			sprintf(suffix, "%d", n);
			
			filename[0] = '\0';
			strcat(filename, proj_data_file);
			strcat(filename, "_");
			strcat(filename, suffix);
			
			array_read(&data_size, filename, &data[0][0]);
			
			if(n==0)
			{
				for (i=0; i<b_xsize; i++) b_xvec[i] = data[0][i+1];
				for (l=0; l<b_lsize; l++) b_lvec[l] = (int)data[l+1][0];
				
			}
	
			for(l=0; l<b_lsize; l++) beta[l][n][0] = 0.0;
			
			for(i=0; i<b_xsize; i++)
			{
				for(l=0; l<b_lsize; l++)
				{
					beta[l][n][i] = data[l+1][i+1];
				}
			}
			
		}
	}

	else
	{
		for(n=0;n<pmax_prim-1;n++)
		{
			char suffix[3] = "";
			sprintf(suffix, "%d", n);
			
			filename[0] = '\0';
			strcat(filename, proj_data_file);
			strcat(filename, "_");
			strcat(filename, suffix);
			
			array_read(&data_size, filename, &data[0][0]);
			
			if(n==0)
			{
				for (i=0; i<b_xsize; i++) b_xvec[i] = data[0][i+1];
				for (l=0; l<b_lsize; l++) b_lvec[l] = (int)data[l+1][0];
			}
	
			for(l=0; l<b_lsize; l++) beta[l][n][0] = 0.0;
			
			for(i=0; i<b_xsize; i++)
			{
				for(l=0; l<b_lsize; l++)
				{
					beta[l][n][i] = data[l+1][i+1];
				}
			}
		}
			
		filename[0] = '\0';
		strcat(filename, proj_data_file);
		strcat(filename, "_l1");
			
		array_read(&data_size, filename, &data[0][0]);
	
		for(l=0; l<b_lsize; l++) beta[l][pmax_prim-1][0] = 0.0;
			
		for(i=0; i<b_xsize; i++){
			for(l=0; l<b_lsize; l++){
				beta[l][pmax_prim-1][i] = data[l+1][i+1];
			}
		}
			
		filename[0] = '\0';
		strcat(filename, proj_data_file);
		strcat(filename, "_l2");
			
		array_read(&data_size, filename, &data[0][0]);
	
		for(l=0; l<b_lsize; l++) beta[l][pmax_prim][0] = 0.0;
			
		for(i=0; i<b_xsize; i++){
			for(l=0; l<b_lsize; l++){
				beta[l][pmax_prim][i] = data[l+1][i+1];
			}
		}
	}
	destroy_array(data);

	// jb offload the values of the globals to the MIC
	int beta_flat_size = b_lsize * (b_xsize) * (pmax_prim+1);
	int offload_target;
    MPI_Comm_rank(MPI_COMM_WORLD, &offload_target);
	#pragma offload_transfer target(mic:offload_target) in(b_xvec[0:b_xsize] : ALLOC RETAIN) \
	                                            in(beta_flat[0:beta_flat_size] : ALLOC RETAIN) \
	                                            in(b_lsize) \
                                                in(b_xsize) 
	return 0;
}

_OFFLOADABLE int get_b_lsize()
{
	return b_lsize;
}

_OFFLOADABLE int get_b_xsize()
{
	return b_xsize;
}

int get_b_lvec(int *vec)
{
	int i;
	for(i=0;i<b_lsize;i++){
		vec[i] = b_lvec[i];
	}
	return 0;
}

_OFFLOADABLE int get_b_xvec(double *vec)
{
	int i;
	for(i=0;i<b_xsize;i++)
	{
		vec[i] = b_xvec[i];
	}
	return 0;
}

int create_gamma()
{
	int n = terms_prim;
	gamma_ = (double **)create_array(n,n);
	
	return 0;
}

int update_gamma(double *results)
{
	int i,j;
	
	i = results[0];
	j = results[1];
	gamma_[i][j] = results[2];
	
	return 0;
}

int output_gamma()
{
	int size = terms_prim;
	int gamma_total_size = terms_prim;
	array_write(&gamma_total_size, gamma_data_file, &gamma_[0][0]);
	
	return 0;
}

double get_orthol(int i,int j)
{
	double value = orthol[i][j];
	return value;
}
	
int read_orthol(char *directory)
{
	int size = terms_late;
	int orthol_total_size = size*size;
	orthol = (double **)create_array(size,size);
	array_read(&orthol_total_size, orthol_data_file, &orthol[0][0]);

	return 0;
}

double get_lambdal(int i,int j)
{
	double value = lambdal[i][j];
	return value;
}
	
int read_lambdal(char *directory)
{
	int size = terms_late;
	int lambdal_total_size = size*size;
	lambdal = (double **)create_array(size,size);
	array_read(&lambdal_total_size, lambdal_data_file, &lambdal[0][0]);

	return 0;
}

_OFFLOADABLE double permsix(int i, int j, int k)
{
	double a;
	if(i==k&&i==j){
		a=1.0;
	}else if(i==j||j==k||i==k){
		a=3.0;
	}else{
		a=6.0;
	}
	return a;
}
