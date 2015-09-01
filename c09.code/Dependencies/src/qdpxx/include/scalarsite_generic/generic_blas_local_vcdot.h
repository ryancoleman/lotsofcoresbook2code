// $Id: generic_blas_local_vcdot.h,v 1.6 2009-09-15 20:48:41 bjoo Exp $

/*! @file
 *  @brief Generic Scalar, CDOT  routine
 *
 */

#ifndef QDP_GENERIC_BLAS_LOCAL_VCDOT
#define QDP_GENERIC_BLAS_LOCAL_VCDOT


namespace QDP {

// Out = < V1, V2 > = V1^{dagger} V2
// Out is complex: Out_re, Out_im 
// V1 V2 are complex vectors of length 3*n_3vec
//volatile

inline
void l_vcdot(DOUBLE *Out_re, DOUBLE *Out_im, REAL *V1, REAL *V2, int n_3vec)
{

  double result_re=0;
  double result_im=0;

  double v1_r;
  double v1_i;
  double v2_r;
  double v2_i;

  int counter;

  if( n_3vec > 0 )  {

    int len = 24*n_3vec;	// 12*(re,im)
    
#pragma omp parallel for reduction(+:result_re,result_im) private(v1_r,v1_i,v2_r,v2_i)
    for(counter=0; counter < len; counter+=2)  {
      v1_r = (DOUBLE)V1[counter];
      v1_i = (DOUBLE)V1[counter+1];
      v2_r = (DOUBLE)V2[counter];
      v2_i = (DOUBLE)V2[counter+1];
      result_re += v1_r*v2_r + v1_i*v2_i;
      result_im += v1_r*v2_i - v1_i*v2_r;
    }
  }
  
  *Out_re=(DOUBLE)result_re;
  *Out_im=(DOUBLE)result_im;
}

} // namespace QDP;

#endif // guard
