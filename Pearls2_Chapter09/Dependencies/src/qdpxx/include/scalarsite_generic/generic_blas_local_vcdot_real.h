// $Id: generic_blas_local_vcdot_real.h,v 1.6 2009-09-15 20:48:41 bjoo Exp $

/*! @file
 *  @brief Generic Scalar, CDOT  routine
 *
 */

#ifndef QDP_GENERIC_BLAS_LOCAL_VCDOT_REAL
#define QDP_GENERIC_BLAS_LOCAL_VCDOT_REAL


namespace QDP {

// Out = Re (< V1, V2 >) = Re(V1^{dagger} V2)
// Out  REAL
// V1 V2 are complex vectors of length 3*n_3vec
// volatile
inline
void l_vcdot_real(DOUBLE *Out, REAL *V1, REAL *V2, int n_3vec)
{

  // routine cleaned up and threaded by Jacques
  double result=0;
  
  double v1_r;
  double v1_i;
  double v2_r;
  double v2_i;

  int counter;

  if( n_3vec > 0 )
  { 
    int len = 24*n_3vec;	// 12*(re,im)

#pragma omp parallel for reduction(+:result) private(v1_r,v1_i,v2_r,v2_i)
    for(counter=0; counter < len; counter+=2)
    {
	    v1_r = (DOUBLE)V1[counter];
	    v1_i = (DOUBLE)V1[counter+1];
  	  v2_r = (DOUBLE)V2[counter];
  	  v2_i = (DOUBLE)V2[counter+1];
  	  
      result += v1_r*v2_r + v1_i*v2_i;
    }
  }
  
  *Out=(DOUBLE)result;
}


} // namespace QDP;

#endif // guard
