// $Id: generic_blas_local_sumsq.h,v 1.3 2009-09-15 20:48:41 bjoo Exp $

/*! @file
 *  @brief Generic Scalar, local sum squared routine
 *
 */

#ifndef QDP_GENERIC_BLAS_LOCAL_SUMSQ
#define QDP_GENERIC_BLAS_LOCAL_SUMSQ


namespace QDP {

// (Double) (*out) = || (Vector) In ||^2 (local to node)
inline
void local_sumsq(DOUBLE *Out, REAL  *In, int n_3vec)
{

  double result;
  double i1;

  int len = 24*n_3vec;
  
  // QDPIO::cout << "Len = " << len << endl;
  int counter;
  result = 0;

  if( n_3vec > 0 ) { 
#pragma omp parallel for reduction(+:result) private(i1)
    for(counter=0; counter < len; counter++) {
      i1 = (double)In[counter];
      result += i1*i1;
    }
  }
  
  *Out=(DOUBLE)result;
}


} // namespace QDP;

#endif // guard
