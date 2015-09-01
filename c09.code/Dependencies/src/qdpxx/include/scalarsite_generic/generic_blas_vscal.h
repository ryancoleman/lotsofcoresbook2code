// $Id: generic_blas_vscal.h,v 1.3 2009-09-15 20:48:42 bjoo Exp $

/*! @file
 *  @brief Generic Scalar VSCAL routine
 *
 */

#ifndef QDP_GENERIC_BLAS_VSCAL
#define QDP_GENERIC_BLAS_VSCAL

namespace QDP {

// (Vector) out = (Scalar) (*scalep) * (Vector) In
inline
void vscal(REAL *Out, REAL *scalep, REAL *In, int n_3vec)
{
  REAL a = (*scalep);
  int len = 24*n_3vec;
  for(int i=0; i < len; i++) { 
    Out[i] = a*In[i];
  }
#if 0 
   double a = *scalep;

   double i0r;
   double i0i;
   double i1r;
   double i1i;
   double i2r;
   double i2i;

   double o0r;
   double o0i;
   double o1r;
   double o1i;
   double o2r;
   double o2i;

   int counter=0;
   int inptr=0;
   int outptr=0;

  if( n_3vec > 0 ) {
    i0r = (double)In[inptr++];
    i0i = (double)In[inptr++];
    i1r = (double)In[inptr++];
    int len = 4*n_3vec;
    for(counter = 0; counter < len-1 ; counter++) {
      o0r = a*i0r;
      Out[outptr++] = (REAL)o0r;
      
      i1i = (double)In[inptr++];
      i2r = (double)In[inptr++];
      o0i = a*i0i;
      Out[outptr++] = (REAL)o0i;
      
      i2i = (double)In[inptr++];
      i0r = (double)In[inptr++];
      o1r = a*i1r;
      Out[outptr++] = (REAL)o1r;
      
      i0i = (double)In[inptr++];
      i1r = (double)In[inptr++]; // Last prefetched
      
      o1i = a*i1i;
      Out[outptr++] = (REAL)o1i;
      
      o2r= a*i2r;
      Out[outptr++] = (REAL)o2r;
      
      o2i= a*i2i;
      Out[outptr++] = (REAL)o2i;
    }

    o0r = a*i0r;
    Out[outptr++] =(REAL) o0r;
    
    i1i = (double)In[inptr++];
    i2r = (double)In[inptr++];
    o0i = a*i0i;
    Out[outptr++] = (REAL)o0i;
    
    i2i = (double)In[inptr++];
    o1r = a*i1r;
    Out[outptr++] = (REAL)o1r;
    
    o1i = a*i1i;
    Out[outptr++] = (REAL)o1i;
    
    o2r= a*i2r;
    Out[outptr++] = (REAL)o2r;
    
    o2i= a*i2i;
    Out[outptr++] = (REAL)o2i;
    
  }
#endif

}  

} // namespace QDP;

#endif // guard
