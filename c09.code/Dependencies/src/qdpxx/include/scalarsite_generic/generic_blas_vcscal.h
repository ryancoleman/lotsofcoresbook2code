// $Id: generic_blas_vcscal.h,v 1.3 2009-09-15 20:48:42 bjoo Exp $

/*! @file
 *  @brief Generic Scalar VSCAL routine
 *
 */

#ifndef QDP_GENERIC_BLAS_VCSCAL
#define QDP_GENERIC_BLAS_VCSCAL

namespace QDP {

inline
void vcscal(REAL *Out, REAL *scalep, REAL *In, int n_3vec)
{
   double a_r =(REAL)*scalep;
   double a_i =(REAL)*(scalep+1);

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
    // Prefetch whole 3 vec
    i0r = (double)In[inptr++];
    i0i = (double)In[inptr++];
    i1r = (double)In[inptr++];
    i1i = (double)In[inptr++];
    i2r = (double)In[inptr++];
    i2i = (double)In[inptr++];

    int len = 4*n_3vec;
    for(counter = 0; counter < len-1 ; counter++) {
      o0r  = a_r * i0r;
      o0i  = a_i * i0r;
      i0r  = (double)In[inptr++]; // Done with real part get next
      o0r -= a_i * i0i;
      Out[outptr++] = (REAL)o0r;
      o0i += a_r * i0i;
      Out[outptr++] = (REAL)o0i;
      i0i  = (double)In[inptr++];

      o1r  = a_r * i1r;
      o1i  = a_i * i1r;
      i1r  = (double)In[inptr++];
      o1r -= a_i * i1i;
      Out[outptr++] = (REAL)o1r;
      o1i += a_r * i1i;
      Out[outptr++] = (REAL)o1i;
      i1i  = (double)In[inptr++];

      o2r  = a_r * i2r;
      o2i  = a_i * i2r;
      i2r  = (double)In[inptr++];
      o2r -= a_i * i2i;
      Out[outptr++] = (REAL)o2r;
      o2i += a_r * i2i;
      Out[outptr++] = (REAL)o2i;
      i2i  = (double)In[inptr++];
    }

    o0r  = a_r * i0r;
    o0i  = a_i * i0r;
    o0r -= a_i * i0i;
    Out[outptr++] = (REAL)o0r;
    o0i += a_r * i0i;
    Out[outptr++] = (REAL)o0i;
    
    o1r  = a_r * i1r;
    o1i  = a_i * i1r;
    o1r -= a_i * i1i;
    Out[outptr++] = (REAL)o1r;
    o1i += a_r * i1i;
    Out[outptr++] = (REAL)o1i;
    
    o2r  = a_r * i2r;
    o2i  = a_i * i2r;
    o2r -= a_i * i2i;
    Out[outptr++] = (REAL)o2r;
    o2i += a_r * i2i;
    Out[outptr++] = (REAL)o2i;    
  }
}  

} // namespace QDP;

#endif // guard
