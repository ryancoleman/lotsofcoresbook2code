// $Id: generic_blas_vaxmby3_norm.h,v 1.2 2007-06-10 14:32:10 edwards Exp $

/*! @file
 *  @brief Generic Scalar VAXPY NORM  routine
 *
 */

#ifndef QDP_GENERIC_BLAS_VAXMBY3_NORM
#define QDP_GENERIC_BLAS_VAXMBY3_NORM

namespace QDP {

// (Vector) Out = Scalar(*ap) * (Vector) xp - (Scalar)(*bp) * (Vector)yp
// *norm = local sumsq(Out);

inline
void vaxmby3_norm(REAL *Out,REAL *ap,REAL *xp, REAL *bp, REAL *yp, 
		  int n_3vec, REAL *norm)
{
   double a;
   double b;
   double x0r;
   double x0i;
  
   double x1r;
   double x1i;
  
   double x2r;
   double x2i;
  
   double y0r;
   double y0i;
  
   double y1r;
   double y1i;
  
   double y2r;
   double y2i;
  
   double z0r;
   double z0i;
  
   double z1r;
   double z1i;
  
   double z2r;
   double z2i;
   double norm_out=0;

  a = *ap;
  b = *bp;
   int index_x = 0;
   int index_y = 0;
   int index_z = 0;
  
   int counter;
  
  for( counter = 0; counter < n_3vec; counter++) {
    x0r = (double)xp[index_x++];
    y0r = (double)yp[index_y++];
    z0r = a*x0r;
    z0r = z0r - b*y0r;
    Out[index_z++] =(REAL) z0r;
    norm_out += z0r*z0r;

    x0i = (double)xp[index_x++];
    y0i = (double)yp[index_y++];
    z0i = a*x0i;
    z0i = z0i - b*y0i;
    Out[index_z++] =(REAL) z0i;
    norm_out += z0i*z0i;    

    x1r = (double)xp[index_x++];
    y1r = (double)yp[index_y++];
    z1r = a*x1r;
    z1r = z1r - b*y1r;
    Out[index_z++] = (REAL)z1r;
    norm_out += z1r * z1r; 
    
    x1i = (double)xp[index_x++];
    y1i = (double)yp[index_y++];
    z1i = a*x1i;
    z1i = z1i - b*y1i;
    Out[index_z++] = (REAL)z1i;
    norm_out += z1i*z1i;
    
    x2r = (double)xp[index_x++];     
    y2r = (double)yp[index_y++];
    z2r = a*x2r;
    z2r = z2r - b*y2r;
    Out[index_z++] = (REAL)z2r;
    norm_out += z2r*z2r;
    
    x2i = (double)xp[index_x++];
    y2i = (double)yp[index_y++];
    z2i = a*x2i;
    z2i = z2i - b*y2i;  
    Out[index_z++] = (REAL)z2i;
    norm_out += z2i*z2i;    
  }
  *norm=(REAL)norm_out;
}

} // namespace QDP;

#endif // guard
