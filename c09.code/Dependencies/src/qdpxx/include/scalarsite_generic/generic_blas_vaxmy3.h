// $Id: generic_blas_vaxmy3.h,v 1.3 2009-09-15 20:48:41 bjoo Exp $

/*! @file
 *  @brief Generic Scalar VAXMY routine
 *
 */

#ifndef QDP_GENERIC_BLAS_VAXMY3
#define QDP_GENERIC_BLAS_VAXMY3

namespace QDP {

// (Vector) Out = (Scalar) (*scalep) * (Vector) InScale - (Vector) Add
inline
void vaxmy3(REAL *Out,REAL *scalep,REAL *InScale, REAL *Sub,int n_3vec)
{

  REAL a = (*scalep);
  int len = n_3vec * 24;
  for(int i=0; i < len; i++) { 
    Out[i] = a*InScale[i] - Sub[i];
  }

#if 0
   double a;
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
  
  a = *scalep;
  
   int index_x = 0;
   int index_y = 0;
   int index_z = 0;
  
   int counter;
  int len = n_3vec*4; // 6 deep unrolled loop 
                      // multiply by 4 to get 24.

  for( counter = 0; counter < len; counter++) {
    x0r = (double)InScale[index_x++];
    y0r = (double)Sub[index_y++];
    z0r = a*x0r - y0r;
    Out[index_z++] = (REAL)z0r;
    
    x0i = (double)InScale[index_x++];
    y0i = (double)Sub[index_y++];
    z0i = a*x0i - y0i;
    Out[index_z++] = (REAL)z0i;
    
    x1r = (double)InScale[index_x++];
    y1r = (double)Sub[index_y++];
    z1r = a*x1r - y1r;
    Out[index_z++] = (REAL)z1r;
    
    x1i = (double)InScale[index_x++];
    y1i = (double)Sub[index_y++];
    z1i = a*x1i - y1i;
    Out[index_z++] = (REAL)z1i;
    
    x2r = (double)InScale[index_x++];     
    y2r = (double)Sub[index_y++];
    z2r = a*x2r - y2r;
    Out[index_z++] = (REAL)z2r;
    
    x2i = (double)InScale[index_x++];
    y2i = (double)Sub[index_y++];
    z2i = a*x2i - y2i;  
    Out[index_z++] = (REAL)z2i;
  }
#endif

}

} // namespace QDP;

#endif // guard
