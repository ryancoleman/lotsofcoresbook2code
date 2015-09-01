// $Id: generic_blas_vcaxpy3.h,v 1.3 2009-09-15 20:48:42 bjoo Exp $

/*! @file
 *  @brief Generic Scalar VAXPY routine
 *
 */

#ifndef QDP_GENERIC_BLAS_VCAXPY3
#define QDP_GENERIC_BLAS_VCAXPY3

namespace QDP {
// (Vector) out = (Complex) (*scalep) * (Vector) InScale + (Vector) Add
inline
void vcaxpy3(REAL *Out,REAL *scalep,REAL *InScale, REAL *Add,int n_3vec)
{
   double a_r;
   double a_i;

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
  
  a_r =(double)(*scalep);
  a_i =(double)*(scalep+1);
  
   int index_x = 0;
   int index_y = 0;
   int index_z = 0;
  
   int counter;

  if( n_3vec > 0 ) { 
    int len = n_3vec*4;

    // Prefetch whole vectors
    x0r = (double)InScale[index_x++];    
    x0i = (double)InScale[index_x++];
    y0r = (double)Add[index_y++];
    y0i = (double)Add[index_y++];

    x1r = (double)InScale[index_x++];    
    x1i = (double)InScale[index_x++];
    y1r = (double)Add[index_y++];
    y1i = (double)Add[index_y++];
 
    x2r = (double)InScale[index_x++];    
    x2i = (double)InScale[index_x++];
    y2r = (double)Add[index_y++];

    
    for( counter = 0; counter < len-1; counter++) {
      y2i = (double)Add[index_y++];
      z0r = a_r * x0r + y0r;    
      z0i = a_i * x0r + y0i;
      x0r = (double)InScale[index_x++];    
      y0r = (double)Add[index_y++];     
      z0r -= a_i * x0i;
      Out[index_z++] = (REAL)z0r;
      z0i += a_r * x0i;
      Out[index_z++] = (REAL)z0i;
      x0i = (double)InScale[index_x++];



      z1r = a_r * x1r + y1r;
      y0i = (double)Add[index_y++];      
      z1i = a_i * x1r + y1i;
      x1r = (double)InScale[index_x++];    
      y1r = (double)Add[index_y++];
      z1r -= a_i * x1i;
      Out[index_z++] = (REAL)z1r;
      z1i += a_r * x1i;
      Out[index_z++] = (REAL)z1i;
      x1i = (double)InScale[index_x++];

      z2r = a_r * x2r + y2r;
      y1i = (double)Add[index_y++];
      z2i = a_i * x2r + y2i;
      x2r = (double)InScale[index_x++];    
      y2r = (double)Add[index_y++];
      z2r -= a_i * x2i;
      Out[index_z++] = (REAL)z2r;
      z2i += a_r * x2i;
      Out[index_z++] = (REAL)z2i;
      x2i = (double)InScale[index_x++];
     
    }

    y2i = (double)Add[index_y++];
    z0r = a_r * x0r + y0r;
    z0i = a_i * x0r + y0i;       
    z0r -= a_i * x0i;
    Out[index_z++] = (REAL)z0r;
    z0i += a_r * x0i;
    Out[index_z++] = (REAL)z0i;

    z1r = a_r * x1r + y1r;
    z1i = a_i * x1r + y1i;        
    z1r -= a_i * x1i;
    Out[index_z++]= (REAL)z1r;
    z1i += a_r * x1i;
    Out[index_z++]= (REAL)z1i;

    z2r = a_r * x2r + y2r;
    z2i = a_i * x2r + y2i;
    z2r -= a_i * x2i;
    Out[index_z++]= (REAL)z2r;
    z2i += a_r * x2i;
    Out[index_z++]= (REAL)z2i;
    
  }
}


} // namespace QDP;

#endif // guard
