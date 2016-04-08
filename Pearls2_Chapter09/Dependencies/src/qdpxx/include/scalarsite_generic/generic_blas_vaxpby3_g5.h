// $Id: generic_blas_vaxpby3_g5.h,v 1.2 2007-06-10 14:32:10 edwards Exp $

/*! @file
 *  @brief Generic Scalar VAXPY routine
 *
 */

#ifndef QDP_GENERIC_BLAS_VAXPBY3_G5
#define QDP_GENERIC_BLAS_VAXPBY3_G5

namespace QDP {
// (Vector) out = (Scalar) (*scalep) * (Vector) InScale + (Vector) P{+} Add
inline
void axpbyz_g5ProjPlus(REAL *Out,REAL *scalep,REAL *InScale, REAL *scalep2, REAL *Add,int n_4vec)
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
  
  a = *scalep;
  b = *scalep2;

   int index_x = 0;
   int index_y = 0;
   int index_z = 0;
  
   int counter;
  
  for( counter = 0; counter < n_4vec; counter++) {
    // Spin Component 0 (AXPY3)
    x0r = (double)InScale[index_x++];
    y0r = (double)Add[index_y++];
    z0r = a*x0r ;
    z0r += b*y0r;
    Out[index_z++] =(REAL) z0r;
    
    x0i = (double)InScale[index_x++];
    y0i = (double)Add[index_y++];
    z0i = a*x0i;
    z0i += b*y0i;
    Out[index_z++] =(REAL) z0i;
    
    x1r = (double)InScale[index_x++];
    y1r = (double)Add[index_y++];
    z1r = a*x1r ;
    z1r += b*y1r;
    Out[index_z++] = (REAL)z1r;
    
    x1i = (double)InScale[index_x++];
    y1i = (double)Add[index_y++];
    z1i = a*x1i;
    z1i += b*y1i;
    Out[index_z++] = (REAL)z1i;
    
    x2r = (double)InScale[index_x++];     
    y2r = (double)Add[index_y++];
    z2r = a*x2r ;
    z2r += b*y2r;
    Out[index_z++] = (REAL)z2r;
    
    x2i = (double)InScale[index_x++];
    y2i = (double)Add[index_y++];
    z2i = a*x2i ;
    z2i +=  b*y2i;  
    Out[index_z++] = (REAL)z2i;

    // Spin Component 1
    x0r = (double)InScale[index_x++];
    y0r = (double)Add[index_y++];
    z0r = a*x0r;
    z0r += b*y0r;
    Out[index_z++] =(REAL) z0r;
    
    x0i = (double)InScale[index_x++];
    y0i = (double)Add[index_y++];
    z0i = a*x0i;
    z0i += b*y0i;
    Out[index_z++] =(REAL) z0i;
    
    x1r = (double)InScale[index_x++];
    y1r = (double)Add[index_y++];
    z1r = a*x1r ;
    z1r += b*y1r;
    Out[index_z++] = (REAL)z1r;
    
    x1i = (double)InScale[index_x++];
    y1i = (double)Add[index_y++];
    z1i = a*x1i;
    z1i += b*y1i;
    Out[index_z++] = (REAL)z1i;
    
    x2r = (double)InScale[index_x++];     
    y2r = (double)Add[index_y++];
    z2r = a*x2r;
    z2r +=  b*y2r;
    Out[index_z++] = (REAL)z2r;
    
    x2i = (double)InScale[index_x++];
    y2i = (double)Add[index_y++];
    z2i = a*x2i;
    z2i += b*y2i;  
    Out[index_z++] = (REAL)z2i;

    // Spin Component 2
    index_y+=12;
    x0r = (double)InScale[index_x++];
    z0r = a*x0r;
    Out[index_z++] =(REAL) z0r;
    
    x0i = (double)InScale[index_x++];
    z0i = a*x0i;
    Out[index_z++] =(REAL) z0i;
    
    x1r = (double)InScale[index_x++];
    z1r = a*x1r;
    Out[index_z++] = (REAL)z1r;
    
    x1i = (double)InScale[index_x++];
    z1i = a*x1i;
    Out[index_z++] = (REAL)z1i;
    
    x2r = (double)InScale[index_x++];     
    z2r = a*x2r;
    Out[index_z++] = (REAL)z2r;
    
    x2i = (double)InScale[index_x++];
    z2i = a*x2i;  
    Out[index_z++] = (REAL)z2i;

    // Spin Component 3
    x0r = (double)InScale[index_x++];
    y0r = (double)Add[index_y];      // Pull into cache
    z0r = a*x0r;
    Out[index_z++] =(REAL) z0r;
    
    x0i = (double)InScale[index_x++];
    z0i = a*x0i;
    Out[index_z++] =(REAL) z0i;
    
    x1r = (double)InScale[index_x++];
    z1r = a*x1r;
    Out[index_z++] = (REAL)z1r;
    
    x1i = (double)InScale[index_x++];
    z1i = a*x1i;
    Out[index_z++] = (REAL)z1i;
    
    x2r = (double)InScale[index_x++];     
    z2r = a*x2r;
    Out[index_z++] = (REAL)z2r;
    
    x2i = (double)InScale[index_x++];
    z2i = a*x2i;  
    Out[index_z++] = (REAL)z2i;

    
  }
}

inline
void axpbyz_g5ProjMinus(REAL *Out,REAL *scalep,REAL *InScale, REAL* scalep2, REAL *Add,int n_4vec)
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
  
  a = *scalep;
  b = *scalep2;

   int index_x = 0;
   int index_y = 0;
   int index_z = 0;
  
   int counter;
  
  for( counter = 0; counter < n_4vec; counter++) {

    // Skip y_index
    index_y += 12;

    // Spin Component 0 (AXPY3)
    x0r = (double)InScale[index_x++];
    Out[index_z++] =(REAL) (a*x0r);
    
    x0i = (double)InScale[index_x++];
    Out[index_z++] =(REAL) (a*x0i);
    
    x1r = (double)InScale[index_x++];
    Out[index_z++] = (REAL) (a*x1r);
    
    x1i = (double)InScale[index_x++];
    Out[index_z++] = (REAL)(a*x1i);
    
    x2r = (double)InScale[index_x++];     
    Out[index_z++] = (REAL)(a * x2r);
    
    x2i = (double)InScale[index_x++];
    Out[index_z++] = (REAL)(a * x2i);

    // Spin Component 1
    x0r = (double)InScale[index_x++];
    Out[index_z++] =(REAL) (a*x0r);
    
    x0i = (double)InScale[index_x++];
    Out[index_z++] =(REAL) (a*x0i);
    
    x1r = (double)InScale[index_x++];
    Out[index_z++] = (REAL) (a*x1r);
    
    x1i = (double)InScale[index_x++];
    Out[index_z++] = (REAL)(a*x1i);
    
    x2r = (double)InScale[index_x++];     
    Out[index_z++] = (REAL)(a * x2r);
    
    x2i = (double)InScale[index_x++];
    Out[index_z++] = (REAL)(a * x2i);

    // Spin Component 2 (AXPY3)
    x0r = (double)InScale[index_x++];
    y0r = (double)Add[index_y++];
    z0r = a*x0r ;
    z0r += b*y0r;
    Out[index_z++] =(REAL) z0r;
    
    x0i = (double)InScale[index_x++];
    y0i = (double)Add[index_y++];
    z0i = a*x0i;
    z0i += b*y0i;
    Out[index_z++] =(REAL) z0i;
    
    x1r = (double)InScale[index_x++];
    y1r = (double)Add[index_y++];
    z1r = a*x1r;
    z1r += b*y1r;
    Out[index_z++] = (REAL)z1r;
    
    x1i = (double)InScale[index_x++];
    y1i = (double)Add[index_y++];
    z1i = a*x1i;
    z1i += b*y1i;
    Out[index_z++] = (REAL)z1i;
    
    x2r = (double)InScale[index_x++];     
    y2r = (double)Add[index_y++];
    z2r = a*x2r;
    z2r += b*y2r;
    Out[index_z++] = (REAL)z2r;
    
    x2i = (double)InScale[index_x++];
    y2i = (double)Add[index_y++];
    z2i = a*x2i;
    z2i += b*y2i;  
    Out[index_z++] = (REAL)z2i;

    // Spin Component 3
    x0r = (double)InScale[index_x++];
    y0r = (double)Add[index_y++];
    z0r = a*x0r;
    z0r += b*y0r;
    Out[index_z++] =(REAL) z0r;
    
    x0i = (double)InScale[index_x++];
    y0i = (double)Add[index_y++];
    z0i = a*x0i;
    z0i += b*y0i;
    Out[index_z++] =(REAL) z0i;
    
    x1r = (double)InScale[index_x++];
    y1r = (double)Add[index_y++];
    z1r = a*x1r;
    z1r += b* y1r;
    Out[index_z++] = (REAL)z1r;
    
    x1i = (double)InScale[index_x++];
    y1i = (double)Add[index_y++];
    z1i = a*x1i;
    z1i += b*y1i;
    Out[index_z++] = (REAL)z1i;
    
    x2r = (double)InScale[index_x++];     
    y2r = (double)Add[index_y++];
    z2r = a*x2r;
    z2r += b*y2r;
    Out[index_z++] = (REAL)z2r;
    
    x2i = (double)InScale[index_x++];
    y2i = (double)Add[index_y++];
    z2i = a*x2i;
    z2i += b*y2i;  
    Out[index_z++] = (REAL)z2i;    
  }
}


// (Vector) out = (Scalar) (*scalep) * (Vector) InScale - (*scalep2) * (Vector) P{+} Add
inline
void axmbyz_g5ProjPlus(REAL *Out,REAL *scalep,REAL *InScale, REAL *scalep2, REAL *Add,int n_4vec)
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
  
  a = *scalep;
  b = *scalep2;

   int index_x = 0;
   int index_y = 0;
   int index_z = 0;
  
   int counter;
  
  for( counter = 0; counter < n_4vec; counter++) {
    // Spin Component 0 (AXPY3)
    x0r = (double)InScale[index_x++];
    y0r = (double)Add[index_y++];
    z0r = a*x0r ;
    z0r -= b*y0r;
    Out[index_z++] =(REAL) z0r;
    
    x0i = (double)InScale[index_x++];
    y0i = (double)Add[index_y++];
    z0i = a*x0i;
    z0i -= b*y0i;
    Out[index_z++] =(REAL) z0i;
    
    x1r = (double)InScale[index_x++];
    y1r = (double)Add[index_y++];
    z1r = a*x1r ;
    z1r -= b*y1r;
    Out[index_z++] = (REAL)z1r;
    
    x1i = (double)InScale[index_x++];
    y1i = (double)Add[index_y++];
    z1i = a*x1i;
    z1i -= b*y1i;
    Out[index_z++] = (REAL)z1i;
    
    x2r = (double)InScale[index_x++];     
    y2r = (double)Add[index_y++];
    z2r = a*x2r ;
    z2r -= b*y2r;
    Out[index_z++] = (REAL)z2r;
    
    x2i = (double)InScale[index_x++];
    y2i = (double)Add[index_y++];
    z2i = a*x2i ;
    z2i -=  b*y2i;  
    Out[index_z++] = (REAL)z2i;

    // Spin Component 1
    x0r = (double)InScale[index_x++];
    y0r = (double)Add[index_y++];
    z0r = a*x0r;
    z0r -= b*y0r;
    Out[index_z++] =(REAL) z0r;
    
    x0i = (double)InScale[index_x++];
    y0i = (double)Add[index_y++];
    z0i = a*x0i;
    z0i -= b*y0i;
    Out[index_z++] =(REAL) z0i;
    
    x1r = (double)InScale[index_x++];
    y1r = (double)Add[index_y++];
    z1r = a*x1r ;
    z1r -= b*y1r;
    Out[index_z++] = (REAL)z1r;
    
    x1i = (double)InScale[index_x++];
    y1i = (double)Add[index_y++];
    z1i = a*x1i;
    z1i -= b*y1i;
    Out[index_z++] = (REAL)z1i;
    
    x2r = (double)InScale[index_x++];     
    y2r = (double)Add[index_y++];
    z2r = a*x2r;
    z2r -=  b*y2r;
    Out[index_z++] = (REAL)z2r;
    
    x2i = (double)InScale[index_x++];
    y2i = (double)Add[index_y++];
    z2i = a*x2i;
    z2i -= b*y2i;  
    Out[index_z++] = (REAL)z2i;

    // Spin Component 2
    index_y+=12;
    x0r = (double)InScale[index_x++];
    z0r = a*x0r;
    Out[index_z++] =(REAL) z0r;
    
    x0i = (double)InScale[index_x++];
    z0i = a*x0i;
    Out[index_z++] =(REAL) z0i;
    
    x1r = (double)InScale[index_x++];
    z1r = a*x1r;
    Out[index_z++] = (REAL)z1r;
    
    x1i = (double)InScale[index_x++];
    z1i = a*x1i;
    Out[index_z++] = (REAL)z1i;
    
    x2r = (double)InScale[index_x++];     
    z2r = a*x2r;
    Out[index_z++] = (REAL)z2r;
    
    x2i = (double)InScale[index_x++];
    z2i = a*x2i;  
    Out[index_z++] = (REAL)z2i;

    // Spin Component 3
    x0r = (double)InScale[index_x++];
    y0r = (double)Add[index_y];      // Pull into cache
    z0r = a*x0r;
    Out[index_z++] =(REAL) z0r;
    
    x0i = (double)InScale[index_x++];
    z0i = a*x0i;
    Out[index_z++] =(REAL) z0i;
    
    x1r = (double)InScale[index_x++];
    z1r = a*x1r;
    Out[index_z++] = (REAL)z1r;
    
    x1i = (double)InScale[index_x++];
    z1i = a*x1i;
    Out[index_z++] = (REAL)z1i;
    
    x2r = (double)InScale[index_x++];     
    z2r = a*x2r;
    Out[index_z++] = (REAL)z2r;
    
    x2i = (double)InScale[index_x++];
    z2i = a*x2i;  
    Out[index_z++] = (REAL)z2i;

    
  }
}

// OUt = scalep*InScale - scalep2*Add
inline
void axmbyz_g5ProjMinus(REAL *Out,REAL *scalep,REAL *InScale, REAL* scalep2, REAL *Add,int n_4vec)
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
  
  a = *scalep;
  b = *scalep2;

   int index_x = 0;
   int index_y = 0;
   int index_z = 0;
  
   int counter;
  
  for( counter = 0; counter < n_4vec; counter++) {

    // Skip y_index
    index_y += 12;

    // Spin Component 0 (AXPY3)
    x0r = (double)InScale[index_x++];
    Out[index_z++] =(REAL) (a*x0r);
    
    x0i = (double)InScale[index_x++];
    Out[index_z++] =(REAL) (a*x0i);
    
    x1r = (double)InScale[index_x++];
    Out[index_z++] = (REAL) (a*x1r);
    
    x1i = (double)InScale[index_x++];
    Out[index_z++] = (REAL)(a*x1i);
    
    x2r = (double)InScale[index_x++];     
    Out[index_z++] = (REAL)(a * x2r);
    
    x2i = (double)InScale[index_x++];
    Out[index_z++] = (REAL)(a * x2i);

    // Spin Component 1
    x0r = (double)InScale[index_x++];
    Out[index_z++] =(REAL) (a*x0r);
    
    x0i = (double)InScale[index_x++];
    Out[index_z++] =(REAL) (a*x0i);
    
    x1r = (double)InScale[index_x++];
    Out[index_z++] = (REAL) (a*x1r);
    
    x1i = (double)InScale[index_x++];
    Out[index_z++] = (REAL)(a*x1i);
    
    x2r = (double)InScale[index_x++];     
    Out[index_z++] = (REAL)(a * x2r);
    
    x2i = (double)InScale[index_x++];
    Out[index_z++] = (REAL)(a * x2i);

    // Spin Component 2 (AXPY3)
    x0r = (double)InScale[index_x++];
    y0r = (double)Add[index_y++];
    z0r = a*x0r ;
    z0r -= b*y0r;
    Out[index_z++] =(REAL) z0r;
    
    x0i = (double)InScale[index_x++];
    y0i = (double)Add[index_y++];
    z0i = a*x0i;
    z0i -= b*y0i;
    Out[index_z++] =(REAL) z0i;
    
    x1r = (double)InScale[index_x++];
    y1r = (double)Add[index_y++];
    z1r = a*x1r;
    z1r -= b*y1r;
    Out[index_z++] = (REAL)z1r;
    
    x1i = (double)InScale[index_x++];
    y1i = (double)Add[index_y++];
    z1i = a*x1i;
    z1i -= b*y1i;
    Out[index_z++] = (REAL)z1i;
    
    x2r = (double)InScale[index_x++];     
    y2r = (double)Add[index_y++];
    z2r = a*x2r;
    z2r -= b*y2r;
    Out[index_z++] = (REAL)z2r;
    
    x2i = (double)InScale[index_x++];
    y2i = (double)Add[index_y++];
    z2i = a*x2i;
    z2i -= b*y2i;  
    Out[index_z++] = (REAL)z2i;

    // Spin Component 3
    x0r = (double)InScale[index_x++];
    y0r = (double)Add[index_y++];
    z0r = a*x0r;
    z0r -= b*y0r;
    Out[index_z++] =(REAL) z0r;
    
    x0i = (double)InScale[index_x++];
    y0i = (double)Add[index_y++];
    z0i = a*x0i;
    z0i -= b*y0i;
    Out[index_z++] =(REAL) z0i;
    
    x1r = (double)InScale[index_x++];
    y1r = (double)Add[index_y++];
    z1r = a*x1r;
    z1r -= b* y1r;
    Out[index_z++] = (REAL)z1r;
    
    x1i = (double)InScale[index_x++];
    y1i = (double)Add[index_y++];
    z1i = a*x1i;
    z1i -= b*y1i;
    Out[index_z++] = (REAL)z1i;
    
    x2r = (double)InScale[index_x++];     
    y2r = (double)Add[index_y++];
    z2r = a*x2r;
    z2r -= b*y2r;
    Out[index_z++] = (REAL)z2r;
    
    x2i = (double)InScale[index_x++];
    y2i = (double)Add[index_y++];
    z2i = a*x2i;
    z2i -= b*y2i;  
    Out[index_z++] = (REAL)z2i;    
  }
}




} // namespace QDP;

#endif // guard
