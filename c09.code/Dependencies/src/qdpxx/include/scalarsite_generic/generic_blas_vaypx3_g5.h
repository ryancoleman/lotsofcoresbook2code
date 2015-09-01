// $Id: generic_blas_vaypx3_g5.h,v 1.2 2007-06-10 14:32:10 edwards Exp $

/*! @file
 *  @brief Generic Scalar VAXPY routine
 *
 */

#ifndef QDP_GENERIC_BLAS_VAYPX3_G5
#define QDP_GENERIC_BLAS_VAYPX3_G5

namespace QDP {

// (Vector) out = (Vector) Add + (Scalar) (*scalep) * (Vector) P{+} InScale 
inline
void xpayz_g5ProjPlus(REAL *Out,REAL *scalep,REAL *Add, REAL *InScale,int n_4vec)
{
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
  
  for( counter = 0; counter < n_4vec; counter++) {
    // Spin Component 0 (AYPX)
    x0r = (double)Add[index_x++];
    y0r = (double)InScale[index_y++];
    z0r = x0r + a*y0r;
    Out[index_z++] =(REAL) z0r;
  
    x0i = (double)Add[index_x++];  
    y0i = (double)InScale[index_y++];
    z0i = x0i + a*y0i;
    Out[index_z++] =(REAL) z0i;

    x1r = (double)Add[index_x++];    
    y1r = (double)InScale[index_y++];
    z1r = x1r + a*y1r;
    Out[index_z++] = (REAL)z1r;

    x1i = (double)Add[index_x++];    
    y1i = (double)InScale[index_y++];
    z1i = x1i + a*y1i;
    Out[index_z++] = (REAL)z1i;
    
    x2r = (double)Add[index_x++];
    y2r = (double)InScale[index_y++];     
    z2r = x2r + a*y2r;
    Out[index_z++] = (REAL)z2r;
   
    x2i = (double)Add[index_x++]; 
    y2i = (double)InScale[index_y++];
    z2i = x2i + a*y2i;  
    Out[index_z++] = (REAL)z2i;

    // Spin Component 1 (AYPX)
    x0r = (double)Add[index_x++];
    y0r = (double)InScale[index_y++];
    z0r = x0r + a*y0r;
    Out[index_z++] =(REAL) z0r;
  
    x0i = (double)Add[index_x++];  
    y0i = (double)InScale[index_y++];
    z0i = x0i + a*y0i;
    Out[index_z++] =(REAL) z0i;

    x1r = (double)Add[index_x++];    
    y1r = (double)InScale[index_y++];
    z1r = x1r + a*y1r;
    Out[index_z++] = (REAL)z1r;

    x1i = (double)Add[index_x++];    
    y1i = (double)InScale[index_y++];
    z1i = x1i + a*y1i;
    Out[index_z++] = (REAL)z1i;
    
    x2r = (double)Add[index_x++];
    y2r = (double)InScale[index_y++];     
    z2r = x2r + a*y2r;
    Out[index_z++] = (REAL)z2r;
   
    x2i = (double)Add[index_x++]; 
    y2i = (double)InScale[index_y++];
    z2i = x2i + a*y2i;  
    Out[index_z++] = (REAL)z2i;

    index_y +=12; // Skip 2 color vecs (3x2x2=12)
    
    // 2 spin components: Copy...

    x0r = (double)Add[index_x++];
    x0i = (double)Add[index_x++];  
    x1r = (double)Add[index_x++];    
    Out[index_z++] =(REAL) x0r;
    x1i = (double)Add[index_x++];   
    Out[index_z++] =(REAL) x0i;
    x2r = (double)Add[index_x++];   
    Out[index_z++] =(REAL) x1r;
    x2i = (double)Add[index_x++];
    Out[index_z++] =(REAL) x1i; 
    x0r = (double)Add[index_x++];    
    Out[index_z++] =(REAL) x2r;
    x0i = (double)Add[index_x++];  
    Out[index_z++] =(REAL) x2i;
    x1r = (double)Add[index_x++];
    Out[index_z++] =(REAL) x0r;
    x1i = (double)Add[index_x++];    
    Out[index_z++] =(REAL) x0i;
    x2r = (double)Add[index_x++];      
    Out[index_z++] =(REAL) x1r;
    x2i = (double)Add[index_x++];
    Out[index_z++] =(REAL) x1i;
   
    // Cache pulls
    x0r = (double)Add[index_x];
    Out[index_z++] =(REAL) x2r;
    x0i = (double)Add[index_x+1];
    Out[index_z++] =(REAL) x2i;

    
  }
}

inline
void xpayz_g5ProjMinus(REAL *Out,REAL *scalep,REAL *Add, REAL *InScale,int n_4vec)
{
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
  
  for( counter = 0; counter < n_4vec; counter++) {
    index_y += 12;

    // 2 spin components: Copy...
    x0r = (double)Add[index_x++];
    x0i = (double)Add[index_x++];  
    x1r = (double)Add[index_x++];    
    Out[index_z++] =(REAL) x0r;
    x1i = (double)Add[index_x++];   
    Out[index_z++] =(REAL) x0i;
    x2r = (double)Add[index_x++];   
    Out[index_z++] =(REAL) x1r;
    x2i = (double)Add[index_x++];
    Out[index_z++] =(REAL) x1i; 
    x0r = (double)Add[index_x++];    
    Out[index_z++] =(REAL) x2r;
    x0i = (double)Add[index_x++];  
    Out[index_z++] =(REAL) x2i;
    x1r = (double)Add[index_x++];
    Out[index_z++] =(REAL) x0r;
    x1i = (double)Add[index_x++];    
    Out[index_z++] =(REAL) x0i;
    x2r = (double)Add[index_x++];      
    Out[index_z++] =(REAL) x1r;
    x2i = (double)Add[index_x++];
    Out[index_z++] =(REAL) x1i;
   
    // Prefetch
    x0r = (double)Add[index_x++];
    Out[index_z++] =(REAL) x2r;
    x0i = (double)Add[index_x++];
    Out[index_z++] =(REAL) x2i;

    // Spin Component 0 (AYPX)
    y0r = (double)InScale[index_y++];
    z0r = x0r + a*y0r;
    Out[index_z++] =(REAL) z0r;
  
    y0i = (double)InScale[index_y++];
    z0i = x0i + a*y0i;
    Out[index_z++] =(REAL) z0i;

    x1r = (double)Add[index_x++];    
    y1r = (double)InScale[index_y++];
    z1r = x1r + a*y1r;
    Out[index_z++] = (REAL)z1r;

    x1i = (double)Add[index_x++];    
    y1i = (double)InScale[index_y++];
    z1i = x1i + a*y1i;
    Out[index_z++] = (REAL)z1i;
    
    x2r = (double)Add[index_x++];
    y2r = (double)InScale[index_y++];     
    z2r = x2r + a*y2r;
    Out[index_z++] = (REAL)z2r;
   
    x2i = (double)Add[index_x++]; 
    y2i = (double)InScale[index_y++];
    z2i = x2i + a*y2i;  
    Out[index_z++] = (REAL)z2i;

    // Spin Component 1 (AYPX)
    x0r = (double)Add[index_x++];
    y0r = (double)InScale[index_y++];
    z0r = x0r + a*y0r;
    Out[index_z++] =(REAL) z0r;
  
    x0i = (double)Add[index_x++];  
    y0i = (double)InScale[index_y++];
    z0i = x0i + a*y0i;
    Out[index_z++] =(REAL) z0i;

    x1r = (double)Add[index_x++];    
    y1r = (double)InScale[index_y++];
    z1r = x1r + a*y1r;
    Out[index_z++] = (REAL)z1r;

    x1i = (double)Add[index_x++];    
    y1i = (double)InScale[index_y++];
    z1i = x1i + a*y1i;
    Out[index_z++] = (REAL)z1i;
    
    x2r = (double)Add[index_x++];
    y2r = (double)InScale[index_y++];     
    z2r = x2r + a*y2r;
    Out[index_z++] = (REAL)z2r;
   
    x2i = (double)Add[index_x++]; 
    y2i = (double)InScale[index_y++];
    z2i = x2i + a*y2i;  
    Out[index_z++] = (REAL)z2i;


  }
}


// (Vector) out = (Vector) Add + (Scalar) (*scalep) * (Vector) P{+} InScale 
inline
void xmayz_g5ProjPlus(REAL *Out,REAL *scalep,REAL *Add, REAL *InScale,int n_4vec)
{
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
  
  for( counter = 0; counter < n_4vec; counter++) {
    // Spin Component 0 (AYPX)
    x0r = (double)Add[index_x++];
    y0r = (double)InScale[index_y++];
    z0r = x0r - a*y0r;
    Out[index_z++] =(REAL) z0r;
  
    x0i = (double)Add[index_x++];  
    y0i = (double)InScale[index_y++];
    z0i = x0i - a*y0i;
    Out[index_z++] =(REAL) z0i;

    x1r = (double)Add[index_x++];    
    y1r = (double)InScale[index_y++];
    z1r = x1r - a*y1r;
    Out[index_z++] = (REAL)z1r;

    x1i = (double)Add[index_x++];    
    y1i = (double)InScale[index_y++];
    z1i = x1i - a*y1i;
    Out[index_z++] = (REAL)z1i;
    
    x2r = (double)Add[index_x++];
    y2r = (double)InScale[index_y++];     
    z2r = x2r - a*y2r;
    Out[index_z++] = (REAL)z2r;
   
    x2i = (double)Add[index_x++]; 
    y2i = (double)InScale[index_y++];
    z2i = x2i - a*y2i;  
    Out[index_z++] = (REAL)z2i;

    // Spin Component 1 (AYPX)
    x0r = (double)Add[index_x++];
    y0r = (double)InScale[index_y++];
    z0r = x0r - a*y0r;
    Out[index_z++] =(REAL) z0r;
  
    x0i = (double)Add[index_x++];  
    y0i = (double)InScale[index_y++];
    z0i = x0i - a*y0i;
    Out[index_z++] =(REAL) z0i;

    x1r = (double)Add[index_x++];    
    y1r = (double)InScale[index_y++];
    z1r = x1r - a*y1r;
    Out[index_z++] = (REAL)z1r;

    x1i = (double)Add[index_x++];    
    y1i = (double)InScale[index_y++];
    z1i = x1i - a*y1i;
    Out[index_z++] = (REAL)z1i;
    
    x2r = (double)Add[index_x++];
    y2r = (double)InScale[index_y++];     
    z2r = x2r - a*y2r;
    Out[index_z++] = (REAL)z2r;
   
    x2i = (double)Add[index_x++]; 
    y2i = (double)InScale[index_y++];
    z2i = x2i - a*y2i;  
    Out[index_z++] = (REAL)z2i;

    index_y +=12; // Skip 2 color vecs (3x2x2=12)
    
    // 2 spin components: Copy...

    x0r = (double)Add[index_x++];
    x0i = (double)Add[index_x++];  
    x1r = (double)Add[index_x++];    
    Out[index_z++] =(REAL) x0r;
    x1i = (double)Add[index_x++];   
    Out[index_z++] =(REAL) x0i;
    x2r = (double)Add[index_x++];   
    Out[index_z++] =(REAL) x1r;
    x2i = (double)Add[index_x++];
    Out[index_z++] =(REAL) x1i; 
    x0r = (double)Add[index_x++];    
    Out[index_z++] =(REAL) x2r;
    x0i = (double)Add[index_x++];  
    Out[index_z++] =(REAL) x2i;
    x1r = (double)Add[index_x++];
    Out[index_z++] =(REAL) x0r;
    x1i = (double)Add[index_x++];    
    Out[index_z++] =(REAL) x0i;
    x2r = (double)Add[index_x++];      
    Out[index_z++] =(REAL) x1r;
    x2i = (double)Add[index_x++];
    Out[index_z++] =(REAL) x1i;
   
    // Cache pulls
    x0r = (double)Add[index_x];
    Out[index_z++] =(REAL) x2r;
    x0i = (double)Add[index_x+1];
    Out[index_z++] =(REAL) x2i;

    
  }
}

inline
void xmayz_g5ProjMinus(REAL *Out,REAL *scalep,REAL *Add, REAL *InScale,int n_4vec)
{
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
  
  for( counter = 0; counter < n_4vec; counter++) {
    index_y += 12;

    // 2 spin components: Copy...
    x0r = (double)Add[index_x++];
    x0i = (double)Add[index_x++];  
    x1r = (double)Add[index_x++];    
    Out[index_z++] =(REAL) x0r;
    x1i = (double)Add[index_x++];   
    Out[index_z++] =(REAL) x0i;
    x2r = (double)Add[index_x++];   
    Out[index_z++] =(REAL) x1r;
    x2i = (double)Add[index_x++];
    Out[index_z++] =(REAL) x1i; 
    x0r = (double)Add[index_x++];    
    Out[index_z++] =(REAL) x2r;
    x0i = (double)Add[index_x++];  
    Out[index_z++] =(REAL) x2i;
    x1r = (double)Add[index_x++];
    Out[index_z++] =(REAL) x0r;
    x1i = (double)Add[index_x++];    
    Out[index_z++] =(REAL) x0i;
    x2r = (double)Add[index_x++];      
    Out[index_z++] =(REAL) x1r;
    x2i = (double)Add[index_x++];
    Out[index_z++] =(REAL) x1i;
   
    // Prefetch
    x0r = (double)Add[index_x++];
    Out[index_z++] =(REAL) x2r;
    x0i = (double)Add[index_x++];
    Out[index_z++] =(REAL) x2i;

    // Spin Component 0 (AYPX)
    y0r = (double)InScale[index_y++];
    z0r = x0r - a*y0r;
    Out[index_z++] =(REAL) z0r;
  
    y0i = (double)InScale[index_y++];
    z0i = x0i - a*y0i;
    Out[index_z++] =(REAL) z0i;

    x1r = (double)Add[index_x++];    
    y1r = (double)InScale[index_y++];
    z1r = x1r - a*y1r;
    Out[index_z++] = (REAL)z1r;

    x1i = (double)Add[index_x++];    
    y1i = (double)InScale[index_y++];
    z1i = x1i - a*y1i;
    Out[index_z++] = (REAL)z1i;
    
    x2r = (double)Add[index_x++];
    y2r = (double)InScale[index_y++];     
    z2r = x2r - a*y2r;
    Out[index_z++] = (REAL)z2r;
   
    x2i = (double)Add[index_x++]; 
    y2i = (double)InScale[index_y++];
    z2i = x2i - a*y2i;  
    Out[index_z++] = (REAL)z2i;

    // Spin Component 1 (AYPX)
    x0r = (double)Add[index_x++];
    y0r = (double)InScale[index_y++];
    z0r = x0r - a*y0r;
    Out[index_z++] =(REAL) z0r;
  
    x0i = (double)Add[index_x++];  
    y0i = (double)InScale[index_y++];
    z0i = x0i - a*y0i;
    Out[index_z++] =(REAL) z0i;

    x1r = (double)Add[index_x++];    
    y1r = (double)InScale[index_y++];
    z1r = x1r - a*y1r;
    Out[index_z++] = (REAL)z1r;

    x1i = (double)Add[index_x++];    
    y1i = (double)InScale[index_y++];
    z1i = x1i - a*y1i;
    Out[index_z++] = (REAL)z1i;
    
    x2r = (double)Add[index_x++];
    y2r = (double)InScale[index_y++];     
    z2r = x2r - a*y2r;
    Out[index_z++] = (REAL)z2r;
   
    x2i = (double)Add[index_x++]; 
    y2i = (double)InScale[index_y++];
    z2i = x2i - a*y2i;  
    Out[index_z++] = (REAL)z2i;


  }
}




} // namespace QDP;

#endif // guard
