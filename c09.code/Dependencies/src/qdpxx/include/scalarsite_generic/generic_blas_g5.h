#ifndef QDP_GENERIC_BLAS_G5
#define QDP_GENERIC_BLAS_G5

namespace QDP {

// (Vector) out = (Scalar) (*scalep) * (Vector) In
inline
void scal_g5(REAL *Out, REAL *scalep, REAL *In, int n_4vec)
{
   double a;
   double x0r;
   double x0i;
  
   double x1r;
   double x1i;
  
   double x2r;
   double x2i;
  
   double z0r;
   double z0i;
  
   double z1r;
   double z1i;
  
   double z2r;
   double z2i;
  
  a = *scalep;
  
   int index_x = 0;
   int index_z = 0;
  
   int counter;
  
  for( counter = 0; counter < n_4vec; counter++) {
    // Spin Component 0
    x0r = (double)In[index_x++];
    z0r = a*x0r;
    Out[index_z++] =(REAL) z0r;
    
    x0i = (double)In[index_x++];
    z0i = a*x0i;
    Out[index_z++] =(REAL) z0i;
    
    x1r = (double)In[index_x++];
    z1r = a*x1r;
    Out[index_z++] = (REAL)z1r;
    
    x1i = (double)In[index_x++];
    z1i = a*x1i;
    Out[index_z++] = (REAL)z1i;
    
    x2r = (double)In[index_x++];     
    z2r = a*x2r;
    Out[index_z++] = (REAL)z2r;
    
    x2i = (double)In[index_x++];
    z2i = a*x2i;
    Out[index_z++] = (REAL)z2i;

    // Spin Component 1
    x0r = (double)In[index_x++];
    z0r = a*x0r;
    Out[index_z++] =(REAL) z0r;
    
    x0i = (double)In[index_x++];
    z0i = a*x0i;
    Out[index_z++] =(REAL) z0i;
    
    x1r = (double)In[index_x++];
    z1r = a*x1r;
    Out[index_z++] = (REAL)z1r;
    
    x1i = (double)In[index_x++];
    z1i = a*x1i;
    Out[index_z++] = (REAL)z1i;
    
    x2r = (double)In[index_x++];     
    z2r = a*x2r;
    Out[index_z++] = (REAL)z2r;
    
    x2i = (double)In[index_x++];
    z2i = a*x2i;
    Out[index_z++] = (REAL)z2i;

    // Spin Component 2
    x0r = (double)In[index_x++];
    z0r = a*x0r;
    Out[index_z++] =-(REAL) z0r;
    
    x0i = (double)In[index_x++];
    z0i = a*x0i;
    Out[index_z++] =-(REAL) z0i;
    
    x1r = (double)In[index_x++];
    z1r = a*x1r;
    Out[index_z++] =-(REAL)z1r;
    
    x1i = (double)In[index_x++];
    z1i = a*x1i;
    Out[index_z++] =-(REAL)z1i;
    
    x2r = (double)In[index_x++];     
    z2r = a*x2r;
    Out[index_z++] =-(REAL)z2r;
    
    x2i = (double)In[index_x++];
    z2i = a*x2i;
    Out[index_z++] =-(REAL)z2i;

    // Spin Component 3
    x0r = (double)In[index_x++];
    z0r = a*x0r;
    Out[index_z++] =-(REAL) z0r;
    
    x0i = (double)In[index_x++];
    z0i = a*x0i;
    Out[index_z++] =-(REAL) z0i;
    
    x1r = (double)In[index_x++];
    z1r = a*x1r;
    Out[index_z++] =-(REAL)z1r;
    
    x1i = (double)In[index_x++];
    z1i = a*x1i;
    Out[index_z++] =-(REAL)z1i;
    
    x2r = (double)In[index_x++];     
    z2r = a*x2r;
    Out[index_z++] =-(REAL)z2r;
    
    x2i = (double)In[index_x++];
    z2i = a*x2i;
    Out[index_z++] =-(REAL)z2i;
  }
}  

// (Vector) out = (Scalar) (*scalep) * (Vector) InScale + (scalep2)*g5*vector)Add)
inline
void axpbyz_g5(REAL *Out,REAL *scalep,REAL *InScale, REAL *scalep2, REAL *Add,int n_4vec)
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
  }
}

// (Vector) out = (Vector) Add + (Scalar) (*scalep) * (Vector) P{+} InScale 
inline
void xmayz_g5(REAL *Out,REAL *scalep,REAL *Add, REAL *InScale,int n_4vec)
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

    // Spin Component 2 (AYPX)
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


    
  }
}

// (Vector) out = (Scalar) (*scalep) * (Vector) InScale + (scalep2)*g5*vector)Add)
inline
void g5_axmbyz(REAL *Out,REAL *scalep,REAL *InScale, REAL *scalep2, REAL *Add,int n_4vec)
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

    // Spin Component 2 (AXPY3)
    x0r = (double)InScale[index_x++];
    y0r = (double)Add[index_y++];
    z0r = b*y0r;
    z0r -= a*x0r ;
    
    Out[index_z++] =(REAL) z0r;
    
    x0i = (double)InScale[index_x++];
    y0i = (double)Add[index_y++];
    z0i = b*y0i;
    z0i -= a*x0i;   
    Out[index_z++] =(REAL) z0i;
    
    x1r = (double)InScale[index_x++];
    y1r = (double)Add[index_y++];
    z1r = b*y1r;
    z1r -= a*x1r ;
    Out[index_z++] = (REAL)z1r;
    
    x1i = (double)InScale[index_x++];
    y1i = (double)Add[index_y++];
    z1i = b*y1i;
    z1i -= a*x1i;
    Out[index_z++] = (REAL)z1i;
    
    x2r = (double)InScale[index_x++];     
    y2r = (double)Add[index_y++];
    z2r = b*y2r;
    z2r -= a*x2r ;
    Out[index_z++] = (REAL)z2r;
    
    x2i = (double)InScale[index_x++];
    y2i = (double)Add[index_y++];
    z2i =  b*y2i;  
    z2i -= a*x2i ;

    Out[index_z++] = (REAL)z2i;

    // Spin Component 2 (AXPY3)
    x0r = (double)InScale[index_x++];
    y0r = (double)Add[index_y++];
    z0r = b*y0r;
    z0r -= a*x0r ;
    
    Out[index_z++] =(REAL) z0r;

    // Spin Component 3
    x0i = (double)InScale[index_x++];
    y0i = (double)Add[index_y++];
    z0i = b*y0i;
    z0i -= a*x0i;   
    Out[index_z++] =(REAL) z0i;
    
    x1r = (double)InScale[index_x++];
    y1r = (double)Add[index_y++];
    z1r = b*y1r;
    z1r -= a*x1r ;
    Out[index_z++] = (REAL)z1r;
    
    x1i = (double)InScale[index_x++];
    y1i = (double)Add[index_y++];
    z1i = b*y1i;
    z1i -= a*x1i;
    Out[index_z++] = (REAL)z1i;
    
    x2r = (double)InScale[index_x++];     
    y2r = (double)Add[index_y++];
    z2r = b*y2r;
    z2r -= a*x2r ;
    Out[index_z++] = (REAL)z2r;
    
    x2i = (double)InScale[index_x++];
    y2i = (double)Add[index_y++];
    z2i =  b*y2i;  
    z2i -= a*x2i ;

    Out[index_z++] = (REAL)z2i;
  }
}

// (Vector) out = (Scalar) (*scalep) * (Vector) InScale + (scalep2)*ig5*vector)Add)
inline
void axpbyz_ig5(REAL *Out,REAL *scalep,REAL *InScale, REAL *scalep2, REAL *Add,int n_4vec)
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
    x0i = (double)InScale[index_x++];
    y0r = (double)Add[index_y++];
    y0i = (double)Add[index_y++];

    z0r = a*x0r ;
    z0r -= b*y0i;
    Out[index_z++] =(REAL) z0r;
    z0i = a*x0i;
    z0i += b*y0r;
    Out[index_z++] =(REAL) z0i;
    
    x1r = (double)InScale[index_x++];
    x1i = (double)InScale[index_x++];
    y1r = (double)Add[index_y++];
    y1i = (double)Add[index_y++];

    z1r = a*x1r ;
    z1r -= b*y1i;
    Out[index_z++] = (REAL)z1r;
    z1i = a*x1i;
    z1i += b*y1r;
    Out[index_z++] = (REAL)z1i;
    
    x2r = (double)InScale[index_x++];     
    x2i = (double)InScale[index_x++];
    y2r = (double)Add[index_y++];
    y2i = (double)Add[index_y++];

    z2r = a*x2r ;
    z2r -= b*y2i;
    Out[index_z++] = (REAL)z2r;
    z2i = a*x2i ;
    z2i +=  b*y2r;  
    Out[index_z++] = (REAL)z2i;

    // Spin Component 1
    x0r = (double)InScale[index_x++];
    x0i = (double)InScale[index_x++];
    y0r = (double)Add[index_y++];
    y0i = (double)Add[index_y++];

    z0r = a*x0r ;
    z0r -= b*y0i;
    Out[index_z++] =(REAL) z0r;
    z0i = a*x0i;
    z0i += b*y0r;
    Out[index_z++] =(REAL) z0i;
    
    x1r = (double)InScale[index_x++];
    x1i = (double)InScale[index_x++];
    y1r = (double)Add[index_y++];
    y1i = (double)Add[index_y++];

    z1r = a*x1r ;
    z1r -= b*y1i;
    Out[index_z++] = (REAL)z1r;
    z1i = a*x1i;
    z1i += b*y1r;
    Out[index_z++] = (REAL)z1i;
    
    x2r = (double)InScale[index_x++];     
    x2i = (double)InScale[index_x++];
    y2r = (double)Add[index_y++];
    y2i = (double)Add[index_y++];

    z2r = a*x2r ;
    z2r -= b*y2i;
    Out[index_z++] = (REAL)z2r;
    z2i = a*x2i ;
    z2i +=  b*y2r;  
    Out[index_z++] = (REAL)z2i;

    // Spin Component 2
    x0r = (double)InScale[index_x++];
    x0i = (double)InScale[index_x++];
    y0r = (double)Add[index_y++];
    y0i = (double)Add[index_y++];

    z0r = a*x0r ;
    z0r += b*y0i;
    Out[index_z++] =(REAL) z0r;
    z0i = a*x0i;
    z0i -= b*y0r;
    Out[index_z++] =(REAL) z0i;
    
    x1r = (double)InScale[index_x++];
    x1i = (double)InScale[index_x++];
    y1r = (double)Add[index_y++];
    y1i = (double)Add[index_y++];

    z1r = a*x1r ;
    z1r += b*y1i;
    Out[index_z++] = (REAL)z1r;
    z1i = a*x1i;
    z1i -= b*y1r;
    Out[index_z++] = (REAL)z1i;
    
    x2r = (double)InScale[index_x++];     
    x2i = (double)InScale[index_x++];
    y2r = (double)Add[index_y++];
    y2i = (double)Add[index_y++];

    z2r = a*x2r ;
    z2r += b*y2i;
    Out[index_z++] = (REAL)z2r;
    z2i = a*x2i ;
    z2i -=  b*y2r;  
    Out[index_z++] = (REAL)z2i;

    // Spin Component 3
    x0r = (double)InScale[index_x++];
    x0i = (double)InScale[index_x++];
    y0r = (double)Add[index_y++];
    y0i = (double)Add[index_y++];

    z0r = a*x0r ;
    z0r += b*y0i;
    Out[index_z++] =(REAL) z0r;
    z0i = a*x0i;
    z0i -= b*y0r;
    Out[index_z++] =(REAL) z0i;
    
    x1r = (double)InScale[index_x++];
    x1i = (double)InScale[index_x++];
    y1r = (double)Add[index_y++];
    y1i = (double)Add[index_y++];

    z1r = a*x1r ;
    z1r += b*y1i;
    Out[index_z++] = (REAL)z1r;
    z1i = a*x1i;
    z1i -= b*y1r;
    Out[index_z++] = (REAL)z1i;
    
    x2r = (double)InScale[index_x++];     
    x2i = (double)InScale[index_x++];
    y2r = (double)Add[index_y++];
    y2i = (double)Add[index_y++];


    z2r = a*x2r ;
    z2r += b*y2i;
    Out[index_z++] = (REAL)z2r;
    z2i = a*x2i ;
    z2i -=  b*y2r;  
    Out[index_z++] = (REAL)z2i;


  }
}


// (Vector) out = (Scalar) (*scalep) * (Vector) InScale - (scalep2)*ig5*vector)Add)
inline
void axmbyz_ig5(REAL *Out,REAL *scalep,REAL *InScale, REAL *scalep2, REAL *Add,int n_4vec)
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
    x0i = (double)InScale[index_x++];
    y0r = (double)Add[index_y++];
    y0i = (double)Add[index_y++];

    z0r = a*x0r ;
    z0r += b*y0i;
    Out[index_z++] =(REAL) z0r;
    z0i = a*x0i;
    z0i -= b*y0r;
    Out[index_z++] =(REAL) z0i;
    
    x1r = (double)InScale[index_x++];
    x1i = (double)InScale[index_x++];
    y1r = (double)Add[index_y++];
    y1i = (double)Add[index_y++];

    z1r = a*x1r ;
    z1r += b*y1i;
    Out[index_z++] = (REAL)z1r;
    z1i = a*x1i;
    z1i -= b*y1r;
    Out[index_z++] = (REAL)z1i;
    
    x2r = (double)InScale[index_x++];     
    x2i = (double)InScale[index_x++];
    y2r = (double)Add[index_y++];
    y2i = (double)Add[index_y++];

    z2r = a*x2r ;
    z2r += b*y2i;
    Out[index_z++] = (REAL)z2r;
    z2i = a*x2i ;
    z2i -=  b*y2r;  
    Out[index_z++] = (REAL)z2i;

    // Spin Component 1
    x0r = (double)InScale[index_x++];
    x0i = (double)InScale[index_x++];
    y0r = (double)Add[index_y++];
    y0i = (double)Add[index_y++];

    z0r = a*x0r ;
    z0r += b*y0i;
    Out[index_z++] =(REAL) z0r;
    z0i = a*x0i;
    z0i -= b*y0r;
    Out[index_z++] =(REAL) z0i;
    
    x1r = (double)InScale[index_x++];
    x1i = (double)InScale[index_x++];
    y1r = (double)Add[index_y++];
    y1i = (double)Add[index_y++];

    z1r = a*x1r ;
    z1r += b*y1i;
    Out[index_z++] = (REAL)z1r;
    z1i = a*x1i;
    z1i -= b*y1r;
    Out[index_z++] = (REAL)z1i;
    
    x2r = (double)InScale[index_x++];     
    x2i = (double)InScale[index_x++];
    y2r = (double)Add[index_y++];
    y2i = (double)Add[index_y++];

    z2r = a*x2r ;
    z2r += b*y2i;
    Out[index_z++] = (REAL)z2r;
    z2i = a*x2i ;
    z2i -=  b*y2r;  
    Out[index_z++] = (REAL)z2i;

    // Spin Component 2
    x0r = (double)InScale[index_x++];
    x0i = (double)InScale[index_x++];
    y0r = (double)Add[index_y++];
    y0i = (double)Add[index_y++];

    z0r = a*x0r ;
    z0r -= b*y0i;
    Out[index_z++] =(REAL) z0r;
    z0i = a*x0i;
    z0i += b*y0r;
    Out[index_z++] =(REAL) z0i;
    
    x1r = (double)InScale[index_x++];
    x1i = (double)InScale[index_x++];
    y1r = (double)Add[index_y++];
    y1i = (double)Add[index_y++];

    z1r = a*x1r ;
    z1r -= b*y1i;
    Out[index_z++] = (REAL)z1r;
    z1i = a*x1i;
    z1i += b*y1r;
    Out[index_z++] = (REAL)z1i;
    
    x2r = (double)InScale[index_x++];     
    x2i = (double)InScale[index_x++];
    y2r = (double)Add[index_y++];
    y2i = (double)Add[index_y++];

    z2r = a*x2r ;
    z2r -= b*y2i;
    Out[index_z++] = (REAL)z2r;
    z2i = a*x2i ;
    z2i +=  b*y2r;  
    Out[index_z++] = (REAL)z2i;

    // Spin Component 3
    x0r = (double)InScale[index_x++];
    x0i = (double)InScale[index_x++];
    y0r = (double)Add[index_y++];
    y0i = (double)Add[index_y++];

    z0r = a*x0r ;
    z0r -= b*y0i;
    Out[index_z++] =(REAL) z0r;
    z0i = a*x0i;
    z0i += b*y0r;
    Out[index_z++] =(REAL) z0i;
    
    x1r = (double)InScale[index_x++];
    x1i = (double)InScale[index_x++];
    y1r = (double)Add[index_y++];
    y1i = (double)Add[index_y++];

    z1r = a*x1r ;
    z1r -= b*y1i;
    Out[index_z++] = (REAL)z1r;
    z1i = a*x1i;
    z1i += b*y1r;
    Out[index_z++] = (REAL)z1i;
    
    x2r = (double)InScale[index_x++];     
    x2i = (double)InScale[index_x++];
    y2r = (double)Add[index_y++];
    y2i = (double)Add[index_y++];


    z2r = a*x2r ;
    z2r -= b*y2i;
    Out[index_z++] = (REAL)z2r;
    z2i = a*x2i ;
    z2i +=  b*y2r;  
    Out[index_z++] = (REAL)z2i;
  }
}



// (Vector) out = (Vector) InScale + (scalep)*ig5*vector)Add)
inline
void xpayz_ig5(REAL *Out,REAL *scalep, REAL *InScale, REAL *Add,int n_4vec)
{
   double a;
  //  double b;

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

    // Spin Component 0 (AXPY3)
    x0r = (double)InScale[index_x++];
    x0i = (double)InScale[index_x++];
    y0r = (double)Add[index_y++];
    y0i = (double)Add[index_y++];

    z0r = x0r - a*y0i ;
    Out[index_z++] =(REAL) z0r;
    z0i = x0i + a*y0r ;
    Out[index_z++] =(REAL) z0i;
    
    x1r = (double)InScale[index_x++];
    x1i = (double)InScale[index_x++];
    y1r = (double)Add[index_y++];
    y1i = (double)Add[index_y++];

    z1r = x1r - a*y1i ;
    Out[index_z++] = (REAL)z1r;
    z1i = x1i + a*y1r;
    Out[index_z++] = (REAL)z1i;
    
    x2r = (double)InScale[index_x++];     
    x2i = (double)InScale[index_x++];
    y2r = (double)Add[index_y++];
    y2i = (double)Add[index_y++];

    z2r = x2r - a*y2i;
    Out[index_z++] = (REAL)z2r;
    z2i = x2i +  a*y2r;  
    Out[index_z++] = (REAL)z2i;

    // Spin Component 1
    x0r = (double)InScale[index_x++];
    x0i = (double)InScale[index_x++];
    y0r = (double)Add[index_y++];
    y0i = (double)Add[index_y++];

    z0r = x0r - a*y0i ;
    Out[index_z++] =(REAL) z0r;
    z0i = x0i + a*y0r ;
    Out[index_z++] =(REAL) z0i;
    
    x1r = (double)InScale[index_x++];
    x1i = (double)InScale[index_x++];
    y1r = (double)Add[index_y++];
    y1i = (double)Add[index_y++];

    z1r = x1r - a*y1i ;
    Out[index_z++] = (REAL)z1r;
    z1i = x1i + a*y1r;
    Out[index_z++] = (REAL)z1i;
    
    x2r = (double)InScale[index_x++];     
    x2i = (double)InScale[index_x++];
    y2r = (double)Add[index_y++];
    y2i = (double)Add[index_y++];

    z2r = x2r - a*y2i;
    Out[index_z++] = (REAL)z2r;
    z2i = x2i +  a*y2r;  
    Out[index_z++] = (REAL)z2i;

    // Spin Component 2
    x0r = (double)InScale[index_x++];
    x0i = (double)InScale[index_x++];
    y0r = (double)Add[index_y++];
    y0i = (double)Add[index_y++];

    z0r = x0r + a*y0i ;
    Out[index_z++] =(REAL) z0r;
    z0i = x0i - a*y0r ;
    Out[index_z++] =(REAL) z0i;
    
    x1r = (double)InScale[index_x++];
    x1i = (double)InScale[index_x++];
    y1r = (double)Add[index_y++];
    y1i = (double)Add[index_y++];

    z1r = x1r + a*y1i ;
    Out[index_z++] = (REAL)z1r;
    z1i = x1i - a*y1r;
    Out[index_z++] = (REAL)z1i;
    
    x2r = (double)InScale[index_x++];     
    x2i = (double)InScale[index_x++];
    y2r = (double)Add[index_y++];
    y2i = (double)Add[index_y++];

    z2r = x2r + a*y2i;
    Out[index_z++] = (REAL)z2r;
    z2i = x2i -  a*y2r;  
    Out[index_z++] = (REAL)z2i;

    // Spin Component 3
    x0r = (double)InScale[index_x++];
    x0i = (double)InScale[index_x++];
    y0r = (double)Add[index_y++];
    y0i = (double)Add[index_y++];

    z0r = x0r + a*y0i ;
    Out[index_z++] =(REAL) z0r;
    z0i = x0i - a*y0r ;
    Out[index_z++] =(REAL) z0i;
    
    x1r = (double)InScale[index_x++];
    x1i = (double)InScale[index_x++];
    y1r = (double)Add[index_y++];
    y1i = (double)Add[index_y++];

    z1r = x1r + a*y1i ;
    Out[index_z++] = (REAL)z1r;
    z1i = x1i - a*y1r;
    Out[index_z++] = (REAL)z1i;
    
    x2r = (double)InScale[index_x++];     
    x2i = (double)InScale[index_x++];
    y2r = (double)Add[index_y++];
    y2i = (double)Add[index_y++];

    z2r = x2r + a*y2i;
    Out[index_z++] = (REAL)z2r;
    z2i = x2i -  a*y2r;  
    Out[index_z++] = (REAL)z2i;
  }
}

// (Vector) out = (Vector) InScale - (scalep)*ig5*vector)Add)
inline
void xmayz_ig5(REAL *Out,REAL *scalep, REAL *InScale, REAL *Add,int n_4vec)
{
   double a;
  //  double b;

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

    // Spin Component 0 (AXPY3)
    x0r = (double)InScale[index_x++];
    x0i = (double)InScale[index_x++];
    y0r = (double)Add[index_y++];
    y0i = (double)Add[index_y++];

    z0r = x0r + a*y0i ;
    Out[index_z++] =(REAL) z0r;
    z0i = x0i - a*y0r ;
    Out[index_z++] =(REAL) z0i;
    
    x1r = (double)InScale[index_x++];
    x1i = (double)InScale[index_x++];
    y1r = (double)Add[index_y++];
    y1i = (double)Add[index_y++];

    z1r = x1r + a*y1i ;
    Out[index_z++] = (REAL)z1r;
    z1i = x1i - a*y1r;
    Out[index_z++] = (REAL)z1i;
    
    x2r = (double)InScale[index_x++];     
    x2i = (double)InScale[index_x++];
    y2r = (double)Add[index_y++];
    y2i = (double)Add[index_y++];

    z2r = x2r + a*y2i;
    Out[index_z++] = (REAL)z2r;
    z2i = x2i -  a*y2r;  
    Out[index_z++] = (REAL)z2i;

    // Spin Component 1
    x0r = (double)InScale[index_x++];
    x0i = (double)InScale[index_x++];
    y0r = (double)Add[index_y++];
    y0i = (double)Add[index_y++];

    z0r = x0r + a*y0i ;
    Out[index_z++] =(REAL) z0r;
    z0i = x0i - a*y0r ;
    Out[index_z++] =(REAL) z0i;
    
    x1r = (double)InScale[index_x++];
    x1i = (double)InScale[index_x++];
    y1r = (double)Add[index_y++];
    y1i = (double)Add[index_y++];

    z1r = x1r + a*y1i ;
    Out[index_z++] = (REAL)z1r;
    z1i = x1i - a*y1r;
    Out[index_z++] = (REAL)z1i;
    
    x2r = (double)InScale[index_x++];     
    x2i = (double)InScale[index_x++];
    y2r = (double)Add[index_y++];
    y2i = (double)Add[index_y++];

    z2r = x2r + a*y2i;
    Out[index_z++] = (REAL)z2r;
    z2i = x2i -  a*y2r;  
    Out[index_z++] = (REAL)z2i;

    // Spin Component 2
    x0r = (double)InScale[index_x++];
    x0i = (double)InScale[index_x++];
    y0r = (double)Add[index_y++];
    y0i = (double)Add[index_y++];

    z0r = x0r - a*y0i ;
    Out[index_z++] =(REAL) z0r;
    z0i = x0i + a*y0r ;
    Out[index_z++] =(REAL) z0i;
    
    x1r = (double)InScale[index_x++];
    x1i = (double)InScale[index_x++];
    y1r = (double)Add[index_y++];
    y1i = (double)Add[index_y++];

    z1r = x1r - a*y1i ;
    Out[index_z++] = (REAL)z1r;
    z1i = x1i + a*y1r;
    Out[index_z++] = (REAL)z1i;
    
    x2r = (double)InScale[index_x++];     
    x2i = (double)InScale[index_x++];
    y2r = (double)Add[index_y++];
    y2i = (double)Add[index_y++];

    z2r = x2r - a*y2i;
    Out[index_z++] = (REAL)z2r;
    z2i = x2i +  a*y2r;  
    Out[index_z++] = (REAL)z2i;

    // Spin Component 3
    x0r = (double)InScale[index_x++];
    x0i = (double)InScale[index_x++];
    y0r = (double)Add[index_y++];
    y0i = (double)Add[index_y++];

    z0r = x0r - a*y0i ;
    Out[index_z++] =(REAL) z0r;
    z0i = x0i + a*y0r ;
    Out[index_z++] =(REAL) z0i;
    
    x1r = (double)InScale[index_x++];
    x1i = (double)InScale[index_x++];
    y1r = (double)Add[index_y++];
    y1i = (double)Add[index_y++];

    z1r = x1r - a*y1i ;
    Out[index_z++] = (REAL)z1r;
    z1i = x1i + a*y1r;
    Out[index_z++] = (REAL)z1i;
    
    x2r = (double)InScale[index_x++];     
    x2i = (double)InScale[index_x++];
    y2r = (double)Add[index_y++];
    y2i = (double)Add[index_y++];

    z2r = x2r - a*y2i;
    Out[index_z++] = (REAL)z2r;
    z2i = x2i +  a*y2r;  
    Out[index_z++] = (REAL)z2i;
  }
}

} // namespace QDP;

#endif // guard
