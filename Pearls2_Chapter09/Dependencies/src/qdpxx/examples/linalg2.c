// #include <stdio.h>
#include "jlab_sse.h"

#define PREFDIST 2  


// Specialization to optimize the case   
//    LatticeColorMatrix = LatticeColorMatrix * LatticeColorMatrix
void QDP_M_eq_M_times_M_jlab(su3* u3_base, su3* u1_base, su3* u2_base,
			     int size)
{
  su3 *u0;
  su3 *u1;
  su3 *u2;
  su3 *u3;
  su3 *u4;
  su3 *u5;
  su3 *u6;
  su3 *u7;
  su3 *u8;
  sse_float _sse_sgn12 ALIGN ={-1.0f,-1.0f,1.0f,1.0f};
  sse_float _sse_sgn13 ALIGN ={-1.0f,1.0f,-1.0f,1.0f};
  sse_float _sse_sgn14 ALIGN ={-1.0f,1.0f,1.0f,-1.0f};
  sse_float _sse_sgn23 ALIGN ={1.0f,-1.0f,-1.0f,1.0f};
  sse_float _sse_sgn24 ALIGN ={1.0f,-1.0f,1.0f,-1.0f};
  sse_float _sse_sgn34 ALIGN ={1.0f,1.0f,-1.0f,-1.0f};
  sse_float _sse_sgn1234 ALIGN ={-1.0f,-1.0f,-1.0f,-1.0f};

  int ix;

//  fprintf(stdout,"call QDP_M_eq_M_times_M\n");

  for(ix=0; ix < size; ix+=2) 
  {
    {
      _prefetch_single(u1_base + (ix+PREFDIST));
      _prefetch_single(u2_base + (ix+PREFDIST));
      _prefetch_single(u3_base + (ix+PREFDIST));
    }
    u4 = u2_base + (ix);
    _sse_pair_load_c1_c2((*(u4)));   /*load 3 colors, first two rows*/
    u3 = u1_base + (ix);
    u5 = u3_base + (ix);

    _sse_su3_multiply(*(u3));
    _sse_pair_store_up_c1_c2((*(u5)));
    {
      _prefetch_single(u1_base + (ix+1+PREFDIST));
      _prefetch_single(u2_base + (ix+1+PREFDIST));
      _prefetch_single(u3_base + (ix+1+PREFDIST));
    }
    u7 = u2_base + (ix+1);
    _sse_pair_load_c3_c1((*(u4)),(*(u7)));
    u6 = u1_base + (ix+1);
    u8 = u3_base + (ix+1);

    _sse_su3_multiply_3x1_2sites(*(u3),*(u6));
    _sse_pair_store_up_c3_c1((*(u5)),(*(u8)));
    _sse_pair_load_c2_c3((*(u7)));
    
    _sse_su3_multiply(*(u6));
    _sse_pair_store_up_c2_c3((*(u8)));
  }
}
