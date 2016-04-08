// This recoded in generic_mat_vec.h

#if 0
#define _inline_generic_mult_adj_su3_mat_vec(aa,bb,cc) \
{ \
  cc.elem(0) = adjMultiply(aa.elem(0,0),bb.elem(0))  \
             + adjMultiply(aa.elem(1,0),bb.elem(1))  \
             + adjMultiply(aa.elem(2,0),bb.elem(2)); \
  cc.elem(1) = adjMultiply(aa.elem(0,1),bb.elem(0))  \
             + adjMultiply(aa.elem(1,1),bb.elem(1))  \
             + adjMultiply(aa.elem(2,1),bb.elem(2)); \
  cc.elem(2) = adjMultiply(aa.elem(0,2),bb.elem(0))  \
             + adjMultiply(aa.elem(1,2),bb.elem(1))  \
             + adjMultiply(aa.elem(2,2),bb.elem(2)); \
}
#endif
