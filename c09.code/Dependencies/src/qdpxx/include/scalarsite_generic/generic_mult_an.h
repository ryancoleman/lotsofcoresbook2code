#define _inline_generic_mult_su3_an(aa,bb,cc) \
{ \
  cc.elem(0,0) = adjMultiply(aa.elem(0,0),bb.elem(0,0))  \
               + adjMultiply(aa.elem(1,0),bb.elem(1,0))  \
               + adjMultiply(aa.elem(2,0),bb.elem(2,0)); \
  cc.elem(1,0) = adjMultiply(aa.elem(0,1),bb.elem(0,0))  \
               + adjMultiply(aa.elem(1,1),bb.elem(1,0))  \
               + adjMultiply(aa.elem(2,1),bb.elem(2,0)); \
  cc.elem(2,0) = adjMultiply(aa.elem(0,2),bb.elem(0,0))  \
               + adjMultiply(aa.elem(1,2),bb.elem(1,0))  \
               + adjMultiply(aa.elem(2,2),bb.elem(2,0)); \
  \
  cc.elem(0,1) = adjMultiply(aa.elem(0,0),bb.elem(0,1))  \
               + adjMultiply(aa.elem(1,0),bb.elem(1,1))  \
               + adjMultiply(aa.elem(2,0),bb.elem(2,1)); \
  cc.elem(1,1) = adjMultiply(aa.elem(0,1),bb.elem(0,1))  \
               + adjMultiply(aa.elem(1,1),bb.elem(1,1))  \
               + adjMultiply(aa.elem(2,1),bb.elem(2,1)); \
  cc.elem(2,1) = adjMultiply(aa.elem(0,2),bb.elem(0,1))  \
               + adjMultiply(aa.elem(1,2),bb.elem(1,1))  \
               + adjMultiply(aa.elem(2,2),bb.elem(2,1)); \
  \
  cc.elem(0,2) = adjMultiply(aa.elem(0,0),bb.elem(0,2))  \
               + adjMultiply(aa.elem(1,0),bb.elem(1,2))  \
               + adjMultiply(aa.elem(2,0),bb.elem(2,2)); \
  cc.elem(1,2) = adjMultiply(aa.elem(0,1),bb.elem(0,2))  \
               + adjMultiply(aa.elem(1,1),bb.elem(1,2))  \
               + adjMultiply(aa.elem(2,1),bb.elem(2,2)); \
  cc.elem(2,2) = adjMultiply(aa.elem(0,2),bb.elem(0,2))  \
               + adjMultiply(aa.elem(1,2),bb.elem(1,2))  \
               + adjMultiply(aa.elem(2,2),bb.elem(2,2)); \
}
