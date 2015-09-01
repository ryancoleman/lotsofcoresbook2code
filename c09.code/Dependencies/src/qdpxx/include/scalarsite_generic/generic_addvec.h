#define _inline_generic_add_su3_vector(aa,bb,cc) \
{ \
  cc.elem(0) = aa.elem(0) + bb.elem(0); \
  cc.elem(1) = aa.elem(1) + bb.elem(1); \
  cc.elem(2) = aa.elem(2) + bb.elem(2); \
}
