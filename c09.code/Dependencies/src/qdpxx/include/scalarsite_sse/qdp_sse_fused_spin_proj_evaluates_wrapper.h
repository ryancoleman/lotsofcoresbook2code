#ifndef QDP_SSE_FUSED_SPIN_PROJ_EVALUATES_WRAPPER_H
#define QDP_SSE_FUSED_SPIN_PROJ_EVALUATES_WRAPPER_H

////////////////////////////////
// Threading evaluates wrappers
//
// by Xu Guo, EPCC, 17 Oct, 2008
////////////////////////////////

// user arg for evaluate having order
struct ordered_sse_fused_spin_proj_user_arg{
  const OLattice< SU3Mat32 >& u;
  const OLattice< FVec32 >& a;
  OLattice< HVec32 >& d;
  int base;
  void (*func)(const REAL32*, REAL32*, unsigned int);
};


// user func for evaluate having order
inline
void ordered_sse_fused_spin_proj_evaluate_function (int lo, int hi, int myId, ordered_sse_fused_spin_proj_user_arg* arg){

  const OLattice< SU3Mat32 >& u = arg->u;
  const OLattice< FVec32 >& a = arg->a ;
  OLattice< HVec32 >& d = arg->d;
  int base = arg->base;
  void (*func)(const REAL32*, REAL32*,unsigned int) = arg->func; 

  int low = lo + base;
  int high = hi + base;

  for (int site = low; site < high; site++){
    HVec32 tmp ; 
    func( (REAL32 *)&(a.elem(site).elem(0).elem(0).real()),
		     (REAL32 *)&(tmp.elem(0).elem(0).real()),
		     1);
      
    su3_matrixf* um = (su3_matrixf *)&(u.elem(site).elem().elem(0,0).real());
    half_wilson_vectorf *tmph = (half_wilson_vectorf *)&( tmp.elem(0).elem(0).real());
    half_wilson_vectorf *dh = (half_wilson_vectorf *)&( d.elem(site).elem(0).elem(0).real());
      
    intrin_sse_mult_adj_su3_mat_hwvec(um, tmph, dh);
  }

}


// user arg for evaluate NOT having order
struct unordered_sse_fused_spin_proj_user_arg{
  const OLattice< SU3Mat32 >& u;
  const OLattice< FVec32 >& a;
  OLattice< HVec32 >& d;
  const int *tab;
  void (*func)(const REAL32*, REAL32*, unsigned int);
};


// user func for evaluate NOT having order
inline
void unordered_sse_fused_spin_proj_evaluate_function (int lo, int hi, int myId, unordered_sse_fused_spin_proj_user_arg* arg){

  const OLattice< SU3Mat32 >& u = arg->u;
  const OLattice< FVec32 >& a = arg->a ;
  OLattice< HVec32 >& d = arg->d;
  const int *tab = arg->tab;
  void (*func)(const REAL32*, REAL32*, unsigned int) = arg->func;

  for (int j = lo; j < hi; j++){
    int site=tab[j];
    HVec32 tmp ; 
    func( (REAL32 *)&(a.elem(site).elem(0).elem(0).real()),
	  (REAL32 *)&(tmp.elem(0).elem(0).real()),
	  1);
      
    su3_matrixf* um = (su3_matrixf *)&(u.elem(site).elem().elem(0,0).real());
    half_wilson_vectorf *tmph = (half_wilson_vectorf *)&( tmp.elem(0).elem(0).real());
    half_wilson_vectorf *dh = (half_wilson_vectorf *)&( d.elem(site).elem(0).elem(0).real());
      
    intrin_sse_mult_adj_su3_mat_hwvec(um, tmph, dh);
  }

}



#endif
