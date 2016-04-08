#ifndef QDP_GENERIC_FUSED_SPIN_PROJ_EVALUATES_WRAPPER_H
#define QDP_GENERIC_FUSED_SPIN_PROJ_EVALUATES_WRAPPER_H

////////////////////////////////
// Threading evaluates wrappers
//
// by Xu Guo, EPCC, 28 August, 2008
////////////////////////////////

// user arg for evaluate having order
struct ordered_fused_spin_proj_user_arg{
  const OLattice< SU3Mat >& u;
  const OLattice< FVec >& a;
  OLattice< HVec >& d;
  int base;
  void (*func)(const REAL*, REAL*, unsigned int);
};


// user func for evaluate having order
inline
void ordered_fused_spin_proj_evaluate_function (int lo, int hi, int myId, ordered_fused_spin_proj_user_arg* arg){

  const OLattice< SU3Mat >& u = arg->u;
  const OLattice< FVec >& a = arg->a ;
  OLattice< HVec >& d = arg->d;
  int base = arg->base;
  void (*func)(const REAL*, REAL*,unsigned int) = arg->func; 

  int low = lo + base;
  int high = hi + base;

  for (int site = low; site < high; site++){
     HVec tmp;
     func( (REAL *)&(a.elem(site).elem(0).elem(0).real()),(REAL *)&(tmp.elem(0).elem(0).real()), 1);
      
     _inline_mult_adj_su3_mat_vec(u.elem(site).elem(), tmp.elem(0), d.elem(site).elem(0));
     _inline_mult_adj_su3_mat_vec(u.elem(site).elem(), tmp.elem(1), d.elem(site).elem(1));   
  }

}


// user arg for evaluate NOT having order
struct unordered_fused_spin_proj_user_arg{
  const OLattice< SU3Mat >& u;
  const OLattice< FVec >& a;
  OLattice< HVec >& d;
  const int *tab;
  void (*func)(const REAL*, REAL*,unsigned int);
};


// user func for evaluate NOT having order
inline
void unordered_fused_spin_proj_evaluate_function (int lo, int hi, int myId, unordered_fused_spin_proj_user_arg* arg){

  const OLattice< SU3Mat >& u = arg->u;
  const OLattice< FVec >& a = arg->a ;
  OLattice< HVec >& d = arg->d;
  const int *tab = arg->tab;
  void (*func)(const REAL*, REAL*, unsigned int) = arg->func;

  for (int j = lo; j < hi; j++){
    int site = tab[j];
      
    HVec tmp;
    func( (REAL *)&(a.elem(site).elem(0).elem(0).real()), (REAL *)&(tmp.elem(0).elem(0).real()), 1);
      
      
    _inline_mult_adj_su3_mat_vec(u.elem(site).elem(), tmp.elem(0), d.elem(site).elem(0));
    _inline_mult_adj_su3_mat_vec(u.elem(site).elem(), tmp.elem(1), d.elem(site).elem(1));   
  }

}




#endif
