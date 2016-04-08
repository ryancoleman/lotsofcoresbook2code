#ifndef QDP_SCALARSITE_GENERIC_LINALG_WRAPPER_H
#define QDP_SCALARSITE_GENERIC_LINALG_WRAPPER_H

////////////////////////////////
// Threading evaluates wrappers
//
// by Xu Guo, EPCC, 07 August, 2008
////////////////////////////////

// typedef
typedef OLattice<PScalar<PColorMatrix<RComplexFloat, 3> > >       C;
typedef OLattice<PSpinVector<PColorVector<RComplexFloat, 3>, 2> > H;


//! user argument for the evaluate function in the ordered situation
//
struct ordered_linalg_user_arg{
  OLattice<PSpinVector<PColorVector<RComplexFloat, 3>, 2> >& d;
  const C& l;
  const H& r;
  int base;
  ordered_linalg_user_arg(  
			  OLattice<PSpinVector<PColorVector<RComplexFloat, 3>, 2> >& d_,
			  const C& l_,
			  const H& r_,
			  int base_) : d(d_), l(l_),r(r_), base(base_) {}
			  
};

//! user function for the evaluate function in the ordered situation
//
inline
void ordered_linalg_evaluate_userfunc(int lo, int hi, int myId, ordered_linalg_user_arg* a)
{

  OLattice<PSpinVector<PColorVector<RComplexFloat, 3>, 2> >& d = a->d;
  const C& l = a->l;
  const H& r = a->r;
  int base = a->base;

  int low = lo + base;
  int high = hi + base;

   // Ordered Way - loop through sites and save a table lookup
   for(int i=low; i < high; i++) { 
      
     _inline_generic_mult_su3_mat_vec(l.elem(i).elem(),
				      r.elem(i).elem(0),
				      d.elem(i).elem(0));
     _inline_generic_mult_su3_mat_vec(l.elem(i).elem(),
				      r.elem(i).elem(1),
				      d.elem(i).elem(1));
   }
  
}

//! user argument for the evaluate function in the unordered situation
//
struct unordered_linalg_user_arg{
  OLattice<PSpinVector<PColorVector<RComplexFloat, 3>, 2> >& d;
  const C& l;
  const H& r;
  const int* tab;

  unordered_linalg_user_arg(
			    OLattice<PSpinVector<PColorVector<RComplexFloat, 3>, 2> >& d_,
			    const C& l_,
			    const H& r_,
			    const int* tab_) : d(d_), l(l_), r(r_), tab(tab_) {}


};

//! user function for the evaluate function in the unordered situation
//
//template<>
inline
void unordered_linalg_evaluate_userfunc(int lo, int hi, int myId,  unordered_linalg_user_arg* a)
{

  OLattice<PSpinVector<PColorVector<RComplexFloat, 3>, 2> >& d = a->d;
  const C& l = a->l;
  const H& r = a->r;
  const int* tab = a->tab;

   // Unordered Way - do a site table lookup
   for(int j=lo; j < hi; j++) { 
     int i = tab[j];
      
     _inline_generic_mult_su3_mat_vec(l.elem(i).elem(),
				      r.elem(i).elem(0),
				      d.elem(i).elem(0));
     _inline_generic_mult_su3_mat_vec(l.elem(i).elem(),
				      r.elem(i).elem(1),
				      d.elem(i).elem(1));
   } 
}


#endif
