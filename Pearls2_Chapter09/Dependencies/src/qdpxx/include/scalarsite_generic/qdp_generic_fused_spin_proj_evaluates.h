#ifndef QDP_GENERIC_FUSED_SPIN_PROJ_EVALUATES_H
#define QDP_GENERIC_FUSED_SPIN_PROJ_EVALUATES_H


/* Evaluates for things like adj(u)*spinProjectDir0Plus(y) */
using namespace QDP;
namespace QDP {

typedef PScalar< PColorMatrix< RComplex<REAL>, 3> > SU3Mat;



////////////////////////////////
// Threading evaluates
//
// by Xu Guo, EPCC, 28 August, 2008
////////////////////////////////

// ther wrappers for the functions to be threaded
#include "qdp_generic_fused_spin_proj_evaluates_wrapper.h"


// HalfVec = adj(u)*SpinProjectDir0Plus(Vec);
template<>
inline
void evaluate(OLattice< HVec >& d,
              const OpAssign& op,
              const QDPExpr<
	              BinaryNode< 
	                FnAdjMultSprojDir0Plus,
	               
	                UnaryNode< OpIdentity, 
	                  Reference< QDPType< SU3Mat, OLattice< SU3Mat > > > >,
	               
	                  UnaryNode< OpIdentity,
                          Reference< QDPType< FVec,   OLattice< FVec > > > >
                      >,
	              OLattice< HVec > 
                    >&rhs,
	      const Subset& s)
{
  const OLattice< SU3Mat >& u = static_cast< const OLattice< SU3Mat >& >(rhs.expression().left().child());
  const OLattice< FVec >& a = static_cast< const OLattice< FVec >& >(rhs.expression().right().child());


  if( s.hasOrderedRep() ) { 

    int totalSize = s.end() - s.start() +1;

    ordered_fused_spin_proj_user_arg arg = {u, a, d, s.start(), inlineSpinProjDir0Plus};

    dispatch_to_threads(totalSize, arg, ordered_fused_spin_proj_evaluate_function);

    ///////////////////
    // Original code
    ////////////////////
    /*
    for(int site = s.start(); site <= s.end(); site++) { 
      HVec tmp;
      inlineSpinProjDir0Plus( (REAL *)&(a.elem(site).elem(0).elem(0).real()),
			      (REAL *)&(tmp.elem(0).elem(0).real()),
			      1);
      
      
      _inline_mult_adj_su3_mat_vec(u.elem(site).elem(), tmp.elem(0), d.elem(site).elem(0));
      _inline_mult_adj_su3_mat_vec(u.elem(site).elem(), tmp.elem(1), d.elem(site).elem(1));   
 
      }*/
  }
  else { 
    const int *tab = s.siteTable().slice();

    int totalSize = s.numSiteTable();

    unordered_fused_spin_proj_user_arg arg = {u, a, d, tab, inlineSpinProjDir0Plus};

    dispatch_to_threads(totalSize, arg, unordered_fused_spin_proj_evaluate_function);

    ///////////////////
    // Original code
    ////////////////////
    /*
    for(int j=0; j < s.numSiteTable(); j++) { 
      int site = tab[j];
      
      HVec tmp;
      inlineSpinProjDir0Plus( (REAL *)&(a.elem(site).elem(0).elem(0).real()),
			      (REAL *)&(tmp.elem(0).elem(0).real()),
			      1);
      
      
      _inline_mult_adj_su3_mat_vec(u.elem(site).elem(), tmp.elem(0), d.elem(site).elem(0));
      _inline_mult_adj_su3_mat_vec(u.elem(site).elem(), tmp.elem(1), d.elem(site).elem(1));   
 
      }*/
  }
}

// HalfVec = adj(u)*SpinProjectDir0Minus(Vec);
template<>
inline
void evaluate(OLattice< HVec >& d,
              const OpAssign& op,
              const QDPExpr<
	              BinaryNode< 
	                FnAdjMultSprojDir0Minus,
	               
	                UnaryNode< OpIdentity, 
	                  Reference< QDPType< SU3Mat, OLattice< SU3Mat > > > >,
	               
	                  UnaryNode< OpIdentity,
                          Reference< QDPType< FVec,   OLattice< FVec > > > >
                      >,
	              OLattice< HVec > 
                    >&rhs,
	      const Subset& s)
{
  const OLattice< SU3Mat >& u = static_cast< const OLattice< SU3Mat >& >(rhs.expression().left().child());
  const OLattice< FVec >& a = static_cast< const OLattice< FVec >& >(rhs.expression().right().child());

  if( s.hasOrderedRep() ) { 
    int totalSize = s.end() - s.start() +1;

    ordered_fused_spin_proj_user_arg arg = {u, a, d, s.start(), inlineSpinProjDir0Minus};

    dispatch_to_threads(totalSize, arg, ordered_fused_spin_proj_evaluate_function);

    ///////////////////
    // Original code
    ////////////////////
    /*
    for(int site = s.start(); site <= s.end(); site++) { 
      HVec tmp;
      inlineSpinProjDir0Minus( (REAL *)&(a.elem(site).elem(0).elem(0).real()),
			       (REAL *)&(tmp.elem(0).elem(0).real()),
			       1);
      
      
      _inline_mult_adj_su3_mat_vec(u.elem(site).elem(), tmp.elem(0), d.elem(site).elem(0));
      _inline_mult_adj_su3_mat_vec(u.elem(site).elem(), tmp.elem(1), d.elem(site).elem(1));   
      
      }*/
  }
  else { 

    const int *tab = s.siteTable().slice();

    int totalSize = s.numSiteTable();

    unordered_fused_spin_proj_user_arg arg = {u, a, d, tab, inlineSpinProjDir0Minus};

    dispatch_to_threads(totalSize, arg, unordered_fused_spin_proj_evaluate_function);

    ///////////////////
    // Original code
    ////////////////////
    /*
    for(int j=0; j < s.numSiteTable(); j++) { 
      int site = tab[j];
      
      
      HVec tmp;
      inlineSpinProjDir0Minus( (REAL *)&(a.elem(site).elem(0).elem(0).real()),
			       (REAL *)&(tmp.elem(0).elem(0).real()),
			       1);
      
      
      _inline_mult_adj_su3_mat_vec(u.elem(site).elem(), tmp.elem(0), d.elem(site).elem(0));
      _inline_mult_adj_su3_mat_vec(u.elem(site).elem(), tmp.elem(1), d.elem(site).elem(1));   
      
      }*/
  }

}

// HalfVec = adj(u)*SpinProjectDir1Plus(Vec);
template<>
inline
void evaluate(OLattice< HVec >& d,
              const OpAssign& op,
              const QDPExpr<
	              BinaryNode< 
	                FnAdjMultSprojDir1Plus,
	               
	                UnaryNode< OpIdentity, 
	                  Reference< QDPType< SU3Mat, OLattice< SU3Mat > > > >,
	               
	                  UnaryNode< OpIdentity,
                          Reference< QDPType< FVec,   OLattice< FVec > > > >
                      >,
	              OLattice< HVec > 
                    >&rhs,
	      const Subset& s)
{
  const OLattice< SU3Mat >& u = static_cast< const OLattice< SU3Mat >& >(rhs.expression().left().child());
  const OLattice< FVec >& a = static_cast< const OLattice< FVec >& >(rhs.expression().right().child());

  if( s.hasOrderedRep() ) { 
    int totalSize = s.end() - s.start() +1;

    ordered_fused_spin_proj_user_arg arg = {u, a, d, s.start(), inlineSpinProjDir1Plus};

    dispatch_to_threads(totalSize, arg, ordered_fused_spin_proj_evaluate_function);

    ///////////////////
    // Original code
    ////////////////////
    /*
    for(int site = s.start(); site <= s.end(); site++) { 
     HVec tmp;
      inlineSpinProjDir1Plus( (REAL *)&(a.elem(site).elem(0).elem(0).real()),
			      (REAL *)&(tmp.elem(0).elem(0).real()),
			      1);
      
      
      _inline_mult_adj_su3_mat_vec(u.elem(site).elem(), tmp.elem(0), d.elem(site).elem(0));
      _inline_mult_adj_su3_mat_vec(u.elem(site).elem(), tmp.elem(1), d.elem(site).elem(1));   
      
      
      }*/
  }
  else { 

    const int *tab = s.siteTable().slice();

   int totalSize = s.numSiteTable();

   unordered_fused_spin_proj_user_arg arg = {u, a, d, tab, inlineSpinProjDir1Plus};

    dispatch_to_threads(totalSize, arg, unordered_fused_spin_proj_evaluate_function);

    ///////////////////
    // Original code
    ////////////////////
    /*
    for(int j=0; j < s.numSiteTable(); j++) { 
      int site = tab[j];
      
      
      HVec tmp;
      inlineSpinProjDir1Plus( (REAL *)&(a.elem(site).elem(0).elem(0).real()),
			      (REAL *)&(tmp.elem(0).elem(0).real()),
			      1);
      
      
      _inline_mult_adj_su3_mat_vec(u.elem(site).elem(), tmp.elem(0), d.elem(site).elem(0));
      _inline_mult_adj_su3_mat_vec(u.elem(site).elem(), tmp.elem(1), d.elem(site).elem(1));   
      
      }*/
  }

}

// HalfVec = adj(u)*SpinProjectDir1Minus(Vec);
template<>
inline
void evaluate(OLattice< HVec >& d,
              const OpAssign& op,
              const QDPExpr<
	              BinaryNode< 
	                FnAdjMultSprojDir1Minus,
	               
	                UnaryNode< OpIdentity, 
	                  Reference< QDPType< SU3Mat, OLattice< SU3Mat > > > >,
	               
	                  UnaryNode< OpIdentity,
                          Reference< QDPType< FVec,   OLattice< FVec > > > >
                      >,
	              OLattice< HVec > 
                    >&rhs,
	      const Subset& s)
{
  const OLattice< SU3Mat >& u = static_cast< const OLattice< SU3Mat >& >(rhs.expression().left().child());
  const OLattice< FVec >& a = static_cast< const OLattice< FVec >& >(rhs.expression().right().child());

  if( s.hasOrderedRep() ) { 
    int totalSize = s.end() - s.start() +1;

    ordered_fused_spin_proj_user_arg arg = {u, a, d, s.start(), inlineSpinProjDir1Minus};

    dispatch_to_threads(totalSize, arg, ordered_fused_spin_proj_evaluate_function);

    ///////////////////
    // Original code
    ////////////////////
    /*
    for(int site = s.start(); site <= s.end(); site++) { 
      HVec tmp;
      inlineSpinProjDir1Minus( (REAL *)&(a.elem(site).elem(0).elem(0).real()),
			       (REAL *)&(tmp.elem(0).elem(0).real()),
			       1);
      
      
      _inline_mult_adj_su3_mat_vec(u.elem(site).elem(), tmp.elem(0), d.elem(site).elem(0));
      _inline_mult_adj_su3_mat_vec(u.elem(site).elem(), tmp.elem(1), d.elem(site).elem(1));   
      
      }*/
  }
  else { 
    
    const int *tab = s.siteTable().slice();

    int totalSize = s.numSiteTable();

    unordered_fused_spin_proj_user_arg arg = {u, a, d, tab, inlineSpinProjDir1Minus};

    dispatch_to_threads(totalSize, arg, unordered_fused_spin_proj_evaluate_function);

    ///////////////////
    // Original code
    ////////////////////
    /*
    for(int j=0; j < s.numSiteTable(); j++) { 
      int site = tab[j];
      
      HVec tmp;
      inlineSpinProjDir1Minus( (REAL *)&(a.elem(site).elem(0).elem(0).real()),
			       (REAL *)&(tmp.elem(0).elem(0).real()),
			       1);
      
      
      _inline_mult_adj_su3_mat_vec(u.elem(site).elem(), tmp.elem(0), d.elem(site).elem(0));
      _inline_mult_adj_su3_mat_vec(u.elem(site).elem(), tmp.elem(1), d.elem(site).elem(1));   
      
      }*/
  }
}


// HalfVec = adj(u)*SpinProjectDir2Plus(Vec);
template<>
inline
void evaluate(OLattice< HVec >& d,
              const OpAssign& op,
              const QDPExpr<
	              BinaryNode< 
	                FnAdjMultSprojDir2Plus,
	               
	                UnaryNode< OpIdentity, 
	                  Reference< QDPType< SU3Mat, OLattice< SU3Mat > > > >,
	               
	                  UnaryNode< OpIdentity,
                          Reference< QDPType< FVec,   OLattice< FVec > > > >
                      >,
	              OLattice< HVec > 
                    >&rhs,
	      const Subset& s)
{
  const OLattice< SU3Mat >& u = static_cast< const OLattice< SU3Mat >& >(rhs.expression().left().child());
  const OLattice< FVec >& a = static_cast< const OLattice< FVec >& >(rhs.expression().right().child());


  if( s.hasOrderedRep() ) { 
    int totalSize = s.end() - s.start() +1;

    ordered_fused_spin_proj_user_arg arg = {u, a, d, s.start(), inlineSpinProjDir2Plus};

    dispatch_to_threads(totalSize, arg, ordered_fused_spin_proj_evaluate_function);

    ///////////////////
    // Original code
    ////////////////////
    /*
    for(int site = s.start(); site <= s.end(); site++) { 
      HVec tmp;
      inlineSpinProjDir2Plus( (REAL *)&(a.elem(site).elem(0).elem(0).real()),
			      (REAL *)&(tmp.elem(0).elem(0).real()),
			      1);
      
      
      _inline_mult_adj_su3_mat_vec(u.elem(site).elem(), tmp.elem(0), d.elem(site).elem(0));
      _inline_mult_adj_su3_mat_vec(u.elem(site).elem(), tmp.elem(1), d.elem(site).elem(1));   
      
      }*/
  }
  else { 

    const int *tab = s.siteTable().slice();

    int totalSize = s.numSiteTable();

    unordered_fused_spin_proj_user_arg arg = {u, a, d, tab, inlineSpinProjDir2Plus};

    dispatch_to_threads(totalSize, arg, unordered_fused_spin_proj_evaluate_function);

    ///////////////////
    // Original code
    ////////////////////
    /*
    for(int j=0; j < s.numSiteTable(); j++) { 
      int site = tab[j];
      
      HVec tmp;
      inlineSpinProjDir2Plus( (REAL *)&(a.elem(site).elem(0).elem(0).real()),
			      (REAL *)&(tmp.elem(0).elem(0).real()),
			      1);
      
      
    _inline_mult_adj_su3_mat_vec(u.elem(site).elem(), tmp.elem(0), d.elem(site).elem(0));
    _inline_mult_adj_su3_mat_vec(u.elem(site).elem(), tmp.elem(1), d.elem(site).elem(1));   
    
    }*/
  }

}

// HalfVec = adj(u)*SpinProjectDir2Minus(Vec);
template<>
inline
void evaluate(OLattice< HVec >& d,
              const OpAssign& op,
              const QDPExpr<
	              BinaryNode< 
	                FnAdjMultSprojDir2Minus,
	               
	                UnaryNode< OpIdentity, 
	                  Reference< QDPType< SU3Mat, OLattice< SU3Mat > > > >,
	               
	                  UnaryNode< OpIdentity,
                          Reference< QDPType< FVec,   OLattice< FVec > > > >
                      >,
	              OLattice< HVec > 
                    >&rhs,
	      const Subset& s)
{
  const OLattice< SU3Mat >& u = static_cast< const OLattice< SU3Mat >& >(rhs.expression().left().child());
  const OLattice< FVec >& a = static_cast< const OLattice< FVec >& >(rhs.expression().right().child());

  if( s.hasOrderedRep() ) { 
    int totalSize = s.end() - s.start() +1;

    ordered_fused_spin_proj_user_arg arg = {u, a, d, s.start(), inlineSpinProjDir2Minus};

    dispatch_to_threads(totalSize, arg, ordered_fused_spin_proj_evaluate_function);

    ///////////////////
    // Original code
    ////////////////////
    /*
    for(int site = s.start(); site <= s.end(); site++) { 

      HVec tmp;
      inlineSpinProjDir2Minus( (REAL *)&(a.elem(site).elem(0).elem(0).real()),
			       (REAL *)&(tmp.elem(0).elem(0).real()),
			       1);
      
      
      _inline_mult_adj_su3_mat_vec(u.elem(site).elem(), tmp.elem(0), d.elem(site).elem(0));
      _inline_mult_adj_su3_mat_vec(u.elem(site).elem(), tmp.elem(1), d.elem(site).elem(1));   
      
      }*/
  }
  else { 

    const int *tab = s.siteTable().slice();

    int totalSize = s.numSiteTable();

    unordered_fused_spin_proj_user_arg arg = {u, a, d, tab, inlineSpinProjDir2Minus};

    dispatch_to_threads(totalSize, arg, unordered_fused_spin_proj_evaluate_function);

    ///////////////////
    // Original code
    ////////////////////
    /*
    for(int j=0; j < s.numSiteTable(); j++) { 
      int site = tab[j];
      
      
      HVec tmp;
      inlineSpinProjDir2Minus( (REAL *)&(a.elem(site).elem(0).elem(0).real()),
			       (REAL *)&(tmp.elem(0).elem(0).real()),
			       1);
      
      
      _inline_mult_adj_su3_mat_vec(u.elem(site).elem(), tmp.elem(0), d.elem(site).elem(0));
      _inline_mult_adj_su3_mat_vec(u.elem(site).elem(), tmp.elem(1), d.elem(site).elem(1));   
      
      }*/
  }

}
// HalfVec = adj(u)*SpinProjectDir3Plus(Vec);
template<>
inline
void evaluate(OLattice< HVec >& d,
              const OpAssign& op,
              const QDPExpr<
	              BinaryNode< 
	                FnAdjMultSprojDir3Plus,
	               
	                UnaryNode< OpIdentity, 
	                  Reference< QDPType< SU3Mat, OLattice< SU3Mat > > > >,
	               
	                  UnaryNode< OpIdentity,
                          Reference< QDPType< FVec,   OLattice< FVec > > > >
                      >,
	              OLattice< HVec > 
                    >&rhs,
	      const Subset& s)
{
  const OLattice< SU3Mat >& u = static_cast< const OLattice< SU3Mat >& >(rhs.expression().left().child());
  const OLattice< FVec >& a = static_cast< const OLattice< FVec >& >(rhs.expression().right().child());

  

  if( s.hasOrderedRep() ) { 
    int totalSize = s.end() - s.start() +1;

    ordered_fused_spin_proj_user_arg arg = {u, a, d, s.start(), inlineSpinProjDir3Plus};

    dispatch_to_threads(totalSize, arg, ordered_fused_spin_proj_evaluate_function);

    ///////////////////
    // Original code
    ////////////////////
    /*
    for(int site = s.start(); site <= s.end(); site++) { 
      HVec tmp;
      inlineSpinProjDir3Plus( (REAL *)&(a.elem(site).elem(0).elem(0).real()),
			      (REAL *)&(tmp.elem(0).elem(0).real()),
			      1);
      
      
      _inline_mult_adj_su3_mat_vec(u.elem(site).elem(), tmp.elem(0), d.elem(site).elem(0));
      _inline_mult_adj_su3_mat_vec(u.elem(site).elem(), tmp.elem(1), d.elem(site).elem(1));   
      
      }*/
  }
  else { 

    const int *tab = s.siteTable().slice();

    int totalSize = s.numSiteTable();

    unordered_fused_spin_proj_user_arg arg = {u, a, d, tab, inlineSpinProjDir3Plus};

    dispatch_to_threads(totalSize, arg, unordered_fused_spin_proj_evaluate_function);

    ///////////////////
    // Original code
    ////////////////////
    /*
    for(int j=0; j < s.numSiteTable(); j++) { 
      int site = tab[j];
      
      HVec tmp;
      inlineSpinProjDir3Plus( (REAL *)&(a.elem(site).elem(0).elem(0).real()),
			      (REAL *)&(tmp.elem(0).elem(0).real()),
			      1);
      
      
      _inline_mult_adj_su3_mat_vec(u.elem(site).elem(), tmp.elem(0), d.elem(site).elem(0));
      _inline_mult_adj_su3_mat_vec(u.elem(site).elem(), tmp.elem(1), d.elem(site).elem(1));   
      
      }*/
  }
}

// HalfVec = adj(u)*SpinProjectDir3Minus(Vec);
template<>
inline
void evaluate(OLattice< HVec >& d,
              const OpAssign& op,
              const QDPExpr<
	              BinaryNode< 
	                FnAdjMultSprojDir3Minus,
	               
	                UnaryNode< OpIdentity, 
	                  Reference< QDPType< SU3Mat, OLattice< SU3Mat > > > >,
	               
	                  UnaryNode< OpIdentity,
                          Reference< QDPType< FVec,   OLattice< FVec > > > >
                      >,
	              OLattice< HVec > 
                    >&rhs,
	      const Subset& s)
{
  const OLattice< SU3Mat >& u = static_cast< const OLattice< SU3Mat >& >(rhs.expression().left().child());
  const OLattice< FVec >& a = static_cast< const OLattice< FVec >& >(rhs.expression().right().child());

  if( s.hasOrderedRep() ) { 
    
    int totalSize = s.end() - s.start() +1;

    ordered_fused_spin_proj_user_arg arg = {u, a, d, s.start(), inlineSpinProjDir3Minus};

    dispatch_to_threads(totalSize, arg, ordered_fused_spin_proj_evaluate_function);

    ///////////////////
    // Original code
    ////////////////////
    /*

    for(int site = s.start(); site <= s.end(); site++) { 
 
      HVec tmp;
      inlineSpinProjDir3Minus( (REAL *)&(a.elem(site).elem(0).elem(0).real()),
			       (REAL *)&(tmp.elem(0).elem(0).real()),
			       1);
      
      
      _inline_mult_adj_su3_mat_vec(u.elem(site).elem(), tmp.elem(0), d.elem(site).elem(0));
      _inline_mult_adj_su3_mat_vec(u.elem(site).elem(), tmp.elem(1), d.elem(site).elem(1));   
      
      }*/
  }
  else { 

    const int *tab = s.siteTable().slice();

    int totalSize = s.numSiteTable();

    unordered_fused_spin_proj_user_arg arg = {u, a, d, tab, inlineSpinProjDir3Minus};

    dispatch_to_threads(totalSize, arg, unordered_fused_spin_proj_evaluate_function);

    ///////////////////
    // Original code
    ////////////////////
    /*
    for(int j=0; j < s.numSiteTable(); j++) { 
      int site = tab[j];
      
      HVec tmp;
      inlineSpinProjDir3Minus( (REAL *)&(a.elem(site).elem(0).elem(0).real()),
			       (REAL *)&(tmp.elem(0).elem(0).real()),
			       1);
      
      
      _inline_mult_adj_su3_mat_vec(u.elem(site).elem(), tmp.elem(0), d.elem(site).elem(0));
      _inline_mult_adj_su3_mat_vec(u.elem(site).elem(), tmp.elem(1), d.elem(site).elem(1));   
      
      }*/
  }
}

} // namespace QDP;

#endif
