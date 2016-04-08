#ifndef QDP_SSE_FUSED_SPIN_PROJ_H
#define QDP_SSE_FUSED_SPIN_PROJ_H

#include "sse_mult_adj_su3_mat_hwvec.h"

/* Evaluates for things like adj(u)*spinProjectDir0Plus(y) */
using namespace QDP;
namespace QDP {

typedef PScalar< PColorMatrix< RComplex<REAL32>, 3> > SU3Mat32;


////////////////////////////////
// Threading evaluates
//
// by Xu Guo, EPCC, 20 October, 2008
////////////////////////////////

// the wrappers for the functions to be threaded
#include "qdp_sse_fused_spin_proj_evaluates_wrapper.h"


// HalfVec = adj(u)*SpinProjectDir0Plus(Vec);
template<>
inline
void evaluate(OLattice< HVec32 >& d,
              const OpAssign& op,
              const QDPExpr<
	              BinaryNode< 
	                FnAdjMultSprojDir0Plus,
	               
	                UnaryNode< OpIdentity, 
	                  Reference< QDPType< SU3Mat32, OLattice< SU3Mat32 > > > >,
	               
	                  UnaryNode< OpIdentity,
                          Reference< QDPType< FVec32,   OLattice< FVec32 > > > >
                      >,
	              OLattice< HVec32 > 
                    >&rhs,
	      const Subset& s)
{
  const OLattice< SU3Mat32 >& u = static_cast< const OLattice< SU3Mat32 >& >(rhs.expression().left().child());
  const OLattice< FVec32 >& a = static_cast< const OLattice< FVec32 >& >(rhs.expression().right().child());

  if( s.hasOrderedRep() ) {

    int totalSize = s.end() - s.start() +1;

    ordered_sse_fused_spin_proj_user_arg arg = {u, a, d, s.start(), inlineSpinProjDir0Plus};

    dispatch_to_threads(totalSize, arg, ordered_sse_fused_spin_proj_evaluate_function);

    ///////////////////
    // Original code
    ////////////////////
    /*
    for(int site = s.start() ; site <= s.end(); site++) { 
      HVec32 tmp ; 
      inlineSpinProjDir0Plus( (REAL32 *)&(a.elem(site).elem(0).elem(0).real()),
			      (REAL32 *)&(tmp.elem(0).elem(0).real()),
			      1);
      
      su3_matrixf* um = (su3_matrixf *)&(u.elem(site).elem().elem(0,0).real());
      half_wilson_vectorf *tmph = (half_wilson_vectorf *)&( tmp.elem(0).elem(0).real());
      half_wilson_vectorf *dh = (half_wilson_vectorf *)&( d.elem(site).elem(0).elem(0).real());
      
      intrin_sse_mult_adj_su3_mat_hwvec(um, tmph, dh);
      }*/
  }
  else { 
    const int* tab = s.siteTable().slice();

    int totalSize = s.numSiteTable();

    unordered_sse_fused_spin_proj_user_arg arg = {u, a, d, tab, inlineSpinProjDir0Plus};

    dispatch_to_threads(totalSize, arg, unordered_sse_fused_spin_proj_evaluate_function);

    ///////////////////
    // Original code
    ////////////////////
    /*
    for(int j=0; j < s.numSiteTable(); j++) { 
      int site=tab[j];
      HVec32 tmp ; 
      inlineSpinProjDir0Plus( (REAL32 *)&(a.elem(site).elem(0).elem(0).real()),
			      (REAL32 *)&(tmp.elem(0).elem(0).real()),
			      1);
      
      su3_matrixf* um = (su3_matrixf *)&(u.elem(site).elem().elem(0,0).real());
      half_wilson_vectorf *tmph = (half_wilson_vectorf *)&( tmp.elem(0).elem(0).real());
      half_wilson_vectorf *dh = (half_wilson_vectorf *)&( d.elem(site).elem(0).elem(0).real());
      
      intrin_sse_mult_adj_su3_mat_hwvec(um, tmph, dh);
    
      }*/
  }

}

// HalfVec = adj(u)*SpinProjectDir0Minus(Vec);
template<>
inline
void evaluate(OLattice< HVec32 >& d,
              const OpAssign& op,
              const QDPExpr<
	              BinaryNode< 
	                FnAdjMultSprojDir0Minus,
	               
	                UnaryNode< OpIdentity, 
	                  Reference< QDPType< SU3Mat32, OLattice< SU3Mat32 > > > >,
	               
	                  UnaryNode< OpIdentity,
                          Reference< QDPType< FVec32,   OLattice< FVec32 > > > >
                      >,
	              OLattice< HVec32 > 
                    >&rhs,
	      const Subset& s)
{
  const OLattice< SU3Mat32 >& u = static_cast< const OLattice< SU3Mat32 >& >(rhs.expression().left().child());
  const OLattice< FVec32 >& a = static_cast< const OLattice< FVec32 >& >(rhs.expression().right().child());

  if( s.hasOrderedRep() ) {

    int totalSize = s.end() - s.start() +1;

    ordered_sse_fused_spin_proj_user_arg arg = {u, a, d, s.start(), inlineSpinProjDir0Minus};

    dispatch_to_threads(totalSize, arg, ordered_sse_fused_spin_proj_evaluate_function);

    ///////////////////
    // Original code
    ////////////////////
    /*
    for(int site = s.start() ; site <= s.end(); site++) { 
      HVec32 tmp ;
      inlineSpinProjDir0Minus( (REAL32 *)&(a.elem(site).elem(0).elem(0).real()),
			       (REAL32 *)&(tmp.elem(0).elem(0).real()),
			       1);
      
      su3_matrixf* um = (su3_matrixf *)&(u.elem(site).elem().elem(0,0).real());
      half_wilson_vectorf *tmph = (half_wilson_vectorf *)&( tmp.elem(0).elem(0).real());
      half_wilson_vectorf *dh = (half_wilson_vectorf *)&( d.elem(site).elem(0).elem(0).real());
      intrin_sse_mult_adj_su3_mat_hwvec(um, tmph, dh);
   
      }*/
  }
  else { 
    const int* tab = s.siteTable().slice();

    int totalSize = s.numSiteTable();

    unordered_sse_fused_spin_proj_user_arg arg = {u, a, d, tab, inlineSpinProjDir0Minus};

    dispatch_to_threads(totalSize, arg, unordered_sse_fused_spin_proj_evaluate_function);

    ///////////////////
    // Original code
    ////////////////////
    /*
    for(int j=0; j < s.numSiteTable(); j++) { 
      int site=tab[j];
      
      HVec32 tmp ;
      inlineSpinProjDir0Minus( (REAL32 *)&(a.elem(site).elem(0).elem(0).real()),
			       (REAL32 *)&(tmp.elem(0).elem(0).real()),
			       1);
      
      su3_matrixf* um = (su3_matrixf *)&(u.elem(site).elem().elem(0,0).real());
      half_wilson_vectorf *tmph = (half_wilson_vectorf *)&( tmp.elem(0).elem(0).real());
      half_wilson_vectorf *dh = (half_wilson_vectorf *)&( d.elem(site).elem(0).elem(0).real());
      intrin_sse_mult_adj_su3_mat_hwvec(um, tmph, dh);
      
      }*/
  }

}

// HalfVec = adj(u)*SpinProjectDir1Plus(Vec);
template<>
inline
void evaluate(OLattice< HVec32 >& d,
              const OpAssign& op,
              const QDPExpr<
	              BinaryNode< 
	                FnAdjMultSprojDir1Plus,
	               
	                UnaryNode< OpIdentity, 
	                  Reference< QDPType< SU3Mat32, OLattice< SU3Mat32 > > > >,
	               
	                  UnaryNode< OpIdentity,
                          Reference< QDPType< FVec32,   OLattice< FVec32 > > > >
                      >,
	              OLattice< HVec32 > 
                    >&rhs,
	      const Subset& s)
{
  const OLattice< SU3Mat32 >& u = static_cast< const OLattice< SU3Mat32 >& >(rhs.expression().left().child());
  const OLattice< FVec32 >& a = static_cast< const OLattice< FVec32 >& >(rhs.expression().right().child());

  
  if( s.hasOrderedRep() ) {

    int totalSize = s.end() - s.start() +1;

    ordered_sse_fused_spin_proj_user_arg arg = {u, a, d, s.start(), inlineSpinProjDir1Plus};

    dispatch_to_threads(totalSize, arg, ordered_sse_fused_spin_proj_evaluate_function);

    ///////////////////
    // Original code
    ////////////////////
    /*
    for(int site = s.start() ; site <= s.end(); site++) { 
      HVec32 tmp ;
      inlineSpinProjDir1Plus( (REAL32 *)&(a.elem(site).elem(0).elem(0).real()),
			      (REAL32 *)&(tmp.elem(0).elem(0).real()),
			      1);
      
      su3_matrixf* um = (su3_matrixf *)&(u.elem(site).elem().elem(0,0).real());
      half_wilson_vectorf *tmph = (half_wilson_vectorf *)&( tmp.elem(0).elem(0).real());
      half_wilson_vectorf *dh = (half_wilson_vectorf *)&( d.elem(site).elem(0).elem(0).real());
      intrin_sse_mult_adj_su3_mat_hwvec(um, tmph, dh);
      
      
      }*/
  }
  else { 
    
    const int* tab = s.siteTable().slice();

    int totalSize = s.numSiteTable();

    unordered_sse_fused_spin_proj_user_arg arg = {u, a, d, tab, inlineSpinProjDir1Plus};

    dispatch_to_threads(totalSize, arg, unordered_sse_fused_spin_proj_evaluate_function);

    ///////////////////
    // Original code
    ////////////////////
    /*
    for(int j=0; j < s.numSiteTable(); j++) { 
      int site=tab[j];
      
      HVec32 tmp ;
      inlineSpinProjDir1Plus( (REAL32 *)&(a.elem(site).elem(0).elem(0).real()),
			      (REAL32 *)&(tmp.elem(0).elem(0).real()),
			      1);
      
      su3_matrixf* um = (su3_matrixf *)&(u.elem(site).elem().elem(0,0).real());
      half_wilson_vectorf *tmph = (half_wilson_vectorf *)&( tmp.elem(0).elem(0).real());
      half_wilson_vectorf *dh = (half_wilson_vectorf *)&( d.elem(site).elem(0).elem(0).real());
      intrin_sse_mult_adj_su3_mat_hwvec(um, tmph, dh);
      
      }*/
  }
}

// HalfVec = adj(u)*SpinProjectDir1Minus(Vec);
template<>
inline
void evaluate(OLattice< HVec32 >& d,
              const OpAssign& op,
              const QDPExpr<
	              BinaryNode< 
	                FnAdjMultSprojDir1Minus,
	               
	                UnaryNode< OpIdentity, 
	                  Reference< QDPType< SU3Mat32, OLattice< SU3Mat32 > > > >,
	               
	                  UnaryNode< OpIdentity,
                          Reference< QDPType< FVec32,   OLattice< FVec32 > > > >
                      >,
	              OLattice< HVec32 > 
                    >&rhs,
	      const Subset& s)
{
  const OLattice< SU3Mat32 >& u = static_cast< const OLattice< SU3Mat32 >& >(rhs.expression().left().child());
  const OLattice< FVec32 >& a = static_cast< const OLattice< FVec32 >& >(rhs.expression().right().child());

  if( s.hasOrderedRep() ) {

    int totalSize = s.end() - s.start() +1;

    ordered_sse_fused_spin_proj_user_arg arg = {u, a, d, s.start(), inlineSpinProjDir1Minus};

    dispatch_to_threads(totalSize, arg, ordered_sse_fused_spin_proj_evaluate_function);

    ///////////////////
    // Original code
    ////////////////////
    /*
    for(int site = s.start() ; site <= s.end(); site++) { 
      HVec32 tmp ;
      inlineSpinProjDir1Minus( (REAL32 *)&(a.elem(site).elem(0).elem(0).real()),
			       (REAL32 *)&(tmp.elem(0).elem(0).real()),
			       1);
      
      su3_matrixf* um = (su3_matrixf *)&(u.elem(site).elem().elem(0,0).real());
      half_wilson_vectorf *tmph = (half_wilson_vectorf *)&( tmp.elem(0).elem(0).real());
      half_wilson_vectorf *dh = (half_wilson_vectorf *)&( d.elem(site).elem(0).elem(0).real());
      intrin_sse_mult_adj_su3_mat_hwvec(um, tmph, dh);
   
      }*/
  }
  else { 

    const int* tab = s.siteTable().slice();

    int totalSize = s.numSiteTable();

    unordered_sse_fused_spin_proj_user_arg arg = {u, a, d, tab, inlineSpinProjDir1Minus};

    dispatch_to_threads(totalSize, arg, unordered_sse_fused_spin_proj_evaluate_function);

    ///////////////////
    // Original code
    ////////////////////
    /*
    for(int j=0; j < s.numSiteTable(); j++) { 
      int site=tab[j];
      
      HVec32 tmp ;
      inlineSpinProjDir1Minus( (REAL32 *)&(a.elem(site).elem(0).elem(0).real()),
			       (REAL32 *)&(tmp.elem(0).elem(0).real()),
			       1);
      
      su3_matrixf* um = (su3_matrixf *)&(u.elem(site).elem().elem(0,0).real());
      half_wilson_vectorf *tmph = (half_wilson_vectorf *)&( tmp.elem(0).elem(0).real());
      half_wilson_vectorf *dh = (half_wilson_vectorf *)&( d.elem(site).elem(0).elem(0).real());
      intrin_sse_mult_adj_su3_mat_hwvec(um, tmph, dh);
      
      }*/
  }
}


// HalfVec = adj(u)*SpinProjectDir2Plus(Vec);
template<>
inline
void evaluate(OLattice< HVec32 >& d,
              const OpAssign& op,
              const QDPExpr<
	              BinaryNode< 
	                FnAdjMultSprojDir2Plus,
	               
	                UnaryNode< OpIdentity, 
	                  Reference< QDPType< SU3Mat32, OLattice< SU3Mat32 > > > >,
	               
	                  UnaryNode< OpIdentity,
                          Reference< QDPType< FVec32,   OLattice< FVec32 > > > >
                      >,
	              OLattice< HVec32 > 
                    >&rhs,
	      const Subset& s)
{
  const OLattice< SU3Mat32 >& u = static_cast< const OLattice< SU3Mat32 >& >(rhs.expression().left().child());
  const OLattice< FVec32 >& a = static_cast< const OLattice< FVec32 >& >(rhs.expression().right().child());

  if( s.hasOrderedRep() ) {

    int totalSize = s.end() - s.start() +1;

    ordered_sse_fused_spin_proj_user_arg arg = {u, a, d, s.start(), inlineSpinProjDir2Plus};

    dispatch_to_threads(totalSize, arg, ordered_sse_fused_spin_proj_evaluate_function);

    ///////////////////
    // Original code
    ////////////////////
    /*
    for(int site = s.start() ; site <= s.end(); site++) { 
      HVec32 tmp ;
      inlineSpinProjDir2Plus( (REAL32 *)&(a.elem(site).elem(0).elem(0).real()),
			      (REAL32 *)&(tmp.elem(0).elem(0).real()),
			      1);
      su3_matrixf* um = (su3_matrixf *)&(u.elem(site).elem().elem(0,0).real());
      half_wilson_vectorf *tmph = (half_wilson_vectorf *)&( tmp.elem(0).elem(0).real());
      half_wilson_vectorf *dh = (half_wilson_vectorf *)&( d.elem(site).elem(0).elem(0).real());
      intrin_sse_mult_adj_su3_mat_hwvec(um, tmph, dh);
      
      }*/
  }
  else { 
    const int* tab = s.siteTable().slice();

    int totalSize = s.numSiteTable();

    unordered_sse_fused_spin_proj_user_arg arg = {u, a, d, tab, inlineSpinProjDir2Plus};

    dispatch_to_threads(totalSize, arg, unordered_sse_fused_spin_proj_evaluate_function);

    ///////////////////
    // Original code
    ////////////////////
    /*
    for(int j=0; j < s.numSiteTable(); j++) { 
      int site=tab[j];
      HVec32 tmp ;
      inlineSpinProjDir2Plus( (REAL32 *)&(a.elem(site).elem(0).elem(0).real()),
			      (REAL32 *)&(tmp.elem(0).elem(0).real()),
			      1);

      su3_matrixf* um = (su3_matrixf *)&(u.elem(site).elem().elem(0,0).real());
      half_wilson_vectorf *tmph = (half_wilson_vectorf *)&( tmp.elem(0).elem(0).real());
      half_wilson_vectorf *dh = (half_wilson_vectorf *)&( d.elem(site).elem(0).elem(0).real());
      intrin_sse_mult_adj_su3_mat_hwvec(um, tmph, dh);
      
      }*/
  }
}

// HalfVec = adj(u)*SpinProjectDir2Minus(Vec);
template<>
inline
void evaluate(OLattice< HVec32 >& d,
              const OpAssign& op,
              const QDPExpr<
	              BinaryNode< 
	                FnAdjMultSprojDir2Minus,
	               
	                UnaryNode< OpIdentity, 
	                  Reference< QDPType< SU3Mat32, OLattice< SU3Mat32 > > > >,
	               
	                  UnaryNode< OpIdentity,
                          Reference< QDPType< FVec32,   OLattice< FVec32 > > > >
                      >,
	              OLattice< HVec32 > 
                    >&rhs,
	      const Subset& s)
{
  const OLattice< SU3Mat32 >& u = static_cast< const OLattice< SU3Mat32 >& >(rhs.expression().left().child());
  const OLattice< FVec32 >& a = static_cast< const OLattice< FVec32 >& >(rhs.expression().right().child());

  if( s.hasOrderedRep() ) {

    int totalSize = s.end() - s.start() +1;

    ordered_sse_fused_spin_proj_user_arg arg = {u, a, d, s.start(), inlineSpinProjDir2Minus};

    dispatch_to_threads(totalSize, arg, ordered_sse_fused_spin_proj_evaluate_function);

    ///////////////////
    // Original code
    ////////////////////
    /*
    for(int site = s.start() ; site <= s.end(); site++) { 
      HVec32 tmp ;
      inlineSpinProjDir2Minus( (REAL32 *)&(a.elem(site).elem(0).elem(0).real()),
			       (REAL32 *)&(tmp.elem(0).elem(0).real()),
			       1);
      su3_matrixf* um = (su3_matrixf *)&(u.elem(site).elem().elem(0,0).real());
      half_wilson_vectorf *tmph = (half_wilson_vectorf *)&( tmp.elem(0).elem(0).real());
      half_wilson_vectorf *dh = (half_wilson_vectorf *)&( d.elem(site).elem(0).elem(0).real());
      intrin_sse_mult_adj_su3_mat_hwvec(um, tmph, dh);
      
      
      }*/
  }
  else { 
    
    const int* tab = s.siteTable().slice();

    int totalSize = s.numSiteTable();

    unordered_sse_fused_spin_proj_user_arg arg = {u, a, d, tab, inlineSpinProjDir2Minus};

    dispatch_to_threads(totalSize, arg, unordered_sse_fused_spin_proj_evaluate_function);

    ///////////////////
    // Original code
    ////////////////////
    /*
    for(int j=0; j < s.numSiteTable(); j++) { 
      int site=tab[j];
      
      HVec32 tmp ;
      inlineSpinProjDir2Minus( (REAL32 *)&(a.elem(site).elem(0).elem(0).real()),
			       (REAL32 *)&(tmp.elem(0).elem(0).real()),
			       1);
      
      
      su3_matrixf* um = (su3_matrixf *)&(u.elem(site).elem().elem(0,0).real());
      half_wilson_vectorf *tmph = (half_wilson_vectorf *)&( tmp.elem(0).elem(0).real());
      half_wilson_vectorf *dh = (half_wilson_vectorf *)&( d.elem(site).elem(0).elem(0).real());
      intrin_sse_mult_adj_su3_mat_hwvec(um, tmph, dh);
      
      
      }*/
  }
}

// HalfVec = adj(u)*SpinProjectDir3Plus(Vec);
template<>
inline
void evaluate(OLattice< HVec32 >& d,
              const OpAssign& op,
              const QDPExpr<
	              BinaryNode< 
	                FnAdjMultSprojDir3Plus,
	               
	                UnaryNode< OpIdentity, 
	                  Reference< QDPType< SU3Mat32, OLattice< SU3Mat32 > > > >,
	               
	                  UnaryNode< OpIdentity,
                          Reference< QDPType< FVec32,   OLattice< FVec32 > > > >
                      >,
	              OLattice< HVec32 > 
                    >&rhs,
	      const Subset& s)
{
  const OLattice< SU3Mat32 >& u = static_cast< const OLattice< SU3Mat32 >& >(rhs.expression().left().child());
  const OLattice< FVec32 >& a = static_cast< const OLattice< FVec32 >& >(rhs.expression().right().child());

  if( s.hasOrderedRep() ) {

    int totalSize = s.end() - s.start() +1;

    ordered_sse_fused_spin_proj_user_arg arg = {u, a, d, s.start(), inlineSpinProjDir3Plus};

    dispatch_to_threads(totalSize, arg, ordered_sse_fused_spin_proj_evaluate_function);

    ///////////////////
    // Original code
    ////////////////////
    /*
    for(int site = s.start() ; site <= s.end(); site++) { 
      HVec32 tmp ;
      inlineSpinProjDir3Plus( (REAL32 *)&(a.elem(site).elem(0).elem(0).real()),
			      (REAL32 *)&(tmp.elem(0).elem(0).real()),
			      1);
      su3_matrixf* um = (su3_matrixf *)&(u.elem(site).elem().elem(0,0).real());
      half_wilson_vectorf *tmph = (half_wilson_vectorf *)&( tmp.elem(0).elem(0).real());
      half_wilson_vectorf *dh = (half_wilson_vectorf *)&( d.elem(site).elem(0).elem(0).real());
      intrin_sse_mult_adj_su3_mat_hwvec(um, tmph, dh);
      
      
      }*/
  }
  else { 

    const int* tab = s.siteTable().slice();

    int totalSize = s.numSiteTable();

    unordered_sse_fused_spin_proj_user_arg arg = {u, a, d, tab, inlineSpinProjDir3Plus};

    dispatch_to_threads(totalSize, arg, unordered_sse_fused_spin_proj_evaluate_function);

    ///////////////////
    // Original code
    ////////////////////
    /*
    for(int j=0; j < s.numSiteTable(); j++) { 
      int site=tab[j];
      HVec32 tmp ;
      inlineSpinProjDir3Plus( (REAL32 *)&(a.elem(site).elem(0).elem(0).real()),
			      (REAL32 *)&(tmp.elem(0).elem(0).real()),
			      1);
      
      su3_matrixf* um = (su3_matrixf *)&(u.elem(site).elem().elem(0,0).real());
      half_wilson_vectorf *tmph = (half_wilson_vectorf *)&( tmp.elem(0).elem(0).real());
      half_wilson_vectorf *dh = (half_wilson_vectorf *)&( d.elem(site).elem(0).elem(0).real());
      intrin_sse_mult_adj_su3_mat_hwvec(um, tmph, dh);
      
      
    }*/
  }
}

// HalfVec = adj(u)*SpinProjectDir3Minus(Vec);
template<>
inline
void evaluate(OLattice< HVec32 >& d,
              const OpAssign& op,
              const QDPExpr<
	              BinaryNode< 
	                FnAdjMultSprojDir3Minus,
	               
	                UnaryNode< OpIdentity, 
	                  Reference< QDPType< SU3Mat32, OLattice< SU3Mat32 > > > >,
	               
	                  UnaryNode< OpIdentity,
                          Reference< QDPType< FVec32,   OLattice< FVec32 > > > >
                      >,
	              OLattice< HVec32 > 
                    >&rhs,
	      const Subset& s)
{
  const OLattice< SU3Mat32 >& u = static_cast< const OLattice< SU3Mat32 >& >(rhs.expression().left().child());
  const OLattice< FVec32 >& a = static_cast< const OLattice< FVec32 >& >(rhs.expression().right().child());


  if( s.hasOrderedRep() ) {

    int totalSize = s.end() - s.start() +1;

    ordered_sse_fused_spin_proj_user_arg arg = {u, a, d, s.start(), inlineSpinProjDir3Minus};

    dispatch_to_threads(totalSize, arg, ordered_sse_fused_spin_proj_evaluate_function);

    ///////////////////
    // Original code
    ////////////////////
    /*
    for(int site = s.start() ; site <= s.end(); site++) { 
      HVec32 tmp ;
      inlineSpinProjDir3Minus( (REAL32 *)&(a.elem(site).elem(0).elem(0).real()),
			       (REAL32 *)&(tmp.elem(0).elem(0).real()),
			       1);
      
      su3_matrixf* um = (su3_matrixf *)&(u.elem(site).elem().elem(0,0).real());
      half_wilson_vectorf *tmph = (half_wilson_vectorf *)&( tmp.elem(0).elem(0).real());
      half_wilson_vectorf *dh = (half_wilson_vectorf *)&( d.elem(site).elem(0).elem(0).real());
      intrin_sse_mult_adj_su3_mat_hwvec(um, tmph, dh);
      

      }*/
  }
  else { 
    
    const int* tab = s.siteTable().slice();

    int totalSize = s.numSiteTable();

    unordered_sse_fused_spin_proj_user_arg arg = {u, a, d, tab, inlineSpinProjDir3Minus};

    dispatch_to_threads(totalSize, arg, unordered_sse_fused_spin_proj_evaluate_function);

    ///////////////////
    // Original code
    ////////////////////
    /*
    for(int j=0; j < s.numSiteTable(); j++) { 
      int site=tab[j];
      
      HVec32 tmp ;
      inlineSpinProjDir3Minus( (REAL32 *)&(a.elem(site).elem(0).elem(0).real()),
			       (REAL32 *)&(tmp.elem(0).elem(0).real()),
			       1);
      
      su3_matrixf* um = (su3_matrixf *)&(u.elem(site).elem().elem(0,0).real());
      half_wilson_vectorf *tmph = (half_wilson_vectorf *)&( tmp.elem(0).elem(0).real());
      half_wilson_vectorf *dh = (half_wilson_vectorf *)&( d.elem(site).elem(0).elem(0).real());
      intrin_sse_mult_adj_su3_mat_hwvec(um, tmph, dh);
      
      }*/
  }

}

} // namespace QDP;

#endif
