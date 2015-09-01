// $Id: qdp_scalarsite_linalg_double.cc,v 1.4 2009-07-14 20:08:42 bjoo Exp $

/*! @file
 * @brief Intel SSE optimizations
 * 
 * SSE optimizations of basic operations
 */


#include "qdp.h"


// These SSE asm instructions are only supported under GCC/G++
#include "scalarsite_sse/qdp_scalarsite_sse_linalg_double.h"

namespace QDP {



//-------------------------------------------------------------------
// Specialization to optimize the case   
//    LatticeColorMatrix[ Subset] = LatticeColorMatrix * LatticeColorMatrix

// Threading hoist
  struct ordered_ssed_3mat_args {
    REAL64* lm;
    REAL64* rm;
    REAL64* dm;
    void (*func)(REAL64*, REAL64*, REAL64*, int);
  };

  inline 
  void ordered_ssed_3mat_evaluate_func(int lo, int hi, int myId, 
				       ordered_ssed_3mat_args *a)
  {
    void (*func)(REAL64*, REAL64*, REAL64*, int) = a->func;
    
    int n_3mat = (hi - lo);
    int index = lo*18; // 3*3*2 doubles
    REAL64* lm = &(a->lm[index]);
    REAL64* rm = &(a->rm[index]);
    REAL64* dm = &(a->dm[index]);

    func(dm,lm,rm, n_3mat);
  }

  struct unordered_ssed_3mat_args { 
    unordered_ssed_3mat_args(    const OLattice<DCol>& l_,
				 const OLattice<DCol>& r_,
				 OLattice<DCol>& d_,
				 const int *tab_,
				 void (*func_)(REAL64*, REAL64*, REAL64*, int) ) :
      l(l_),r(r_),d(d_),tab(tab_),func(func_) {}

    const OLattice<DCol>& l;
    const OLattice<DCol>& r;
    OLattice<DCol>& d;
    const int *tab;
    void (*func)(REAL64*, REAL64*, REAL64*, int);
  };

  inline 
  void unordered_ssed_3mat_evaluate_func(int lo, int hi, int myId,
					     unordered_ssed_3mat_args* a)
  {
    void (*func)(REAL64*, REAL64*, REAL64*, int) = a->func;
    for(int site=lo; site < hi; site++) { 
      int i = a->tab[site];
      REAL64 *lm = (REAL64 *)&(a->l.elem(i).elem().elem(0,0).real());
      REAL64 *rm = (REAL64 *)&(a->r.elem(i).elem().elem(0,0).real());
      REAL64 *dm = (REAL64 *)&(a->d.elem(i).elem().elem(0,0).real());
      func(dm,lm, rm, 1);
    }
  }
      


template<>
void evaluate(OLattice< DCol >& d, 
	      const OpAssign& op, 
	      const QDPExpr<BinaryNode<OpMultiply, 
	                    Reference<QDPType< DCol, OLattice< DCol > > >, 
	                    Reference<QDPType< DCol, OLattice< DCol > > > >,
	                    OLattice< DCol > >& rhs,
	      const Subset& s)
{
//  cout << "call single site QDP_M_eq_M_times_M" << endl;

  typedef OLattice< DCol >    C;

  const C& l = static_cast<const C&>(rhs.expression().left());
  const C& r = static_cast<const C&>(rhs.expression().right());

  if( s.hasOrderedRep() ) { 
    int n_3mat=s.end() - s.start() + 1;
    REAL64 *lm = (REAL64 *)&(l.elem(s.start()).elem().elem(0,0).real());
    REAL64 *rm = (REAL64 *)&(r.elem(s.start()).elem().elem(0,0).real());
    REAL64 *dm = (REAL64 *)&(d.elem(s.start()).elem().elem(0,0).real());
    
    ordered_ssed_3mat_args a = {lm,rm,dm,ssed_m_eq_mm};
    dispatch_to_threads(n_3mat, a, ordered_ssed_3mat_evaluate_func);

    //ssed_m_eq_mm(dm,lm,rm, n_3mat);
  }
  else { 
    const int *tab = s.siteTable().slice();
    unordered_ssed_3mat_args a(l,r,d,tab,ssed_m_eq_mm);
    dispatch_to_threads(s.numSiteTable(), a, unordered_ssed_3mat_evaluate_func);

    //for(int j=0; j < s.numSiteTable(); ++j) {
    // int i = tab[j];
    // REAL64 *lm = (REAL64 *)&(l.elem(i).elem().elem(0,0).real());
    // REAL64 *rm = (REAL64 *)&(r.elem(i).elem().elem(0,0).real());
    // REAL64 *dm = (REAL64 *)&(d.elem(i).elem().elem(0,0).real());
    // ssed_m_eq_mm(dm,lm, rm, 1);
    //    }
  }
}


// Specialization to optimize the case   
//    LatticeColorMatrix[Subset] = adj(LatticeColorMatrix) * LatticeColorMatrix
template<>
void evaluate(OLattice< DCol >& d, 
	      const OpAssign& op, 
	      const QDPExpr<BinaryNode<OpAdjMultiply, 
	                    UnaryNode<OpIdentity, Reference<QDPType< DCol, OLattice< DCol > > > >, 
	                    Reference<QDPType< DCol, OLattice< DCol > > > >,
	                    OLattice< DCol > >& rhs,
	      const Subset& s)
{
//  cout << "call single site QDP_M_eq_aM_times_M" << endl;

  typedef OLattice< DCol >    C;

  const C& l = static_cast<const C&>(rhs.expression().left().child());
  const C& r = static_cast<const C&>(rhs.expression().right());

  if( s.hasOrderedRep() ) { 
    int n_3mat = s.end()-s.start()+1;
    REAL64 *lm = (REAL64 *)&(l.elem(s.start()).elem().elem(0,0).real());
    REAL64 *rm = (REAL64 *)&(r.elem(s.start()).elem().elem(0,0).real());
    REAL64 *dm = (REAL64 *)&(d.elem(s.start()).elem().elem(0,0).real());
    ordered_ssed_3mat_args a = {lm,rm,dm,ssed_m_eq_hm};
    dispatch_to_threads(n_3mat, a, ordered_ssed_3mat_evaluate_func);
    
    //ssed_m_eq_hm(dm,lm,rm,n_3mat);
  }
  else { 
    const int *tab = s.siteTable().slice();
    unordered_ssed_3mat_args a(l,r,d,tab,ssed_m_eq_hm);
    dispatch_to_threads(s.numSiteTable(), a, unordered_ssed_3mat_evaluate_func);

    //for(int j=0; j < s.numSiteTable(); ++j) {
    //
    //  int i = tab[j];
    //  REAL64 *lm = (REAL64 *)&(l.elem(i).elem().elem(0,0).real());
    //  REAL64 *rm = (REAL64 *)&(r.elem(i).elem().elem(0,0).real());
    //  REAL64 *dm = (REAL64 *)&(d.elem(i).elem().elem(0,0).real());
    //  ssed_m_eq_hm(dm, lm, rm, 1);
    //    }
  }

}


// Specialization to optimize the case   
//    LatticeColorMatrix[Subset] = LatticeColorMatrix * adj(LatticeColorMatrix)
template<>
void evaluate(OLattice< DCol >& d, 
	      const OpAssign& op, 
	      const QDPExpr<BinaryNode<OpMultiplyAdj, 
	                    Reference<QDPType< DCol, OLattice< DCol > > >, 
	                    UnaryNode<OpIdentity, Reference<QDPType< DCol, OLattice< DCol > > > > >,
	                    OLattice< DCol > >& rhs,
	      const Subset& s)
{
//  cout << "call single site QDP_M_eq_M_times_aM" << endl;

  typedef OLattice< DCol >    C;

  const C& l = static_cast<const C&>(rhs.expression().left());
  const C& r = static_cast<const C&>(rhs.expression().right().child());

  if( s.hasOrderedRep() ) { 
      REAL64 *lm = (REAL64 *)&(l.elem(s.start()).elem().elem(0,0).real());
      REAL64 *rm = (REAL64 *)&(r.elem(s.start()).elem().elem(0,0).real());
      REAL64 *dm = (REAL64 *)&(d.elem(s.start()).elem().elem(0,0).real());
      int n_3mat = s.end()-s.start()+1;
      ordered_ssed_3mat_args a = {lm,rm,dm,ssed_m_eq_mh};
      dispatch_to_threads(n_3mat, a, ordered_ssed_3mat_evaluate_func);

      //ssed_m_eq_mh(dm,lm,rm,n_3mat);
  }
  else { 

    const int *tab = s.siteTable().slice();
    unordered_ssed_3mat_args a(l,r,d,tab,ssed_m_eq_mh);
    dispatch_to_threads(s.numSiteTable(), a, unordered_ssed_3mat_evaluate_func);

    // for(int j=0; j < s.numSiteTable(); ++j) {
    //  int i = tab[j];
    //  REAL64 *lm = (REAL64 *)&(l.elem(i).elem().elem(0,0).real());
    //  REAL64 *rm = (REAL64 *)&(r.elem(i).elem().elem(0,0).real());
    //  REAL64 *dm = (REAL64 *)&(d.elem(i).elem().elem(0,0).real());
    //  ssed_m_eq_mh(dm, lm, rm, 1);
    //
    // }
  }
}


// Specialization to optimize the case   
//    LatticeColorMatrix[Subset] = adj(LatticeColorMatrix) * adj(LatticeColorMatrix)
template<>
void evaluate(OLattice< DCol >& d, 
	      const OpAssign& op, 
	      const QDPExpr<BinaryNode<OpAdjMultiplyAdj, 
	                    UnaryNode<OpIdentity, Reference<QDPType< DCol, OLattice< DCol > > > >,
	                    UnaryNode<OpIdentity, Reference<QDPType< DCol, OLattice< DCol > > > > >,
	                    OLattice< DCol > >& rhs,
	      const Subset& s)
{
//  cout << "call single site QDP_M_eq_Ma_times_Ma" << endl;

  typedef OLattice< DCol >    C;

  const C& l = static_cast<const C&>(rhs.expression().left().child());
  const C& r = static_cast<const C&>(rhs.expression().right().child());

  if( s.hasOrderedRep() ) { 
    int n_3mat = s.end() - s.start() + 1;
    REAL64 *lm = (REAL64 *)&(l.elem(s.start()).elem().elem(0,0).real());
    REAL64 *rm = (REAL64 *)&(r.elem(s.start()).elem().elem(0,0).real());
    REAL64 *dm = (REAL64 *)&(d.elem(s.start()).elem().elem(0,0).real());
    ordered_ssed_3mat_args a = {lm,rm,dm,ssed_m_eq_hh};
    dispatch_to_threads(n_3mat, a, ordered_ssed_3mat_evaluate_func);
    
    //ssed_m_eq_hh(dm,lm,rm,n_3mat);
  }
  else { 
    const int *tab = s.siteTable().slice();
    unordered_ssed_3mat_args a(l,r,d,tab,ssed_m_eq_hh);
    dispatch_to_threads(s.numSiteTable(), a, unordered_ssed_3mat_evaluate_func);

    // for(int j=0; j < s.numSiteTable(); ++j) {
    //  int i = tab[j];
    //  REAL64 *lm = (REAL64 *)&(l.elem(i).elem().elem(0,0).real());
    //  REAL64 *rm = (REAL64 *)&(r.elem(i).elem().elem(0,0).real());
    //  REAL64 *dm = (REAL64 *)&(d.elem(i).elem().elem(0,0).real());
    //
    // ssed_m_eq_hh(dm,lm,rm,1);
    //}
  }
}

//-------------------------------------------------------------------

// Specialization to optimize the case   
//    LatticeColorMatrix[Subset] += LatticeColorMatrix * LatticeColorMatrix
// Threading hoist
  struct ordered_ssed_3mat_const_args {
    REAL64* lm;
    REAL64* rm;
    REAL64* dm;
    REAL64* c;
    void (*func)(REAL64*, REAL64*, REAL64*, REAL64*, int);
  };

  inline 
  void ordered_ssed_3mat_const_evaluate_func(int lo, int hi, int myId, 
					     ordered_ssed_3mat_const_args *a)
  {
    void (*func)(REAL64*, REAL64*, REAL64*, REAL64*, int) = a->func;
    
    int n_3mat = (hi - lo);
    int index = lo*18; // 3*3*2 doubles
    REAL64* lm = &(a->lm[index]);
    REAL64* rm = &(a->rm[index]);
    REAL64* dm = &(a->dm[index]);
    REAL64* c  = a->c; // The const

    func(dm,c,lm,rm, n_3mat);
  }

  struct unordered_ssed_3mat_const_args { 
    unordered_ssed_3mat_const_args(    const OLattice<DCol>& l_,
				       const OLattice<DCol>& r_,
				       OLattice<DCol>& d_,
				       REAL64 *c_,
				       const int *tab_,
				       void (*func_)(REAL64*, REAL64*, REAL64*, REAL64*, int) ) :
      l(l_), r(r_),d(d_),c(c_),tab(tab_),func(func_) {}

    const OLattice<DCol>& l;
    const OLattice<DCol>& r;
    OLattice<DCol>& d;
    REAL64 *c;
    const int *tab;
    void (*func)(REAL64*, REAL64*, REAL64*, REAL64*, int);
  };

  inline 
  void unordered_ssed_3mat_const_evaluate_func(int lo, int hi, int myId,
					       unordered_ssed_3mat_const_args* a)
  {
    void (*func)(REAL64*, REAL64*, REAL64*, REAL64*, int) = a->func;
    for(int site=lo; site < hi; site++) { 
      int i = a->tab[site];
      REAL64 *lm = (REAL64 *)&(a->l.elem(i).elem().elem(0,0).real());
      REAL64 *rm = (REAL64 *)&(a->r.elem(i).elem().elem(0,0).real());
      REAL64 *dm = (REAL64 *)&(a->d.elem(i).elem().elem(0,0).real());
      func(dm,a->c, lm, rm, 1);
    }
  }
      

template<>
void evaluate(OLattice< DCol >& d, 
	      const OpAddAssign& op, 
	      const QDPExpr<BinaryNode<OpMultiply, 
	                    Reference<QDPType< DCol, OLattice< DCol > > >, 
	                    Reference<QDPType< DCol, OLattice< DCol > > > >,
	                    OLattice< DCol > >& rhs,
	      const Subset& s)
{
//  cout << "call single site QDP_M_peq_M_times_M" << endl;

  typedef OLattice< DCol >    C;

  const C& l = static_cast<const C&>(rhs.expression().left());
  const C& r = static_cast<const C&>(rhs.expression().right());

  REAL64 one=(REAL64)1;
  
  if( s.hasOrderedRep() ) { 
    REAL64 *lm = (REAL64 *)&(l.elem(s.start()).elem().elem(0,0).real());
    REAL64 *rm = (REAL64 *)&(r.elem(s.start()).elem().elem(0,0).real());
    REAL64 *dm = (REAL64 *)&(d.elem(s.start()).elem().elem(0,0).real());
    int n_3mat=s.end()-s.start()+1;
    ordered_ssed_3mat_const_args a = {lm,rm,dm,&one,ssed_m_peq_amm};
    dispatch_to_threads(n_3mat, a, ordered_ssed_3mat_const_evaluate_func);
  }
  else { 
    const int *tab = s.siteTable().slice();
    unordered_ssed_3mat_const_args a(l,r,d,&one,tab,ssed_m_peq_amm);
    dispatch_to_threads(s.numSiteTable(), a, unordered_ssed_3mat_const_evaluate_func);  
    // const int *tab = s.siteTable().slice();
    // for(int j=0; j < s.numSiteTable(); ++j) {
    //   int i = tab[j];
    //  REAL64 *lm = (REAL64 *)&(l.elem(i).elem().elem(0,0).real());
    //  REAL64 *rm = (REAL64 *)&(r.elem(i).elem().elem(0,0).real());
    //  REAL64 *dm = (REAL64 *)&(d.elem(i).elem().elem(0,0).real());
    //  ssed_m_peq_amm(dm,&one,lm,rm,1);
    // }
  }
}


// Specialization to optimize the case   
//    LatticeColorMatrix[Subset] += adj(LatticeColorMatrix) * LatticeColorMatrix
template<>
void evaluate(OLattice< DCol >& d, 
	      const OpAddAssign& op, 
	      const QDPExpr<BinaryNode<OpAdjMultiply, 
	                    UnaryNode<OpIdentity, Reference<QDPType< DCol, OLattice< DCol > > > >, 
	                    Reference<QDPType< DCol, OLattice< DCol > > > >,
	                    OLattice< DCol > >& rhs,
	      const Subset& s)
{
//  cout << "call single site QDP_M_peq_aM_times_M" << endl;

  typedef OLattice< DCol >    C;

  const C& l = static_cast<const C&>(rhs.expression().left().child());
  const C& r = static_cast<const C&>(rhs.expression().right());

 
  REAL64 one = (REAL64)1;
  if( s.hasOrderedRep() ) { 
      REAL64 *lm = (REAL64 *)&(l.elem(s.start()).elem().elem(0,0).real());
      REAL64 *rm = (REAL64 *)&(r.elem(s.start()).elem().elem(0,0).real());
      REAL64 *dm = (REAL64 *)&(d.elem(s.start()).elem().elem(0,0).real());

      int n_3mat=s.end()-s.start()+1;
      ordered_ssed_3mat_const_args a = {lm,rm,dm,&one,ssed_m_peq_ahm};
      dispatch_to_threads(n_3mat, a, ordered_ssed_3mat_const_evaluate_func);
      // ssed_m_peq_ahm(dm,&one,lm,rm,n_3mat);
  }
  else { 
    const int *tab = s.siteTable().slice();
    unordered_ssed_3mat_const_args a(l,r,d,&one,tab,ssed_m_peq_ahm);
    dispatch_to_threads(s.numSiteTable(), a, unordered_ssed_3mat_const_evaluate_func);  
    //const int *tab = s.siteTable().slice();
    //for(int j=0; j < s.numSiteTable(); ++j) {
    //  int i = tab[j];
    //  REAL64 *lm = (REAL64 *)&(l.elem(i).elem().elem(0,0).real());
    //  REAL64 *rm = (REAL64 *)&(r.elem(i).elem().elem(0,0).real());
    //  REAL64 *dm = (REAL64 *)&(d.elem(i).elem().elem(0,0).real());
    //  ssed_m_peq_ahm(dm,&one,lm,rm,1);
    //}
  }
}


// Specialization to optimize the case   
//    LatticeColorMatrix[Subset] += LatticeColorMatrix * adj(LatticeColorMatrix)
template<>
void evaluate(OLattice< DCol >& d, 
	      const OpAddAssign& op, 
	      const QDPExpr<BinaryNode<OpMultiplyAdj, 
	                    Reference<QDPType< DCol, OLattice< DCol > > >, 
	                    UnaryNode<OpIdentity, Reference<QDPType< DCol, OLattice< DCol > > > > >,
	                    OLattice< DCol > >& rhs,
	      const Subset& s)
{
//  cout << "call single site QDP_M_peq_M_times_aM" << endl;

  typedef OLattice< DCol >    C;

  const C& l = static_cast<const C&>(rhs.expression().left());
  const C& r = static_cast<const C&>(rhs.expression().right().child());

  REAL64 one=(REAL64)1;

  if( s.hasOrderedRep() ) { 
    REAL64 *lm = (REAL64 *)&(l.elem(s.start()).elem().elem(0,0).real());
    REAL64 *rm = (REAL64 *)&(r.elem(s.start()).elem().elem(0,0).real());
    REAL64 *dm = (REAL64 *)&(d.elem(s.start()).elem().elem(0,0).real());
    int n_3mat=s.end()-s.start()+1;
    ordered_ssed_3mat_const_args a = {lm,rm,dm,&one,ssed_m_peq_amh};
    dispatch_to_threads(n_3mat, a, ordered_ssed_3mat_const_evaluate_func);
    // ssed_m_peq_amh(dm,&one,lm,rm,n_3mat);
  }
  else { 
    const int *tab = s.siteTable().slice();
    unordered_ssed_3mat_const_args a(l,r,d,&one,tab,ssed_m_peq_amh);
    dispatch_to_threads(s.numSiteTable(), a, unordered_ssed_3mat_const_evaluate_func);  
    // const int *tab = s.siteTable().slice();
    // for(int j=0; j < s.numSiteTable(); ++j) {
    //  int i = tab[j];
    //  REAL64 *lm = (REAL64 *)&(l.elem(i).elem().elem(0,0).real());
    //  REAL64 *rm = (REAL64 *)&(r.elem(i).elem().elem(0,0).real());
    //  REAL64 *dm = (REAL64 *)&(d.elem(i).elem().elem(0,0).real());
    //  ssed_m_peq_amh(dm,&one,lm,rm,1);
    // }
  }
}


// Specialization to optimize the case   
//    LatticeColorMatrix[Subset] += adj(LatticeColorMatrix) * adj(LatticeColorMatrix)
template<>
void evaluate(OLattice< DCol >& d, 
	      const OpAddAssign& op, 
	      const QDPExpr<BinaryNode<OpAdjMultiplyAdj, 
	                    UnaryNode<OpIdentity, Reference<QDPType< DCol, OLattice< DCol > > > >,
	                    UnaryNode<OpIdentity, Reference<QDPType< DCol, OLattice< DCol > > > > >,
	                    OLattice< DCol > >& rhs,
	      const Subset& s)
{
//  cout << "call single site QDP_M_peq_Ma_times_Ma" << endl;

  typedef OLattice< DCol >    C;

  const C& l = static_cast<const C&>(rhs.expression().left().child());
  const C& r = static_cast<const C&>(rhs.expression().right().child());
  
  REAL64 one=(REAL64)1;
  

  if( s.hasOrderedRep() ) { 
      
    REAL64 *lm = (REAL64 *)&(l.elem(s.start()).elem().elem(0,0).real());
    REAL64 *rm = (REAL64 *)&(r.elem(s.start()).elem().elem(0,0).real());
    REAL64 *dm = (REAL64 *)&(d.elem(s.start()).elem().elem(0,0).real());
    int n_3mat=s.end()-s.start()+1;
    ordered_ssed_3mat_const_args a = {lm,rm,dm,&one,ssed_m_peq_ahh};
    dispatch_to_threads(n_3mat, a, ordered_ssed_3mat_const_evaluate_func);

    //ssed_m_peq_ahh(dm,&one,lm,rm,n_3mat);
      
  }
  else { 
    const int *tab = s.siteTable().slice();
    unordered_ssed_3mat_const_args a(l,r,d,&one,tab,ssed_m_peq_ahh);
    dispatch_to_threads(s.numSiteTable(), a, unordered_ssed_3mat_const_evaluate_func);  
    // const int *tab = s.siteTable().slice();
    // for(int j=0; j < s.numSiteTable(); ++j) {
    //  int i = tab[j];
    //  REAL64 *lm = (REAL64 *)&(l.elem(i).elem().elem(0,0).real());
    //  REAL64 *rm = (REAL64 *)&(r.elem(i).elem().elem(0,0).real());
    //  REAL64 *dm = (REAL64 *)&(d.elem(i).elem().elem(0,0).real());
    //  ssed_m_peq_ahh(dm,&one,lm,rm,1);
    // }
  }
}

//-------------------------------------------------------------------
// Specialization to optimize the case
//   LatticeColorMatrix = LatticeColorMatrix
//   Implement as m1 = a*m2  (a=1)

  struct ordered_ssed_2mat_const_args {
    REAL64* src;
    REAL64* target;
    REAL64* c;
    void (*func)(REAL64*, REAL64*, REAL64*, int);
  };

  inline 
  void ordered_ssed_2mat_const_evaluate_func(int lo, int hi, int myId, 
					     ordered_ssed_2mat_const_args *a)
  {
    void (*func)(REAL64*, REAL64*, REAL64*, int) = a->func;
    
    int n_3mat = (hi - lo);
    int index = lo*18; // 3*3*2 doubles
    REAL64* src = &(a->src[index]);
    REAL64* target = &(a->target[index]);
    REAL64* c  = a->c; // The const

    func(target,c,src, n_3mat);
  }

  struct unordered_ssed_2mat_const_args { 
    unordered_ssed_2mat_const_args(
				   const OLattice<DCol>& src_,
				   OLattice<DCol>& target_,
				   REAL64*c_,
				   const int* tab_,
				   void  (*func_)(REAL64*, REAL64*, REAL64*, int) ) :
      src(src_), target(target_),c(c_),tab(tab_),func(func_) {}

    const OLattice<DCol>& src;
    OLattice<DCol>& target;
    REAL64*c;
    const int* tab;
    void  (*func)(REAL64*, REAL64*, REAL64*, int);
  };

  inline 
  void unordered_ssed_2mat_const_evaluate_func(int lo, int hi, int myId, 
					     unordered_ssed_2mat_const_args *a)
  {
    void (*func)(REAL64*, REAL64*, REAL64*, int) = a->func;
    REAL64* c  = a->c; // The const    

    
    // Loop through the sites
    for(int j=lo; j < hi; j++) { 
      int i = a->tab[j];
      REAL64* d_ptr =&(a->target.elem(i).elem().elem(0,0).real());
      REAL64* l_ptr =(REAL64*)&(a->src.elem(i).elem().elem(0,0).real());
      func(d_ptr,a->c,l_ptr,1);
    }


  }


template<>
void evaluate(OLattice< DCol >& d, 
	      const OpAssign& op, 
	      const QDPExpr<
	         UnaryNode<OpIdentity, Reference< QDPType< DCol, OLattice< DCol > > > >,
                 OLattice< DCol > >& rhs, 
	      const Subset& s) 
{
  typedef OLattice<DCol> C;
  const C& l = static_cast<const C&>(rhs.expression().child());
  REAL64 one=(REAL64)1;

  if( s.hasOrderedRep() ) {
    int n_3mat=s.end()-s.start()+1;
    REAL64* d_ptr =&(d.elem(s.start()).elem().elem(0,0).real());
    REAL64* l_ptr =(REAL64*)&(l.elem(s.start()).elem().elem(0,0).real());
    ordered_ssed_2mat_const_args a = {l_ptr, d_ptr, &one, ssed_m_eq_scal_m};
    dispatch_to_threads(n_3mat, a, ordered_ssed_2mat_const_evaluate_func);

    // ssed_m_eq_scal_m(d_ptr,&one,l_ptr,n_3mat);
  }
  else {
    // Unordered case 
    const int* tab = s.siteTable().slice();
    unordered_ssed_2mat_const_args a(l,d,&one,tab, ssed_m_eq_scal_m);
    dispatch_to_threads(s.numSiteTable(), a, 
			unordered_ssed_2mat_const_evaluate_func);

    
    // Loop through the sites
    // for(int j=0; j < s.numSiteTable(); j++) { 
    //  int i = tab[j];
    //  REAL64* d_ptr =&(d.elem(i).elem().elem(0,0).real());
    //  REAL64* l_ptr =(REAL64*)&(l.elem(i).elem().elem(0,0).real());
    //
    //  ssed_m_eq_scal_m(d_ptr,&one,l_ptr,1);
    // }
  }
}

//-------------------------------------------------------------------
// Specialization to optimize the case
//   LatticeColorMatrix += LatticeColorMatrix
  struct ordered_ssed_2mat_args {
    REAL64* src;
    REAL64* target;
    void (*func)(REAL64*, REAL64*, int);
  };

  inline 
  void ordered_ssed_2mat_evaluate_func(int lo, int hi, int myId, 
					     ordered_ssed_2mat_args *a)
  {
    void (*func)(REAL64*, REAL64*, int) = a->func;
    
    int n_3mat = (hi - lo);
    int index = lo*18; // 3*3*2 doubles
    REAL64* src = &(a->src[index]);
    REAL64* target = &(a->target[index]);

    func(target,src, n_3mat);
  }

  struct unordered_ssed_2mat_args { 
    unordered_ssed_2mat_args(    const OLattice<DCol>& src_,
				 OLattice<DCol>& target_,
				 const int* tab_,
				 void  (*func_)(REAL64*, REAL64*, int) ) :
      src(src_), target(target_), tab(tab_),func(func_) {}

    const OLattice<DCol>& src;
    OLattice<DCol>& target;
    const int* tab;
    void  (*func)(REAL64*, REAL64*, int);
  };

  inline 
  void unordered_ssed_2mat_evaluate_func(int lo, int hi, int myId, 
					 unordered_ssed_2mat_args *a)
  {
    void (*func)(REAL64*, REAL64*, int) = a->func;

    
    // Loop through the sites
    for(int j=lo; j < hi; j++) { 
      int i = a->tab[j];
      REAL64* d_ptr =&(a->target.elem(i).elem().elem(0,0).real());
      REAL64* l_ptr =(REAL64*)&(a->src.elem(i).elem().elem(0,0).real());
      func(d_ptr,l_ptr,1);
    }


  }

template<>
void evaluate(OLattice< DCol >& d, 
	      const OpAddAssign& op, 
	      const QDPExpr<
	         UnaryNode<OpIdentity, Reference< QDPType< DCol, OLattice< DCol > > > >,
                 OLattice< DCol > >& rhs, 
	      const Subset& s) 
{
  typedef OLattice<DCol> C;
  const C& l = static_cast<const C&>(rhs.expression().child());

  if( s.hasOrderedRep() ) { 
    int n_3mat=s.end()-s.start()+1;
    REAL64* d_ptr =&(d.elem(s.start()).elem().elem(0,0).real());
    REAL64* l_ptr =(REAL64*)&(l.elem(s.start()).elem().elem(0,0).real());
    ordered_ssed_2mat_args a = {l_ptr, d_ptr, ssed_m_peq_m};
    dispatch_to_threads(n_3mat, a, ordered_ssed_2mat_evaluate_func);
    //    ssed_m_peq_m(d_ptr,l_ptr,n_3mat);
  }
  else {
    // Unordered case 
    const int* tab = s.siteTable().slice();
    unordered_ssed_2mat_args a(l,d, tab, ssed_m_peq_m);
    dispatch_to_threads(s.numSiteTable(), a, 
			unordered_ssed_2mat_evaluate_func);
    
    // Loop through the sites
    // for(int j=0; j < s.numSiteTable(); j++) { 
    //  int i = tab[j];
    //   REAL64* d_ptr =&(d.elem(i).elem().elem(0,0).real());
    //  REAL64* l_ptr =(REAL64*)&(l.elem(i).elem().elem(0,0).real());
    //  ssed_m_peq_m(d_ptr,l_ptr,1);
    // }

  }
}

//-------------------------------------------------------------------
// Specialization to optimize the case
//   LatticeColorMatrix -= LatticeColorMatrix
template<>
void evaluate(OLattice< DCol >& d, 
	      const OpSubtractAssign& op, 
	      const QDPExpr<
	         UnaryNode<OpIdentity, Reference< QDPType< DCol, OLattice< DCol > > > >,
                 OLattice< DCol > >& rhs, 
	      const Subset& s) 
{
  typedef OLattice<DCol> C;
  const C& l = static_cast<const C&>(rhs.expression().child());
  if (s.hasOrderedRep()) { 
    int n_3mat=s.end()-s.start()+1;
    REAL64* d_ptr =&(d.elem(s.start()).elem().elem(0,0).real());
    REAL64* l_ptr =(REAL64*)&(l.elem(s.start()).elem().elem(0,0).real());
    ordered_ssed_2mat_args a = {l_ptr, d_ptr, ssed_m_meq_m};
    dispatch_to_threads(n_3mat, a, ordered_ssed_2mat_evaluate_func);

    //   ssed_m_meq_m(d_ptr,l_ptr,n_3mat);
  }
  else {   
    // Unordered case 
    const int* tab = s.siteTable().slice();
    unordered_ssed_2mat_args a(l,d, tab, ssed_m_meq_m);
    dispatch_to_threads(s.numSiteTable(), a, 
			unordered_ssed_2mat_evaluate_func);
    
    // Loop through the sites
    //  for(int j=0; j < s.numSiteTable(); j++) { 
    //  int i = tab[j];
    //  REAL64* d_ptr =&(d.elem(i).elem().elem(0,0).real());
    //  REAL64* l_ptr =(REAL64*)&(l.elem(i).elem().elem(0,0).real());
    //  ssed_m_meq_m(d_ptr,l_ptr,1);
    // }
  }
}

// Specialization to optimize the case   
//    LatticeColorMatrix[Subset] -= LatticeColorMatrix * LatticeColorMatrix
template<>
void evaluate(OLattice< DCol >& d, 
	      const OpSubtractAssign& op, 
	      const QDPExpr<BinaryNode<OpMultiply, 
	                    Reference<QDPType< DCol, OLattice< DCol > > >, 
	                    Reference<QDPType< DCol, OLattice< DCol > > > >,
	                    OLattice< DCol > >& rhs,
	      const Subset& s)
{
//  cout << "call single site QDP_M_meq_M_times_M" << endl;

  typedef OLattice< DCol >    C;

  const C& l = static_cast<const C&>(rhs.expression().left());
  const C& r = static_cast<const C&>(rhs.expression().right());

  REAL64 mone=(REAL64)-1;


  if( s.hasOrderedRep() ) { 

    REAL64 *lm = (REAL64 *)&(l.elem(s.start()).elem().elem(0,0).real());
    REAL64 *rm = (REAL64 *)&(r.elem(s.start()).elem().elem(0,0).real());
    REAL64 *dm = (REAL64 *)&(d.elem(s.start()).elem().elem(0,0).real());
    int n_3mat=s.end()-s.start()+1;
    ordered_ssed_3mat_const_args a = {lm,rm,dm,&mone,ssed_m_peq_amm};
    dispatch_to_threads(n_3mat, a, ordered_ssed_3mat_const_evaluate_func);

    //ssed_m_peq_amm(dm,&mone,lm, rm, n_3mat);
  }
  else { 

    const int *tab = s.siteTable().slice();
    unordered_ssed_3mat_const_args a(l,r,d,&mone,tab,ssed_m_peq_amm);
    dispatch_to_threads(s.numSiteTable(), a, unordered_ssed_3mat_const_evaluate_func);  

    //    for(int j=0; j < s.numSiteTable(); j++) { 
    //      int i=tab[j];
    //      REAL64 *lm = (REAL64 *)&(l.elem(i).elem().elem(0,0).real());
    //      REAL64 *rm = (REAL64 *)&(r.elem(i).elem().elem(0,0).real());
    //       REAL64 *dm = (REAL64 *)&(d.elem(i).elem().elem(0,0).real());
    //
    //      ssed_m_peq_amm(dm,&mone,lm,rm,1);
    //   }
  }
}


// Specialization to optimize the case   
//    LatticeColorMatrix[Subset] -= adj(LatticeColorMatrix) * LatticeColorMatrix
template<>
void evaluate(OLattice< DCol >& d, 
	      const OpSubtractAssign& op, 
	      const QDPExpr<BinaryNode<OpAdjMultiply, 
	                    UnaryNode<OpIdentity, Reference<QDPType< DCol, OLattice< DCol > > > >, 
	                    Reference<QDPType< DCol, OLattice< DCol > > > >,
	                    OLattice< DCol > >& rhs,
	      const Subset& s)
{
//  cout << "call single site QDP_M_meq_aM_times_M" << endl;

  typedef OLattice< DCol >    C;

  const C& l = static_cast<const C&>(rhs.expression().left().child());
  const C& r = static_cast<const C&>(rhs.expression().right());

  REAL64 mone=(REAL64)-1;

  if( s.hasOrderedRep() ) { 
    REAL64 *lm = (REAL64 *)&(l.elem(s.start()).elem().elem(0,0).real());
    REAL64 *rm = (REAL64 *)&(r.elem(s.start()).elem().elem(0,0).real());
    REAL64 *dm = (REAL64 *)&(d.elem(s.start()).elem().elem(0,0).real());
    int n_3mat=s.end()-s.start()+1;
    ordered_ssed_3mat_const_args a = {lm,rm,dm,&mone,ssed_m_peq_ahm};
    dispatch_to_threads(n_3mat, a, ordered_ssed_3mat_const_evaluate_func);

    //ssed_m_peq_ahm(dm,&mone,lm,rm,n_3mat);
  }
  else { 

    const int *tab = s.siteTable().slice();
    unordered_ssed_3mat_const_args a(l,r,d,&mone,tab,ssed_m_peq_ahm);
    dispatch_to_threads(s.numSiteTable(), a, unordered_ssed_3mat_const_evaluate_func);  

    //for(int j=0; j < s.numSiteTable(); j++) { 
    //  int i=tab[j];
    //  REAL64 *lm = (REAL64 *)&(l.elem(i).elem().elem(0,0).real());
    //  REAL64 *rm = (REAL64 *)&(r.elem(i).elem().elem(0,0).real());
    //  REAL64 *dm = (REAL64 *)&(d.elem(i).elem().elem(0,0).real());
    //  ssed_m_peq_ahm(dm,&mone,lm, rm, 1);
    // }
  }
}


// Specialization to optimize the case   
//    LatticeColorMatrix[Subset] -= LatticeColorMatrix * adj(LatticeColorMatrix)
template<>
void evaluate(OLattice< DCol >& d, 
	      const OpSubtractAssign& op, 
	      const QDPExpr<BinaryNode<OpMultiplyAdj, 
	                    Reference<QDPType< DCol, OLattice< DCol > > >, 
	                    UnaryNode<OpIdentity, Reference<QDPType< DCol, OLattice< DCol > > > > >,
	                    OLattice< DCol > >& rhs,
	      const Subset& s)
{
//  cout << "call single site QDP_M_meq_M_times_aM" << endl;

  typedef OLattice< DCol >    C;

  const C& l = static_cast<const C&>(rhs.expression().left());
  const C& r = static_cast<const C&>(rhs.expression().right().child());

  REAL64 mone=(REAL64)-1;

  if( s.hasOrderedRep() ) { 
    REAL64 *lm = (REAL64 *)&(l.elem(s.start()).elem().elem(0,0).real());
    REAL64 *rm = (REAL64 *)&(r.elem(s.start()).elem().elem(0,0).real());
    REAL64 *dm = (REAL64 *)&(d.elem(s.start()).elem().elem(0,0).real());
    int n_3mat=s.end()-s.start()+1;
    ordered_ssed_3mat_const_args a = {lm,rm,dm,&mone,ssed_m_peq_amh};
    dispatch_to_threads(n_3mat, a, ordered_ssed_3mat_const_evaluate_func);

    //   ssed_m_peq_amh(dm,&mone,lm, rm,n_3mat);
  }
  else { 
    const int *tab = s.siteTable().slice();
    unordered_ssed_3mat_const_args a(l,r,d,&mone,tab,ssed_m_peq_amh);
    dispatch_to_threads(s.numSiteTable(), a, unordered_ssed_3mat_const_evaluate_func);  

    //  for(int j=0; j < s.numSiteTable(); j++) { 
    //    int i=tab[j];
    //    REAL64 *lm = (REAL64 *)&(l.elem(i).elem().elem(0,0).real());
    //    REAL64 *rm = (REAL64 *)&(r.elem(i).elem().elem(0,0).real());
    //    REAL64 *dm = (REAL64 *)&(d.elem(i).elem().elem(0,0).real());
    //
    //    ssed_m_peq_amh(dm,&mone,lm, rm,1);
    //  }
  }
}


// Specialization to optimize the case   
//    LatticeColorMatrix[Subset] -= adj(LatticeColorMatrix) * adj(LatticeColorMatrix)
template<>
void evaluate(OLattice< DCol >& d, 
	      const OpSubtractAssign& op, 
	      const QDPExpr<BinaryNode<OpAdjMultiplyAdj, 
	                    UnaryNode<OpIdentity, Reference<QDPType< DCol, OLattice< DCol > > > >,
	                    UnaryNode<OpIdentity, Reference<QDPType< DCol, OLattice< DCol > > > > >,
	                    OLattice< DCol > >& rhs,
	      const Subset& s)
{
//  cout << "call single site QDP_M_meq_Ma_times_Ma" << endl;

  typedef OLattice< DCol >    C;

  const C& l = static_cast<const C&>(rhs.expression().left().child());
  const C& r = static_cast<const C&>(rhs.expression().right().child());
  REAL64 mone=(REAL64)-1;

  if( s.hasOrderedRep() ) { 

    REAL64 *lm = (REAL64 *)&(l.elem(s.start()).elem().elem(0,0).real());
    REAL64 *rm = (REAL64 *)&(r.elem(s.start()).elem().elem(0,0).real());
    REAL64 *dm = (REAL64 *)&(d.elem(s.start()).elem().elem(0,0).real());
    int n_3mat=s.end()-s.start()+1;
    ordered_ssed_3mat_const_args a = {lm,rm,dm,&mone,ssed_m_peq_ahh};
    dispatch_to_threads(n_3mat, a, ordered_ssed_3mat_const_evaluate_func);

    //    ssed_m_peq_ahh(dm,&mone,lm,rm,n_3mat);
  }
  else { 
    const int *tab = s.siteTable().slice();
    unordered_ssed_3mat_const_args a(l,r,d,&mone,tab,ssed_m_peq_ahh);
    dispatch_to_threads(s.numSiteTable(), a, unordered_ssed_3mat_const_evaluate_func);  

    //  for(int j=0; j < s.numSiteTable(); j++) { 
    //    int i=tab[j];
    //    REAL64 *lm = (REAL64 *)&(l.elem(i).elem().elem(0,0).real());
    //    REAL64 *rm = (REAL64 *)&(r.elem(i).elem().elem(0,0).real());
    //    REAL64 *dm = (REAL64 *)&(d.elem(i).elem().elem(0,0).real());
    //    ssed_m_peq_ahh(dm,&mone,lm,rm,1);
    //  }
  }
}





} // namespace QDP;
