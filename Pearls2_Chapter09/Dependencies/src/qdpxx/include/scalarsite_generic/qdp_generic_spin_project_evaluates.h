#ifndef QDP_GENERIC_SPIN_PROJECT_EVALUTATES_H
#define QDP_GENERIC_SPIN_PROJECT_EVALUTATES_H

using namespace QDP;
namespace QDP {

// Typedefs
typedef PSpinVector< PColorVector< RComplex<REAL>, 3>, 2 > HVec;
typedef PSpinVector< PColorVector< RComplex<REAL>, 3>, 4> FVec;

// Four spinor (Ns * Nc * Ncomplex ) Ncomplex fastest
typedef REAL SpinColFull[4][3][2];

// Half spinor (Ns/2 * Nc * Ncomplex ) Ncomplex fastest
typedef REAL SpinColHalf[2][3][2];
// d = SpinProjectDir0Plus(Vec);

////////////////////////////////
// Threading evaluates
//
// by Xu Guo, EPCC, 28 August, 2008
////////////////////////////////

// the wrappers for the functions to be threaded
#include "qdp_generic_spin_project_evaluates_wrapper.h"


template<class A, class B>
inline
void evaluate(OLattice< HVec >& b,
              const OpAssign& op,
              const QDPExpr<                             
	              UnaryNode< FnSpinProjectDir0Plus, 
	      Reference< QDPType<FVec,OLattice< FVec > > > >,
	      OLattice< HVec > > &rhs,
	      const Subset& s) 
{

  //  Get at pointer for 4 vec
  const OLattice< FVec >& a = static_cast<const OLattice< FVec > &>(rhs.expression().child());

  if(s.hasOrderedRep()) { 
    REAL *aptr =(REAL *)&(a.elem(s.start()).elem(0).elem(0).real());
    REAL *bptr =(REAL *)&(b.elem(s.start()).elem(0).elem(0).real());

    int total_n_vec = s.end() - s.start() +1;

    ordered_spin_project_user_arg arg = {aptr, bptr, inlineSpinProjDir0Plus};

    dispatch_to_threads(total_n_vec, arg, ordered_spin_project_evaluate_function);

    ///////////////////
    // Original code
    ////////////////////
    //unsigned int n_vec=s.end() - s.start()+1;

    //inlineSpinProjDir0Plus(aptr, bptr, n_vec);    
  }
  else{ 
    const int* tab = s.siteTable().slice();

    int totalSize = s.numSiteTable();

    unordered_spin_project_user_arg<A, B> arg = {a, b, tab, inlineSpinProjDir0Plus};

    dispatch_to_threads(totalSize, arg, unordered_spin_project_evaluate_function);

    ///////////////////
    // Original code
    ////////////////////
    /*
    for(int j=0; j < s.numSiteTable(); j++) { 
      int i = tab[j];

      REAL *aptr =(REAL *)&(a.elem(i).elem(0).elem(0).real());
      REAL *bptr =(REAL *)&(b.elem(i).elem(0).elem(0).real());
      inlineSpinProjDir0Plus(aptr, bptr, 1);    
      }*/
  }

}

// d = SpinProjectDir1Plus(Vec);
template<class A, class B>
inline
void evaluate(OLattice< HVec >& b,
              const OpAssign& op,
              const QDPExpr<                             
	              UnaryNode< FnSpinProjectDir1Plus, 
	      Reference< QDPType<FVec,OLattice< FVec > > > >,
	      OLattice< HVec > > &rhs,
	      const Subset& s) 
{

  
  const OLattice< FVec >& a = static_cast<const OLattice< FVec > &>(rhs.expression().child());

  if( s.hasOrderedRep() ) { 
    REAL *aptr =(REAL *)&(a.elem(s.start()).elem(0).elem(0).real());
    REAL *bptr =(REAL *)&(b.elem(s.start()).elem(0).elem(0).real());

    
    int total_n_vec = s.end() - s.start() +1;

    ordered_spin_project_user_arg arg = {aptr, bptr, inlineSpinProjDir1Plus};

    dispatch_to_threads(total_n_vec, arg, ordered_spin_project_evaluate_function);
    
    ///////////////////
    // Original code
    ////////////////////
    //unsigned int n_vec=s.end() - s.start()+1;
    //inlineSpinProjDir1Plus(aptr, bptr, n_vec);
  }
  else{ 
    const int* tab = s.siteTable().slice();
    
    int totalSize = s.numSiteTable();

    unordered_spin_project_user_arg<A, B> arg = {a, b, tab,  inlineSpinProjDir1Plus};

    dispatch_to_threads(totalSize, arg, unordered_spin_project_evaluate_function);
    
    ///////////////////
    // Original code
    ////////////////////
    /*
    for(int j=0; j < s.numSiteTable(); j++) { 
      int i = tab[j];

      REAL *aptr =(REAL *)&(a.elem(i).elem(0).elem(0).real());
      REAL *bptr =(REAL *)&(b.elem(i).elem(0).elem(0).real());
      inlineSpinProjDir1Plus(aptr, bptr, 1);    
      }*/
  }

}

// d = SpinProjectDir2Plus(Vec);
template<class A, class B>
inline
void evaluate(OLattice< HVec >& b,
              const OpAssign& op,
              const QDPExpr<                             
	              UnaryNode< FnSpinProjectDir2Plus, 
	      Reference< QDPType<FVec,OLattice< FVec > > > >,
	      OLattice< HVec > > &rhs,
	      const Subset& s) 
{
  
  const OLattice< FVec >& a = static_cast<const OLattice< FVec > &>(rhs.expression().child());

  if( s.hasOrderedRep() ) { 
    REAL *aptr =(REAL *)&(a.elem(s.start()).elem(0).elem(0).real());
    REAL *bptr =(REAL *)&(b.elem(s.start()).elem(0).elem(0).real());
    
    int total_n_vec = s.end() - s.start() +1;

    ordered_spin_project_user_arg arg = {aptr, bptr, inlineSpinProjDir2Plus};

    dispatch_to_threads(total_n_vec, arg, ordered_spin_project_evaluate_function);
    
    ///////////////////
    // Original code
    ////////////////////    
    //unsigned int n_vec=s.end() - s.start()+1;
    //inlineSpinProjDir2Plus(aptr, bptr, n_vec);
  }
  else { 
    const int* tab = s.siteTable().slice();
    
    int totalSize = s.numSiteTable();

    unordered_spin_project_user_arg<A, B> arg = {a, b, tab,  inlineSpinProjDir2Plus};

    dispatch_to_threads(totalSize, arg, unordered_spin_project_evaluate_function);
    
    ///////////////////
    // Original code
    ////////////////////
    /*
    for(int j=0; j < s.numSiteTable(); j++) { 
      int i = tab[j];

      REAL *aptr =(REAL *)&(a.elem(i).elem(0).elem(0).real());
      REAL *bptr =(REAL *)&(b.elem(i).elem(0).elem(0).real());
      inlineSpinProjDir2Plus(aptr, bptr, 1);    
      }*
  }

}

// d = SpinProjectDir3Plus(Vec);
template<class A, class B>
inline
void evaluate(OLattice< HVec >& b,
              const OpAssign& op,
              const QDPExpr<                             
	              UnaryNode< FnSpinProjectDir3Plus, 
	      Reference< QDPType<FVec,OLattice< FVec > > > >,
	      OLattice< HVec > > &rhs,
	      const Subset& s) 
{
  const OLattice< FVec >& a = static_cast<const OLattice< FVec > &>(rhs.expression().child());

  if( s.hasOrderedRep() ) {
    REAL *aptr =(REAL *)&(a.elem(s.start()).elem(0).elem(0).real());
    REAL *bptr =(REAL *)&(b.elem(s.start()).elem(0).elem(0).real());

    int total_n_vec = s.end() - s.start() +1;

    ordered_spin_project_user_arg arg = {aptr, bptr, inlineSpinProjDir3Plus};

    dispatch_to_threads(total_n_vec, arg, ordered_spin_project_evaluate_function);
    
    ///////////////////
    // Original code
    ////////////////////    
    //unsigned int n_vec=s.end() - s.start()+1;
    //inlineSpinProjDir3Plus(aptr, bptr, n_vec);
  }
  else { 
    const int* tab = s.siteTable().slice();
    
    int totalSize = s.numSiteTable();

    unordered_spin_project_user_arg<A, B> arg = {a, b, tab,  inlineSpinProjDir3Plus};

    dispatch_to_threads(totalSize, arg, unordered_spin_project_evaluate_function);
    
    ///////////////////
    // Original code
    ////////////////////
    /*
    for(int j=0; j < s.numSiteTable(); j++) { 
      int i = tab[j];

      REAL *aptr =(REAL *)&(a.elem(i).elem(0).elem(0).real());
      REAL *bptr =(REAL *)&(b.elem(i).elem(0).elem(0).real());
      inlineSpinProjDir3Plus(aptr, bptr, 1);    
      }*
  }

}

// d = SpinProjectDir0Minus(Vec);
template<class A, class B>
inline
void evaluate(OLattice< HVec >& b,
              const OpAssign& op,
              const QDPExpr<                             
	              UnaryNode< FnSpinProjectDir0Minus, 
	      Reference< QDPType<FVec,OLattice< FVec > > > >,
	      OLattice< HVec > > &rhs,
	      const Subset& s) 
{

  //  Get at pointer for 4 vec
  const OLattice< FVec >& a = static_cast<const OLattice< FVec > &>(rhs.expression().child());

  if( s.hasOrderedRep() ) { 
    REAL *aptr =(REAL *)&(a.elem(s.start()).elem(0).elem(0).real());
    REAL *bptr =(REAL *)&(b.elem(s.start()).elem(0).elem(0).real());
    
    int total_n_vec = s.end() - s.start() +1;

    ordered_spin_project_user_arg arg = {aptr, bptr, inlineSpinProjDir0Minus};

    dispatch_to_threads(total_n_vec, arg, ordered_spin_project_evaluate_function);
    
    ///////////////////
    // Original code
    ////////////////////
    //unsigned int n_vec=s.end() - s.start()+1;
    //inlineSpinProjDir0Minus(aptr, bptr, n_vec);
  }
  else { 
    const int* tab = s.siteTable().slice();

    int totalSize = s.numSiteTable();

    unordered_spin_project_user_arg<A, B> arg = {a, b, tab,  inlineSpinProjDir0Minus};

    dispatch_to_threads(totalSize, arg, unordered_spin_project_evaluate_function);
    
    ///////////////////
    // Original code
    ////////////////////
    /*
    for(int j=0; j < s.numSiteTable(); j++) { 
      int i = tab[j];
      
      REAL *aptr =(REAL *)&(a.elem(i).elem(0).elem(0).real());
      REAL *bptr =(REAL *)&(b.elem(i).elem(0).elem(0).real());
      inlineSpinProjDir0Minus(aptr, bptr, 1);    
      }*/
  }

}

// d = SpinProjectDir1Minus(Vec);
template<class A, class B>
inline
void evaluate(OLattice< HVec >& b,
              const OpAssign& op,
              const QDPExpr<                             
	              UnaryNode< FnSpinProjectDir1Minus, 
	      Reference< QDPType<FVec,OLattice< FVec > > > >,
	      OLattice< HVec > > &rhs,
	      const Subset& s) 
{

  const OLattice< FVec >& a = static_cast<const OLattice< FVec > &>(rhs.expression().child());

  if( s.hasOrderedRep() ) { 
    REAL *aptr =(REAL *)&(a.elem(s.start()).elem(0).elem(0).real());
    REAL *bptr =(REAL *)&(b.elem(s.start()).elem(0).elem(0).real());
    
    int total_n_vec = s.end() - s.start() +1;

    ordered_spin_project_user_arg arg = {aptr, bptr, inlineSpinProjDir1Minus};

    dispatch_to_threads(total_n_vec, arg, ordered_spin_project_evaluate_function);
    
    ///////////////////
    // Original code
    ////////////////////
    //unsigned int n_vec=s.end() - s.start()+1;
    //inlineSpinProjDir1Minus(aptr, bptr, n_vec);
  }
  else { 
    const int* tab = s.siteTable().slice();

    
    int totalSize = s.numSiteTable();
 
    unordered_spin_project_user_arg<A, B> arg = {a, b, tab, inlineSpinProjDir1Minus};

    dispatch_to_threads(totalSize, arg, unordered_spin_project_evaluate_function);
    
    ///////////////////
    // Original code
    ////////////////////
    /*
      for(int j=0; j < s.numSiteTable(); j++) { 
      int i = tab[j];
      
      REAL *aptr =(REAL *)&(a.elem(i).elem(0).elem(0).real());
      REAL *bptr =(REAL *)&(b.elem(i).elem(0).elem(0).real());
      inlineSpinProjDir1Minus(aptr, bptr, 1);    
      }*/
  }

}

// d = SpinProjectDir2Minus(Vec);
template<class A, class B>
inline
void evaluate(OLattice< HVec >& b,
              const OpAssign& op,
              const QDPExpr<                             
	              UnaryNode< FnSpinProjectDir2Minus, 
	      Reference< QDPType<FVec,OLattice< FVec > > > >,
	      OLattice< HVec > > &rhs,
	      const Subset& s) 
{

  
  const OLattice< FVec >& a = static_cast<const OLattice< FVec > &>(rhs.expression().child());

  if( s.hasOrderedRep() ) { 
    REAL *aptr =(REAL *)&(a.elem(s.start()).elem(0).elem(0).real());
    REAL *bptr =(REAL *)&(b.elem(s.start()).elem(0).elem(0).real());
    
    int total_n_vec = s.end() - s.start() +1;

    ordered_spin_project_user_arg arg = {aptr, bptr, inlineSpinProjDir2Minus};

    dispatch_to_threads(total_n_vec, arg, ordered_spin_project_evaluate_function);
    
    ///////////////////
    // Original code
    ////////////////////   
    //unsigned int n_vec=s.end() - s.start()+1;
    //inlineSpinProjDir2Minus(aptr, bptr, n_vec);
  }
  else { 
    const int* tab = s.siteTable().slice();

    int totalSize = s.numSiteTable();
 
    unordered_spin_project_user_arg<A, B> arg = {a, b, tab,  inlineSpinProjDir2Minus};

    dispatch_to_threads(totalSize, arg, unordered_spin_project_evaluate_function);
    
    ///////////////////
    // Original code
    ////////////////////
    /*  
    for(int j=0; j < s.numSiteTable(); j++) { 
      int i = tab[j];
      
      REAL *aptr =(REAL *)&(a.elem(i).elem(0).elem(0).real());
      REAL *bptr =(REAL *)&(b.elem(i).elem(0).elem(0).real());
      inlineSpinProjDir2Minus(aptr, bptr, 1);    
      }*/
  }


}

// d = SpinProjectDir3Minus(Vec);
template<class A, class B>
inline
void evaluate(OLattice< HVec >& b,
              const OpAssign& op,
              const QDPExpr<                             
	              UnaryNode< FnSpinProjectDir3Minus, 
	      Reference< QDPType<FVec,OLattice< FVec > > > >,
	      OLattice< HVec > > &rhs,
	      const Subset& s) 
{
  const OLattice< FVec >& a = static_cast<const OLattice< FVec > &>(rhs.expression().child());

  if( s.hasOrderedRep() ) { 

    REAL *aptr =(REAL *)&(a.elem(s.start()).elem(0).elem(0).real());
    REAL *bptr =(REAL *)&(b.elem(s.start()).elem(0).elem(0).real());

    int total_n_vec = s.end() - s.start() +1;

    ordered_spin_project_user_arg arg = {aptr, bptr, inlineSpinProjDir3Minus};

    dispatch_to_threads(total_n_vec, arg, ordered_spin_project_evaluate_function);
    
    ///////////////////
    // Original code
    ////////////////////       
    //unsigned int n_vec=s.end() - s.start()+1;
    //inlineSpinProjDir3Minus(aptr, bptr, n_vec);
  }
  else { 
    const int* tab = s.siteTable().slice();

    int totalSize = s.numSiteTable();
 
    unordered_spin_project_user_arg<A, B> arg = {a, b, tab,  inlineSpinProjDir3Minus};

    dispatch_to_threads(totalSize, arg, unordered_spin_project_evaluate_function);
    
    ///////////////////
    // Original code
    ////////////////////
    /*  
    for(int j=0; j < s.numSiteTable(); j++) { 
      int i = tab[j];
      
      REAL *aptr =(REAL *)&(a.elem(i).elem(0).elem(0).real());
      REAL *bptr =(REAL *)&(b.elem(i).elem(0).elem(0).real());
      inlineSpinProjDir3Minus(aptr, bptr, 1);    
      }*/
  }

}




// d = SpinReconstructDir0Plus(Vec);
template<class A, class B>
inline
void evaluate(OLattice< FVec >& b,
              const OpAssign& op,
              const QDPExpr<                             
	              UnaryNode< FnSpinReconstructDir0Plus, 
	      Reference< QDPType<HVec,OLattice< HVec > > > >,
	      OLattice< FVec > > &rhs,
	      const Subset& s) 
{

  const OLattice< HVec >& a = static_cast<const OLattice< HVec > &>(rhs.expression().child());

  if( s.hasOrderedRep() ) { 

    REAL *aptr =(REAL *)&(a.elem(s.start()).elem(0).elem(0).real());
    REAL *bptr =(REAL *)&(b.elem(s.start()).elem(0).elem(0).real());
    
    int total_n_vec = s.end() - s.start() +1;

    ordered_spin_project_user_arg arg = {aptr, bptr, inlineSpinReconDir0Plus};

    dispatch_to_threads(total_n_vec, arg, ordered_spin_project_evaluate_function);
    
    ///////////////////
    // Original code
    ////////////////////       
    //unsigned int n_vec=s.end() - s.start()+1;
    //inlineSpinReconDir0Plus(aptr, bptr, n_vec);
  }
  else { 
    const int* tab = s.siteTable().slice();
    
    int totalSize = s.numSiteTable();
 
    unordered_spin_project_user_arg<A, B> arg = {a, b, tab,  inlineSpinReconDir0Plus};

    dispatch_to_threads(totalSize, arg, unordered_spin_project_evaluate_function);
    
    ///////////////////
    // Original code
    ////////////////////
    /* 
    for(int j=0; j < s.numSiteTable(); j++) { 
      int i = tab[j];
      
      REAL *aptr =(REAL *)&(a.elem(i).elem(0).elem(0).real());
      REAL *bptr =(REAL *)&(b.elem(i).elem(0).elem(0).real());
      inlineSpinReconDir0Plus(aptr, bptr, 1);    
      }*/
  }


  
}

// d = SpinReconstructDir1Plus(Vec);
template<class A, class B>
inline
void evaluate(OLattice< FVec >& b,
              const OpAssign& op,
              const QDPExpr<                             
	              UnaryNode< FnSpinReconstructDir1Plus, 
	      Reference< QDPType<HVec,OLattice< HVec > > > >,
	      OLattice< FVec > > &rhs,
	      const Subset& s) 
{

  const OLattice< HVec >& a = static_cast<const OLattice< HVec > &>(rhs.expression().child());

  if( s.hasOrderedRep() ) { 
    REAL *aptr =(REAL *)&(a.elem(s.start()).elem(0).elem(0).real());
    REAL *bptr =(REAL *)&(b.elem(s.start()).elem(0).elem(0).real());
    
    int total_n_vec = s.end() - s.start() +1;

    ordered_spin_project_user_arg arg = {aptr, bptr, inlineSpinReconDir1Plus};

    dispatch_to_threads(total_n_vec, arg, ordered_spin_project_evaluate_function);
    
    ///////////////////
    // Original code
    ////////////////////      
    //unsigned int n_vec=s.end() - s.start()+1;
    //inlineSpinReconDir1Plus(aptr, bptr, n_vec);
  }
  else { 
    const int* tab = s.siteTable().slice();
    
    int totalSize = s.numSiteTable();
 
    unordered_spin_project_user_arg<A, B> arg = {a, b, tab,  inlineSpinReconDir1Plus};

    dispatch_to_threads(totalSize, arg, unordered_spin_project_evaluate_function);
    
    ///////////////////
    // Original code
    ////////////////////
    /* 
    for(int j=0; j < s.numSiteTable(); j++) { 
      int i = tab[j];
      
      REAL *aptr =(REAL *)&(a.elem(i).elem(0).elem(0).real());
      REAL *bptr =(REAL *)&(b.elem(i).elem(0).elem(0).real());
      inlineSpinReconDir1Plus(aptr, bptr, 1);    
      }*/
  }

}

// d = SpinReconstructDir2Plus(Vec);
template<class A, class B>
inline
void evaluate(OLattice< FVec >& b,
              const OpAssign& op,
              const QDPExpr<                             
	              UnaryNode< FnSpinReconstructDir2Plus, 
	      Reference< QDPType<HVec,OLattice< HVec > > > >,
	      OLattice< FVec > > &rhs,
	      const Subset& s) 
{


  const OLattice< HVec >& a = static_cast<const OLattice< HVec > &>(rhs.expression().child());

  if( s.hasOrderedRep() ) { 
    REAL *aptr =(REAL *)&(a.elem(s.start()).elem(0).elem(0).real());
    REAL *bptr =(REAL *)&(b.elem(s.start()).elem(0).elem(0).real());
    
    int total_n_vec = s.end() - s.start() +1;

    ordered_spin_project_user_arg arg = {aptr, bptr, inlineSpinReconDir2Plus};

    dispatch_to_threads(total_n_vec, arg, ordered_spin_project_evaluate_function);
    
    ///////////////////
    // Original code
    ////////////////////    
    //unsigned int n_vec=s.end() - s.start()+1;
    //inlineSpinReconDir2Plus(aptr, bptr, n_vec);
  }
  else { 
    const int* tab = s.siteTable().slice();

    int totalSize = s.numSiteTable();
 
    unordered_spin_project_user_arg<A, B> arg = {a, b, tab,  inlineSpinReconDir2Plus};

    dispatch_to_threads(totalSize, arg, unordered_spin_project_evaluate_function);
    
    ///////////////////
    // Original code
    ////////////////////
    /* 
    for(int j=0; j < s.numSiteTable(); j++) { 
      int i = tab[j];
      
      REAL *aptr =(REAL *)&(a.elem(i).elem(0).elem(0).real());
      REAL *bptr =(REAL *)&(b.elem(i).elem(0).elem(0).real());
      inlineSpinReconDir2Plus(aptr, bptr, 1);    
      }*/
  }

}

// d = SpinReconstructDir3Plus(Vec);
template<class A, class B>
inline
void evaluate(OLattice< FVec >& b,
              const OpAssign& op,
              const QDPExpr<                             
	              UnaryNode< FnSpinReconstructDir3Plus, 
	      Reference< QDPType<HVec,OLattice< HVec > > > >,
	      OLattice< FVec > > &rhs,
	      const Subset& s) 
{

  const OLattice< HVec >& a = static_cast<const OLattice< HVec > &>(rhs.expression().child());

  if( s.hasOrderedRep() ) { 
    REAL *aptr =(REAL *)&(a.elem(s.start()).elem(0).elem(0).real());
    REAL *bptr =(REAL *)&(b.elem(s.start()).elem(0).elem(0).real());
    
    int total_n_vec = s.end() - s.start() +1;

    ordered_spin_project_user_arg arg = {aptr, bptr, inlineSpinReconDir3Plus};

    dispatch_to_threads(total_n_vec, arg, ordered_spin_project_evaluate_function);
    
    ///////////////////
    // Original code
    ////////////////////        
    //unsigned int n_vec=s.end() - s.start()+1;
    //inlineSpinReconDir3Plus(aptr, bptr, n_vec);
  }
  else { 
    const int* tab = s.siteTable().slice();
    
    int totalSize = s.numSiteTable();
 
    unordered_spin_project_user_arg<A, B> arg = {a, b, tab,  inlineSpinReconDir3Plus};

    dispatch_to_threads(totalSize, arg, unordered_spin_project_evaluate_function);
    
    ///////////////////
    // Original code
    ////////////////////
    /* 
    for(int j=0; j < s.numSiteTable(); j++) { 
      int i = tab[j];
      
      REAL *aptr =(REAL *)&(a.elem(i).elem(0).elem(0).real());
      REAL *bptr =(REAL *)&(b.elem(i).elem(0).elem(0).real());
      inlineSpinReconDir3Plus(aptr, bptr, 1);    
      }*/
  }


}

// d = SpinReconstructDir0Minus(Vec);
template<class A, class B>
inline
void evaluate(OLattice< FVec >& b,
              const OpAssign& op,
              const QDPExpr<                             
	              UnaryNode< FnSpinReconstructDir0Minus, 
	      Reference< QDPType<HVec,OLattice< HVec > > > >,
	      OLattice< FVec > > &rhs,
	      const Subset& s) 
{

  const OLattice< HVec >& a = static_cast<const OLattice< HVec > &>(rhs.expression().child());

  if( s.hasOrderedRep() ) { 
    REAL *aptr =(REAL *)&(a.elem(s.start()).elem(0).elem(0).real());
    REAL *bptr =(REAL *)&(b.elem(s.start()).elem(0).elem(0).real());
    
    int total_n_vec = s.end() - s.start() +1;

    ordered_spin_project_user_arg arg = {aptr, bptr, inlineSpinReconDir0Minus};

    dispatch_to_threads(total_n_vec, arg, ordered_spin_project_evaluate_function);
    
    ///////////////////
    // Original code
    ////////////////////        
    //unsigned int n_vec=s.end() - s.start()+1;
    //inlineSpinReconDir0Minus(aptr, bptr, n_vec);
  }
  else { 

    const int* tab = s.siteTable().slice();

    int totalSize = s.numSiteTable();
 
    unordered_spin_project_user_arg<A, B> arg = {a, b, tab,  inlineSpinReconDir0Minus};

    dispatch_to_threads(totalSize, arg, unordered_spin_project_evaluate_function);
    
    ///////////////////
    // Original code
    ////////////////////
    /* 
    for(int j=0; j < s.numSiteTable(); j++) { 
      int i = tab[j];
      
      REAL *aptr =(REAL *)&(a.elem(i).elem(0).elem(0).real());
      REAL *bptr =(REAL *)&(b.elem(i).elem(0).elem(0).real());
      inlineSpinReconDir0Minus(aptr, bptr, 1);    
      }*/
  }
  
}

// d = SpinReconstructDir1Minus(Vec);
template<class A, class B>
inline
void evaluate(OLattice< FVec >& b,
              const OpAssign& op,
              const QDPExpr<                             
	              UnaryNode< FnSpinReconstructDir1Minus, 
	      Reference< QDPType<HVec,OLattice< HVec > > > >,
	      OLattice< FVec > > &rhs,
	      const Subset& s) 
{

  const OLattice< HVec >& a = static_cast<const OLattice< HVec > &>(rhs.expression().child());

  if( s.hasOrderedRep() ) { 
    REAL *aptr =(REAL *)&(a.elem(s.start()).elem(0).elem(0).real());
    REAL *bptr =(REAL *)&(b.elem(s.start()).elem(0).elem(0).real());

    
    int total_n_vec = s.end() - s.start() +1;

    ordered_spin_project_user_arg arg = {aptr, bptr, inlineSpinReconDir1Minus};

    dispatch_to_threads(total_n_vec, arg, ordered_spin_project_evaluate_function);
    
    ///////////////////
    // Original code
    ////////////////////       
    //unsigned int n_vec=s.end() - s.start()+1;
    //inlineSpinReconDir1Minus(aptr, bptr, n_vec);
  } 
  else { 

    const int* tab = s.siteTable().slice();

    int totalSize = s.numSiteTable();
 
    unordered_spin_project_user_arg<A, B> arg = {a, b, tab,  inlineSpinReconDir1Minus};

    dispatch_to_threads(totalSize, arg, unordered_spin_project_evaluate_function);
    
    ///////////////////
    // Original code
    ////////////////////
    /* 
    for(int j=0; j < s.numSiteTable(); j++) { 
      int i = tab[j];
      
      REAL *aptr =(REAL *)&(a.elem(i).elem(0).elem(0).real());
      REAL *bptr =(REAL *)&(b.elem(i).elem(0).elem(0).real());
      inlineSpinReconDir1Minus(aptr, bptr, 1);    
      }*/
  }
 
}

// d = SpinReconstructDir2Minus(Vec);
template<class A, class B>
inline
void evaluate(OLattice< FVec >& b,
              const OpAssign& op,
              const QDPExpr<                             
	              UnaryNode< FnSpinReconstructDir2Minus, 
	      Reference< QDPType<HVec,OLattice< HVec > > > >,
	      OLattice< FVec > > &rhs,
	      const Subset& s) 
{


  const OLattice< HVec >& a = static_cast<const OLattice< HVec > &>(rhs.expression().child());

  if( s.hasOrderedRep() ) { 
    REAL *aptr =(REAL *)&(a.elem(s.start()).elem(0).elem(0).real());
    REAL *bptr =(REAL *)&(b.elem(s.start()).elem(0).elem(0).real());
    
    int total_n_vec = s.end() - s.start() +1;

    ordered_spin_project_user_arg arg = {aptr, bptr, inlineSpinReconDir2Minus};

    dispatch_to_threads(total_n_vec, arg, ordered_spin_project_evaluate_function);
    
    ///////////////////
    // Original code
    ////////////////////       
    //unsigned int n_vec=s.end() - s.start()+1;
    //inlineSpinReconDir2Minus(aptr, bptr, n_vec);
  }
  else { 

    const int* tab = s.siteTable().slice();
    
    int totalSize = s.numSiteTable();
 
    unordered_spin_project_user_arg<A, B> arg = {a, b, tab,  inlineSpinReconDir2Minus};

    dispatch_to_threads(totalSize, arg, unordered_spin_project_evaluate_function);
    
    ///////////////////
    // Original code
    ////////////////////
    /* 
    for(int j=0; j < s.numSiteTable(); j++) { 
      int i = tab[j];
      
      REAL *aptr =(REAL *)&(a.elem(i).elem(0).elem(0).real());
      REAL *bptr =(REAL *)&(b.elem(i).elem(0).elem(0).real());
      inlineSpinReconDir2Minus(aptr, bptr, 1);    
      }*/
  }

}

// d = SpinReconstructDir3Minus(Vec);
template<class A, class B>
inline
void evaluate(OLattice< FVec >& b,
              const OpAssign& op,
              const QDPExpr<                             
	              UnaryNode< FnSpinReconstructDir3Minus, 
	      Reference< QDPType<HVec,OLattice< HVec > > > >,
	      OLattice< FVec > > &rhs,
	      const Subset& s) 
{

  const OLattice< HVec >& a = static_cast<const OLattice< HVec > &>(rhs.expression().child());

  if( s.hasOrderedRep()) { 
    REAL *aptr =(REAL *)&(a.elem(s.start()).elem(0).elem(0).real());
    REAL *bptr =(REAL *)&(b.elem(s.start()).elem(0).elem(0).real());
    
    int total_n_vec = s.end() - s.start() +1;

    ordered_spin_project_user_arg arg = {aptr, bptr, inlineSpinReconDir3Minus};

    dispatch_to_threads(total_n_vec, arg, ordered_spin_project_evaluate_function);
    
    ///////////////////
    // Original code
    ////////////////////       
    //unsigned int n_vec=s.end() - s.start()+1;
    //inlineSpinReconDir3Minus(aptr, bptr, n_vec);
  }
  else { 

    const int* tab = s.siteTable().slice();

    int totalSize = s.numSiteTable();
 
    unordered_spin_project_user_arg<A, B> arg = {a, b, tab,  inlineSpinReconDir3Minus};

    dispatch_to_threads(totalSize, arg, unordered_spin_project_evaluate_function);
    
    ///////////////////
    // Original code
    ////////////////////
    /* 
    for(int j=0; j < s.numSiteTable(); j++) { 
      int i = tab[j];
      
      REAL *aptr =(REAL *)&(a.elem(i).elem(0).elem(0).real());
      REAL *bptr =(REAL *)&(b.elem(i).elem(0).elem(0).real());
      inlineSpinReconDir3Minus(aptr, bptr, 1);    
      }*/
  }

}



// d += SpinReconstructDir0Plus(Vec);
template<class A, class B>
inline
void evaluate(OLattice< FVec >& b,
              const OpAddAssign& op,
              const QDPExpr<                             
	              UnaryNode< FnSpinReconstructDir0Plus, 
	      Reference< QDPType<HVec,OLattice< HVec > > > >,
	      OLattice< FVec > > &rhs,
	      const Subset& s) 
{


  const OLattice< HVec >& a = static_cast<const OLattice< HVec > &>(rhs.expression().child());

  if( s.hasOrderedRep() ) { 
    
    REAL *aptr =(REAL *)&(a.elem(s.start()).elem(0).elem(0).real());
    REAL *bptr =(REAL *)&(b.elem(s.start()).elem(0).elem(0).real());
        
    int total_n_vec = s.end() - s.start() +1;

    ordered_spin_project_user_arg arg = {aptr, bptr, inlineAddSpinReconDir0Plus};

    dispatch_to_threads(total_n_vec, arg, ordered_spin_project_evaluate_function);

    ///////////////////
    // Original code
    ////////////////////  
    //unsigned int n_vec=s.end() - s.start()+1;
    //inlineAddSpinReconDir0Plus(aptr, bptr, n_vec);
  }
  else { 
    
    const int* tab = s.siteTable().slice();

    int totalSize = s.numSiteTable();
 
    unordered_spin_project_user_arg<A, B> arg = {a, b, tab,  inlineAddSpinReconDir0Plus};

    dispatch_to_threads(totalSize, arg, unordered_spin_project_evaluate_function);
    
    ///////////////////
    // Original code
    ////////////////////
    /* 
    for(int j=0; j < s.numSiteTable(); j++) { 
      int i = tab[j];
      
      REAL *aptr =(REAL *)&(a.elem(i).elem(0).elem(0).real());
      REAL *bptr =(REAL *)&(b.elem(i).elem(0).elem(0).real());
      inlineAddSpinReconDir0Plus(aptr, bptr, 1);    
      }*/
  }

}

// d += SpinReconstructDir1Plus(Vec);
template<class A, class B>
inline
void evaluate(OLattice< FVec >& b,
              const OpAddAssign& op,
              const QDPExpr<                             
	              UnaryNode< FnSpinReconstructDir1Plus, 
	      Reference< QDPType<HVec,OLattice< HVec > > > >,
	      OLattice< FVec > > &rhs,
	      const Subset& s) 
{

  const OLattice< HVec >& a = static_cast<const OLattice< HVec > &>(rhs.expression().child());

  if( s.hasOrderedRep() ) { 
    REAL *aptr =(REAL *)&(a.elem(s.start()).elem(0).elem(0).real());
    REAL *bptr =(REAL *)&(b.elem(s.start()).elem(0).elem(0).real());

    int total_n_vec = s.end() - s.start() +1;

    ordered_spin_project_user_arg arg = {aptr, bptr, inlineAddSpinReconDir1Plus};

    dispatch_to_threads(total_n_vec, arg, ordered_spin_project_evaluate_function);
    
    ///////////////////
    // Original code
    ////////////////////      
    //unsigned int n_vec=s.end() - s.start()+1;
    //inlineAddSpinReconDir1Plus(aptr, bptr, n_vec);
  }
  else { 
    
    const int* tab = s.siteTable().slice();

    int totalSize = s.numSiteTable();
 
    unordered_spin_project_user_arg<A, B> arg = {a, b, tab,  inlineAddSpinReconDir1Plus};

    dispatch_to_threads(totalSize, arg, unordered_spin_project_evaluate_function);
    
    ///////////////////
    // Original code
    ////////////////////
    /* 
    for(int j=0; j < s.numSiteTable(); j++) { 
      int i = tab[j];
      
      REAL *aptr =(REAL *)&(a.elem(i).elem(0).elem(0).real());
      REAL *bptr =(REAL *)&(b.elem(i).elem(0).elem(0).real());
      inlineAddSpinReconDir1Plus(aptr, bptr, 1);    
      }*/
  }

}

// d += SpinReconstructDir2Plus(Vec);
template<class A, class B>
inline
void evaluate(OLattice< FVec >& b,
              const OpAddAssign& op,
              const QDPExpr<                             
	              UnaryNode< FnSpinReconstructDir2Plus, 
	      Reference< QDPType<HVec,OLattice< HVec > > > >,
	      OLattice< FVec > > &rhs,
	      const Subset& s) 
{
  
  const OLattice< HVec >& a = static_cast<const OLattice< HVec > &>(rhs.expression().child());

  if( s.hasOrderedRep() ) { 
    REAL *aptr =(REAL *)&(a.elem(s.start()).elem(0).elem(0).real());
    REAL *bptr =(REAL *)&(b.elem(s.start()).elem(0).elem(0).real());
        
    int total_n_vec = s.end() - s.start() +1;

    ordered_spin_project_user_arg arg = {aptr, bptr, inlineAddSpinReconDir2Plus};

    dispatch_to_threads(total_n_vec, arg, ordered_spin_project_evaluate_function);
    
    ///////////////////
    // Original code
    ////////////////////  
    //unsigned int n_vec=s.end() - s.start()+1;
    //inlineAddSpinReconDir2Plus(aptr, bptr, n_vec);
  }
  else { 
    
    const int* tab = s.siteTable().slice();

    int totalSize = s.numSiteTable();
 
    unordered_spin_project_user_arg<A, B> arg = {a, b, tab,  inlineAddSpinReconDir2Plus};

    dispatch_to_threads(totalSize, arg, unordered_spin_project_evaluate_function);
    
    ///////////////////
    // Original code
    ////////////////////
    /* 
    for(int j=0; j < s.numSiteTable(); j++) { 
      int i = tab[j];
      
      REAL *aptr =(REAL *)&(a.elem(i).elem(0).elem(0).real());
      REAL *bptr =(REAL *)&(b.elem(i).elem(0).elem(0).real());
      inlineAddSpinReconDir2Plus(aptr, bptr, 1);    
      }*/
  }

}

// d += SpinReconstructDir3Plus(Vec);
template<class A, class B>
inline
void evaluate(OLattice< FVec >& b,
              const OpAddAssign& op,
              const QDPExpr<                             
	              UnaryNode< FnSpinReconstructDir3Plus, 
	      Reference< QDPType<HVec,OLattice< HVec > > > >,
	      OLattice< FVec > > &rhs,
	      const Subset& s) 
{

  const OLattice< HVec >& a = static_cast<const OLattice< HVec > &>(rhs.expression().child());

  if( s.hasOrderedRep() ) {
    REAL *aptr =(REAL *)&(a.elem(s.start()).elem(0).elem(0).real());
    REAL *bptr =(REAL *)&(b.elem(s.start()).elem(0).elem(0).real());

    int total_n_vec = s.end() - s.start() +1;

    ordered_spin_project_user_arg arg = {aptr, bptr, inlineAddSpinReconDir3Plus};

    dispatch_to_threads(total_n_vec, arg, ordered_spin_project_evaluate_function);
    
    ///////////////////
    // Original code
    ////////////////////  
    //unsigned int n_vec=s.end() - s.start()+1;
    //inlineAddSpinReconDir3Plus(aptr, bptr, n_vec);
  }
  else { 
    
    const int* tab = s.siteTable().slice();
    
    int totalSize = s.numSiteTable();
 
    unordered_spin_project_user_arg<A, B> arg = {a, b, tab, inlineAddSpinReconDir3Plus};

    dispatch_to_threads(totalSize, arg, unordered_spin_project_evaluate_function);
    
    ///////////////////
    // Original code
    ////////////////////
    /* 
    for(int j=0; j < s.numSiteTable(); j++) { 
      int i = tab[j];
      
      REAL *aptr =(REAL *)&(a.elem(i).elem(0).elem(0).real());
      REAL *bptr =(REAL *)&(b.elem(i).elem(0).elem(0).real());
      inlineAddSpinReconDir3Plus(aptr, bptr, 1);    
      }*/
  }

}

// d += SpinReconstructDir0Minus(Vec);
template<class A, class B>
inline
void evaluate(OLattice< FVec >& b,
              const OpAddAssign& op,
              const QDPExpr<                             
	              UnaryNode< FnSpinReconstructDir0Minus, 
	      Reference< QDPType<HVec,OLattice< HVec > > > >,
	      OLattice< FVec > > &rhs,
	      const Subset& s) 
{


  const OLattice< HVec >& a = static_cast<const OLattice< HVec > &>(rhs.expression().child());

  if( s.hasOrderedRep() ) {
    REAL *aptr =(REAL *)&(a.elem(s.start()).elem(0).elem(0).real());
    REAL *bptr =(REAL *)&(b.elem(s.start()).elem(0).elem(0).real());
    
    int total_n_vec = s.end() - s.start() +1;

    ordered_spin_project_user_arg arg = {aptr, bptr, inlineAddSpinReconDir0Minus};

    dispatch_to_threads(total_n_vec, arg, ordered_spin_project_evaluate_function);
    
    ///////////////////
    // Original code
    ////////////////////  
    //unsigned int n_vec=s.end() - s.start()+1;
    //inlineAddSpinReconDir0Minus(aptr, bptr, n_vec);
  }
  else { 
    
    const int* tab = s.siteTable().slice();

    int totalSize = s.numSiteTable();
 
    unordered_spin_project_user_arg<A, B> arg = {a, b, tab,  inlineAddSpinReconDir0Minus};

    dispatch_to_threads(totalSize, arg, unordered_spin_project_evaluate_function);
    
    ///////////////////
    // Original code
    ////////////////////
    /* 
    for(int j=0; j < s.numSiteTable(); j++) { 
      int i = tab[j];
      
      REAL *aptr =(REAL *)&(a.elem(i).elem(0).elem(0).real());
      REAL *bptr =(REAL *)&(b.elem(i).elem(0).elem(0).real());
      inlineAddSpinReconDir0Minus(aptr, bptr, 1);    
      }*/
  }

  
}

// d += SpinReconstructDir1Minus(Vec);
template<class A, class B>
inline
void evaluate(OLattice< FVec >& b,
              const OpAddAssign& op,
              const QDPExpr<                             
	              UnaryNode< FnSpinReconstructDir1Minus, 
	      Reference< QDPType<HVec,OLattice< HVec > > > >,
	      OLattice< FVec > > &rhs,
	      const Subset& s) 
{

  const OLattice< HVec >& a = static_cast<const OLattice< HVec > &>(rhs.expression().child());

  if( s.hasOrderedRep()) { 
    REAL *aptr =(REAL *)&(a.elem(s.start()).elem(0).elem(0).real());
    REAL *bptr =(REAL *)&(b.elem(s.start()).elem(0).elem(0).real());
    
    int total_n_vec = s.end() - s.start() +1;

    ordered_spin_project_user_arg arg = {aptr, bptr, inlineAddSpinReconDir1Minus};

    dispatch_to_threads(total_n_vec, arg, ordered_spin_project_evaluate_function);
    
    ///////////////////
    // Original code
    ////////////////////  
    //unsigned int n_vec=s.end() - s.start()+1;
    //inlineAddSpinReconDir1Minus(aptr, bptr, n_vec);
  }
  else { 
    
    const int* tab = s.siteTable().slice();

    int totalSize = s.numSiteTable();
 
    unordered_spin_project_user_arg<A, B> arg = {a, b, tab,  inlineAddSpinReconDir1Minus};

    dispatch_to_threads(totalSize, arg, unordered_spin_project_evaluate_function);
    
    ///////////////////
    // Original code
    ////////////////////
    /* 
    for(int j=0; j < s.numSiteTable(); j++) { 
      int i = tab[j];
      
      REAL *aptr =(REAL *)&(a.elem(i).elem(0).elem(0).real());
      REAL *bptr =(REAL *)&(b.elem(i).elem(0).elem(0).real());
      inlineAddSpinReconDir1Minus(aptr, bptr, 1);    
      }*/
  }

}

// d += SpinReconstructDir2Minus(Vec);
template<class A, class B>
inline
void evaluate(OLattice< FVec >& b,
              const OpAddAssign& op,
              const QDPExpr<                             
	              UnaryNode< FnSpinReconstructDir2Minus, 
	      Reference< QDPType<HVec,OLattice< HVec > > > >,
	      OLattice< FVec > > &rhs,
	      const Subset& s) 
{


  const OLattice< HVec >& a = static_cast<const OLattice< HVec > &>(rhs.expression().child());

  if( s.hasOrderedRep() ) { 
  REAL *aptr =(REAL *)&(a.elem(s.start()).elem(0).elem(0).real());
  REAL *bptr =(REAL *)&(b.elem(s.start()).elem(0).elem(0).real());
  
  
  int total_n_vec = s.end() - s.start() +1;

  ordered_spin_project_user_arg arg = {aptr, bptr, inlineAddSpinReconDir2Minus};

  dispatch_to_threads(total_n_vec, arg, ordered_spin_project_evaluate_function);
  
  ///////////////////
  // Original code
  ////////////////////  
  //unsigned int n_vec=s.end() - s.start()+1;
  //inlineAddSpinReconDir2Minus(aptr, bptr, n_vec);
  }
  else { 
    
    const int* tab = s.siteTable().slice();
    
    int totalSize = s.numSiteTable();
 
    unordered_spin_project_user_arg<A, B> arg = {a, b, tab,  inlineAddSpinReconDir2Minus};

    dispatch_to_threads(totalSize, arg, unordered_spin_project_evaluate_function);
    
    ///////////////////
    // Original code
    ////////////////////
    /* 
    for(int j=0; j < s.numSiteTable(); j++) { 
      int i = tab[j];
      
      REAL *aptr =(REAL *)&(a.elem(i).elem(0).elem(0).real());
      REAL *bptr =(REAL *)&(b.elem(i).elem(0).elem(0).real());
      inlineAddSpinReconDir2Minus(aptr, bptr, 1);    
      }*/
  }

}

// d += SpinReconstructDir3Minus(Vec);
template<class A, class B>
inline
void evaluate(OLattice< FVec >& b,
              const OpAddAssign& op,
              const QDPExpr<                             
	              UnaryNode< FnSpinReconstructDir3Minus, 
	      Reference< QDPType<HVec,OLattice< HVec > > > >,
	      OLattice< FVec > > &rhs,
	      const Subset& s) 
{


  const OLattice< HVec >& a = static_cast<const OLattice< HVec > &>(rhs.expression().child());

  if( s.hasOrderedRep() ) { 
    REAL *aptr =(REAL *)&(a.elem(s.start()).elem(0).elem(0).real());
    REAL *bptr =(REAL *)&(b.elem(s.start()).elem(0).elem(0).real());
    
    int total_n_vec = s.end() - s.start() +1;

    ordered_spin_project_user_arg arg = {aptr, bptr, inlineAddSpinReconDir3Minus};

    dispatch_to_threads(total_n_vec, arg, ordered_spin_project_evaluate_function);
    
    ///////////////////
    // Original code
    ////////////////////    
    //unsigned int n_vec=s.end() - s.start()+1;
    //inlineAddSpinReconDir3Minus(aptr, bptr, n_vec);
  }
  else { 
    
    const int* tab = s.siteTable().slice();

    int totalSize = s.numSiteTable();
 
    unordered_spin_project_user_arg<A, B> arg = {a, b, tab,  inlineAddSpinReconDir3Minus};

    dispatch_to_threads(totalSize, arg, unordered_spin_project_evaluate_function);

    ///////////////////
    // Original code
    ////////////////////
    /*
    for(int j=0; j < s.numSiteTable(); j++) { 
      int i = tab[j];
      
      REAL *aptr =(REAL *)&(a.elem(i).elem(0).elem(0).real());
      REAL *bptr =(REAL *)&(b.elem(i).elem(0).elem(0).real());
      inlineAddSpinReconDir3Minus(aptr, bptr, 1);    
      }*/
  }

}

} // namespace QDP;

#endif
