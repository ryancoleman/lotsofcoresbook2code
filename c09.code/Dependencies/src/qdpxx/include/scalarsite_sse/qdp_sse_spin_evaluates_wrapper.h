#ifndef QDP_SSE_SPIN_EVALUATES_WRAPPER_H
#define QDP_SSE_SPIN_EVALUATES_WRAPPER_H

////////////////////////////////
// Threading evaluates wrappers
//
// by Xu Guo, EPCC, 13 October, 2008
////////////////////////////////


// for inlineSpinProjDir0Plus, inlineSpinProjDir1Plus, inlineSpinProjDir2Plus, inlineSpinProjDir3Plus, inlineSpinProjDir0Minus, inlineSpinProjDir1Minus, inlineSpinProjDir2Minus, inlineSpinProjDir3Minus
// for inlineSpinReconDir0Plus, inlineSpinReconDir1Plus, inlineSpinReconDir2Plus, inlineSpinReconDir3Plus, inlineSpinReconDir0Minus, inlineSpinReconDir1Minus, inlineSpinReconDir2Minus, inlineSpinReconDir3Minus,
// for inlineAddSpinReconDir0Plus, inlineAddSpinReconDir1Plus, inlineAddSpinReconDir2Plus, inlineAddSpinReconDir3Plus, inlineAddSpinReconDir0Minus, inlineAddSpinReconDir1Minus, inlineAddSpinReconDir2Minus, inlineAddSpinReconDir3Minus,

// user arg for evaluate having order
struct ordered_sse_spin_project_user_arg{
  REAL32* aptr;
  REAL32* bptr;
  void (*func)(const REAL32*, REAL32*, unsigned int);
};


// user func for evaluate having order
inline
void ordered_sse_spin_project_evaluate_function (int lo, int hi, int myId, ordered_sse_spin_project_user_arg* arg){

  REAL32* aptr = arg->aptr;
  REAL32* bptr = arg->bptr;
  void (*func)(const REAL32*, REAL32*, unsigned int) = arg->func; 

  unsigned int n_vec = hi - lo;

  func(aptr, bptr, n_vec);   
}


// user arg for evaluate NOT having order
template <class A, class B>
struct unordered_sse_spin_project_user_arg{
  A a;
  B b;
  const int *tab;
  void (*func)(const REAL32*, REAL32*,unsigned int);
};


// user func for evaluate NOT having order
template <class A, class B>
inline
void unordered_sse_spin_project_evaluate_function (int lo, int hi, int myId, unordered_sse_spin_project_user_arg<A, B>* arg){

  A a = arg->a;
  B b = arg->b;
  const int *tab = arg->tab;
  void (*func)(const REAL32*, REAL32*, unsigned int) = arg->func;


  for(int j=lo; j < hi; j++) { 
      int i = tab[j];

      REAL32 *aptr =(REAL *)&(a.elem(i).elem(0).elem(0).real());
      REAL32 *bptr =(REAL *)&(b.elem(i).elem(0).elem(0).real());
      func(aptr, bptr, 1);    
    }

}




#endif
