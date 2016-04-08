#ifndef QDP_SCALARSITE_GENERIC_BLAS_G5_WRAPPER_H
#define QDP_SCALARSITE_GENERIC_BLAS_G5_WRAPPER_H

////////////////////////////////
// Threading evaluates wrappers
//
// by Xu Guo, EPCC, 26 August, 2008
////////////////////////////////

//
// for vaypx3_g5: xpayz_g5ProjPlus,  xpayz_g5ProjMinus,  xmayz_g5ProjPlus,  xmayz_g5ProjMinus
//

// structure for vaypx3_g5  of having order
struct ordered_vaypx3_g5_user_arg{
  REAL* Out;
  REAL* scalep;
  REAL* Add;
  REAL* InScale;
  void (*func)(REAL*, REAL*, REAL*, REAL*, int);
};

// user func for vaypx3_g5 of having order
inline
void ordered_vaypx3_g5_evaluate_function (int lo, int hi, int myId, ordered_vaypx3_g5_user_arg* a){

  REAL* Out = a->Out;
  REAL* scalep = a->scalep;
  REAL* Add = a->Add;
  REAL* InScale = a->InScale;
  void (*func)(REAL*, REAL*, REAL*, REAL*, int) = a->func;

  int n_4vec = hi - lo;

  int index = lo * 24;
  InScale = &InScale[index];
  Add = &Add[index];
  Out = &Out[index];

  func(Out, scalep, Add, InScale, n_4vec);

}

// structure for vaypx3_g5 of NOT having order (with y only)
struct unordered_vaypx3_g5_y_user_arg{
  unordered_vaypx3_g5_y_user_arg(
				 const OLattice< TVec >& x_,
				 OLattice< TVec >& d_,
				 REAL* aptr_,
				 int Ns_,
				 const int* tab_,
				 void (*func_)(REAL*, REAL*, REAL*, REAL*, int)) : x(x_), d(d_), aptr(aptr_), Ns(Ns_), tab(tab_), func(func_) {}


  const OLattice< TVec >& x;
  OLattice< TVec >& d;
  REAL* aptr;
  int Ns;
  const int* tab;
  void (*func)(REAL*, REAL*, REAL*, REAL*, int);
};

// user func for vaypx3_g5 of NOT having order (with y only)
inline
void unordered_vaypx3_g5_y_evaluate_function (int lo, int hi, int myId, unordered_vaypx3_g5_y_user_arg* a){

  const OLattice< TVec >& x = a->x;
  OLattice< TVec >& d = a->d;
  REAL* aptr = a->aptr;
  int Ns = a->Ns;
  const int* tab = a->tab;
  void (*func)(REAL*, REAL*, REAL*, REAL*, int) = a->func;
  
  for(int j=lo; j < hi; j++) { 
    int i=tab[j];
    REAL* xptr = (REAL *)&(x.elem(i).elem(0).elem(0).real());
    REAL* yptr = (REAL *)&(d.elem(i).elem(0).elem(0).real());
    func(yptr, aptr, yptr, xptr, Ns);
  }
}
 

// structure for vaypx3_g5 of NOT having order (with z )
struct unordered_vaypx3_g5_z_user_arg{
  unordered_vaypx3_g5_z_user_arg(
				 const OLattice< TVec >& x_,
				 const OLattice< TVec >& y_,
				 OLattice< TVec >& d_,
				 REAL* aptr_,
				 int Ns_,
				 const int* tab_,
				 void (*func_)(REAL*, REAL*, REAL*, REAL*, int)) : x(x_),y(y_),d(d_),aptr(aptr_), Ns(Ns_), tab(tab_),func(func_) {}
				 
  const OLattice< TVec >& x;
  const OLattice< TVec >& y;
  OLattice< TVec >& d;
  REAL* aptr;
  int Ns;
  const int* tab;
  void (*func)(REAL*, REAL*, REAL*, REAL*, int);
};

// user func for vaypx3_g5 of NOT having order (with z)
inline
void unordered_vaypx3_g5_z_evaluate_function (int lo, int hi, int myId, unordered_vaypx3_g5_z_user_arg* a){

  const OLattice< TVec >& x = a->x;
  const OLattice< TVec >& y = a->y;
  OLattice< TVec >& d = a->d;
  REAL* aptr = a->aptr;
  int Ns = a->Ns;
  const int* tab = a->tab;
  void (*func)(REAL*, REAL*, REAL*, REAL*, int) = a->func;
  
  for(int j=lo; j < hi; j++) { 
    int i=tab[j];

    REAL *xptr = (REAL *) &(x.elem(i).elem(0).elem(0).real());
    REAL *yptr = (REAL *) &(y.elem(i).elem(0).elem(0).real());
    REAL* zptr =  &(d.elem(i).elem(0).elem(0).real());
      
    // Get the no of 3vecs. s.start() and s.end() are inclusive so add +1
    func(zptr, aptr, xptr, yptr, Ns);
  }
}


//
// for vadd3_g5: add_g5ProjPlus, add_g5ProjMinus, sub_g5ProjPlus, sub_g5ProjMinus
//

// structure for vadd3_g5  of having order
struct ordered_vadd3_g5_user_arg{
  REAL* Out;
  REAL* X;
  REAL* Y;
  void (*func)(REAL*, REAL*, REAL*, int);
};

// user func for vadd3_g5 of having order
inline
void ordered_vadd3_g5_evaluate_function (int lo, int hi, int myId, ordered_vadd3_g5_user_arg* a){

  REAL* Out = a->Out;
  REAL* X = a->X;
  REAL* Y = a->Y;
  void (*func)(REAL*, REAL*, REAL*,int) = a->func;

  int n_4vec = hi - lo;

  int index = lo * 24;
  X = &X[index];
  Y = &Y[index];
  Out = &Out[index];

  func(Out, X, Y, n_4vec);

}

// structure for vadd3_g5 of NOT having order
struct unordered_vadd3_g5_user_arg{
  unordered_vadd3_g5_user_arg(  const OLattice< TVec >& x_,
				OLattice< TVec >& d_,
				int Ns_,
				const int* tab_,
				void (*func_)(REAL*, REAL*, REAL*, int)) : x(x_),d(d_),Ns(Ns_),tab(tab_),func(func_) {}

  const OLattice< TVec >& x;
  OLattice< TVec >& d;
  int Ns;
  const int* tab;
  void (*func)(REAL*, REAL*, REAL*, int);
};

// user func for vadd3_g5 of NOT having order
inline
void unordered_vadd3_g5_evaluate_function (int lo, int hi, int myId, unordered_vadd3_g5_user_arg* a){

  const OLattice< TVec >& x = a->x;
  OLattice< TVec >& d = a->d;
  int Ns = a->Ns;
  const int* tab = a->tab;
  void (*func)(REAL*, REAL*, REAL*, int) = a->func;
  
  for(int j=lo; j < hi; j++) { 
    int i=tab[j];
    REAL* xptr = (REAL *)&(x.elem(i).elem(0).elem(0).real());
    REAL* yptr = (REAL *)&(d.elem(i).elem(0).elem(0).real());
    func(yptr, yptr, xptr, Ns);
  }
}



//
// for vaxpy3_g5: axpyz_g5ProjPlus, axpyz_g5ProjMinus, axmyz_g5ProjPlus, axmyz_g5ProjMinus
//

// structure for vaxpy3_g5  of having order
struct ordered_vaxpy3_g5_user_arg{
  REAL* Out;
  REAL* scalep;
  REAL* InScale;
  REAL* Add;
  void (*func)(REAL*, REAL*, REAL*, REAL*, int);
};

// user func for vaxpy3_g5 of having order
inline
void ordered_vaxpy3_g5_evaluate_function (int lo, int hi, int myId, ordered_vaxpy3_g5_user_arg* a){

  REAL* Out = a->Out;
  REAL* scalep = a->scalep;
  REAL* InScale = a->InScale;
  REAL* Add = a->Add;
  void (*func)(REAL*, REAL*, REAL*, REAL*, int) = a->func;

  int n_4vec = hi - lo;

  int index = lo * 24;
  Out = &Out[index];
  InScale = &InScale[index];
  Add = &Add[index];

  func(Out, scalep, InScale, Add, n_4vec);

}

// structure for vaxpy3_g5 of NOT having order
struct unordered_vaxpy3_g5_user_arg{
  unordered_vaxpy3_g5_user_arg(
			       const OLattice< TVec >& x_,
			       const OLattice< TVec >& y_,
			       OLattice< TVec >& d_,
			       REAL* aptr_,
			       int Ns_,
			       const int* tab_,
			       void (*func_)(REAL*, REAL*, REAL*, REAL*, int)) : x(x_), y(y_),d(d_),aptr(aptr_),Ns(Ns_), tab(tab_), func(func_) {}
  
  const OLattice< TVec >& x;
  const OLattice< TVec >& y;
  OLattice< TVec >& d;
  REAL* aptr;
  int Ns;
  const int* tab;
  void (*func)(REAL*, REAL*, REAL*, REAL*, int);
};

// user func for vaxpy3_g5 of NOT having order
inline
void unordered_vaxpy3_g5_evaluate_function (int lo, int hi, int myId, unordered_vaxpy3_g5_user_arg* a){

  const OLattice< TVec >& x = a->x;
  const OLattice< TVec >& y = a->y;
  OLattice< TVec >& d = a->d;
  REAL* aptr = a->aptr;
  int Ns = a->Ns;
  const int* tab = a->tab;
  void (*func)(REAL*, REAL*, REAL*, REAL*, int) = a->func;
  
  for(int j=lo; j < hi; j++) { 
    int i=tab[j];
    REAL *xptr = (REAL *) &(x.elem(i).elem(0).elem(0).real());
    REAL *yptr = (REAL *) &(y.elem(i).elem(0).elem(0).real());
    REAL* zptr =  &(d.elem(i).elem(0).elem(0).real());
      
    func(zptr, aptr, xptr, yptr, Ns);
  }
}




//
// for vscal3_g5: scal_g5ProjPlus, scal_g5ProjMinus
//

// structure for vscal_g5  of having order
struct ordered_vscal_g5_user_arg{
  REAL* Out;
  REAL* scalep;
  REAL* In;
  void (*func)(REAL*, REAL*, REAL*, int);
};

// user func for vscal_g5 of having order
inline
void ordered_vscal_g5_evaluate_function (int lo, int hi, int myId, ordered_vscal_g5_user_arg* a){

  REAL* Out = a->Out;
  REAL* scalep = a->scalep;
  REAL* In = a->In;
  void (*func)(REAL*, REAL*, REAL*, int) = a->func;

  int n_4vec = hi - lo;

  int index = lo * 24;
  Out = &Out[index];
  In = &In[index];
 
  func(Out, scalep, In, n_4vec);

}

// structure for vscal_g5 of NOT having order
struct unordered_vscal_g5_user_arg{
  unordered_vscal_g5_user_arg(const OLattice< TVec >& x_,
			      OLattice< TVec >& d_,
			      REAL* aptr_,
			      int Ns_,
			      const int* tab_,
			      void (*func_)(REAL*, REAL*, REAL*, int))
  : x(x_), d(d_),aptr(aptr_),Ns(Ns_),tab(tab_),func(func_) {}

  const OLattice< TVec >& x;
  OLattice< TVec >& d;
  REAL* aptr;
  int Ns;
  const int* tab;
  void (*func)(REAL*, REAL*, REAL*, int);
};

// user func for vscal_g5 of NOT having order
inline
void unordered_vscal_g5_evaluate_function (int lo, int hi, int myId, unordered_vscal_g5_user_arg* a){

  const OLattice< TVec >& x = a->x;
  OLattice< TVec >& d = a->d;
  REAL* aptr = a->aptr;
  int Ns = a->Ns;
  const int* tab = a->tab;
  void (*func)(REAL*, REAL*, REAL*, int) = a->func;
  
  for(int j=lo; j < hi; j++) { 
    int i=tab[j];
    REAL *xptr = (REAL *) &(x.elem(i).elem(0).elem(0).real());
    REAL* zptr =  &(d.elem(i).elem(0).elem(0).real());
      
    func(zptr, aptr, xptr, Ns);
  }
}





//
// for vaxpby3_g5: axpbyz_g5ProjPlus, axpbyz_g5ProjMinus, axmbyz_g5ProjPlus, axmbyz_g5ProjMinus
//

// structure for vaxpby3_g5  of having order
struct ordered_vaxpby3_g5_user_arg{
  REAL* Out;
  REAL* scalep;
  REAL* InScale;
  REAL* scalep2;
  REAL* Add;
  void (*func)(REAL*, REAL*, REAL*, REAL*, REAL*, int);
};

// user func for vaxpby3_g5 of having order
inline
void ordered_vaxpby3_g5_evaluate_function (int lo, int hi, int myId, ordered_vaxpby3_g5_user_arg* a){

  REAL* Out = a->Out;
  REAL* scalep = a->scalep;
  REAL* InScale = a->InScale;
  REAL* scalep2 = a->scalep2;
  REAL* Add = a->Add;
  void (*func)(REAL*, REAL*, REAL*, REAL*, REAL*, int) = a->func;

  int n_4vec = hi - lo;

  int index = lo * 24;
  Out = &Out[index];
  InScale = &InScale[index];
  Add = &Add[index];

  func(Out, scalep, InScale, scalep2, Add, n_4vec);

}

// structure for vaxpby3_g5 of NOT having order
struct unordered_vaxpby3_g5_user_arg{
  unordered_vaxpby3_g5_user_arg(
				const OLattice< TVec >& x_,
				const OLattice< TVec >& y_,
				OLattice< TVec >& d_,
				REAL* aptr_,
				REAL* bptr_,
				int Ns_,
				const int* tab_,
				void (*func_)(REAL*, REAL*, REAL*, REAL*, REAL*, int)) : x(x_),y(y_),d(d_),aptr(aptr_),bptr(bptr_),Ns(Ns_),tab(tab_),func(func_) {}

  const OLattice< TVec >& x;
  const OLattice< TVec >& y;
  OLattice< TVec >& d;
  REAL* aptr; 
  REAL* bptr;
  int Ns;
  const int* tab;
  void (*func)(REAL*, REAL*, REAL*, REAL*, REAL*, int);
};

// user func for vaxpby3_g5 of NOT having order
inline
void unordered_vaxpby3_g5_evaluate_function (int lo, int hi, int myId, unordered_vaxpby3_g5_user_arg* a){

  const OLattice< TVec >& x = a->x;
  const OLattice< TVec >& y = a->y;
  OLattice< TVec >& d = a->d;
  REAL* aptr = a->aptr;
  REAL* bptr = a->bptr;
  int Ns = a->Ns;
  const int* tab = a->tab;
  void (*func)(REAL*, REAL*, REAL*, REAL*, REAL*, int) = a->func;
  
  for(int j=lo; j < hi; j++) { 
    int i=tab[j];
    REAL *xptr = (REAL *) &(x.elem(i).elem(0).elem(0).real());
    REAL *yptr = (REAL *) &(y.elem(i).elem(0).elem(0).real());
    REAL* zptr =  &(d.elem(i).elem(0).elem(0).real());
      
    func(zptr, aptr, xptr, bptr, yptr, Ns);
  }
}





//
// for scal3_g5
//

// structure for scal_g5  of having order
struct ordered_scal_g5_user_arg{
  REAL* Out;
  REAL* scalep;
  REAL* In;
};

// user func for scal_g5 of having order
inline
void ordered_scal_g5_evaluate_function (int lo, int hi, int myId, ordered_scal_g5_user_arg* a){

  REAL* Out = a->Out;
  REAL* scalep = a->scalep;
  REAL* In = a->In;

  int n_4vec = hi - lo;

  int index = lo * 24;
  Out = &Out[index];
  In = &In[index];
 
  scal_g5(Out, scalep, In, n_4vec);

}

// structure for scal_g5 of NOT having order
struct unordered_scal_g5_user_arg{
  unordered_scal_g5_user_arg(  const OLattice< TVec >& x_,
			       OLattice< TVec >& d_,
			       REAL* aptr_,
			       int Ns_,
			       const int* tab_):x(x_),d(d_),aptr(aptr_),Ns(Ns_),tab(tab_) {}

  const OLattice< TVec >& x;
  OLattice< TVec >& d;
  REAL* aptr;
  int Ns;
  const int* tab;
};

// user func for scal_g5 of NOT having order
inline
void unordered_scal_g5_evaluate_function (int lo, int hi, int myId, unordered_scal_g5_user_arg* a){

  const OLattice< TVec >& x = a->x;
  OLattice< TVec >& d = a->d;
  REAL* aptr = a->aptr;
  int Ns = a->Ns;
  const int* tab = a->tab;
 
  for(int j=lo; j < hi; j++) { 
    int i=tab[j];
    REAL *xptr = (REAL *) &(x.elem(i).elem(0).elem(0).real());
    REAL* zptr =  &(d.elem(i).elem(0).elem(0).real());
      
    scal_g5(zptr, aptr, xptr, Ns);
  }
}





//
// for xOpayz_g5: xpayz_g5, xmayz_g5, xpayz_ig5, xmayz_ig5
//

// structure for x0payz_g5  of having order
struct ordered_xOpayz_g5_user_arg{
  REAL* Out;
  REAL* scalep;
  REAL* InScale;
  REAL* Add;
  void (*func)(REAL*, REAL*, REAL*, REAL*, int);
};

// user func for xOpayz_g5 of having order
inline
void ordered_xOpayz_g5_evaluate_function (int lo, int hi, int myId, ordered_xOpayz_g5_user_arg* a){

  REAL* Out = a->Out;
  REAL* scalep = a->scalep;
  REAL* InScale = a->InScale;
  REAL* Add = a->Add;
  void (*func)(REAL*, REAL*, REAL*, REAL*, int) = a->func;

  int n_4vec = hi - lo;

  int index = lo * 24;
  Out = &Out[index];
  InScale = &InScale[index];
  Add = &Add[index];

  func(Out, scalep, InScale, Add, n_4vec);

}

// structure for xOpayz_g5 of NOT having order
struct unordered_xOpayz_g5_user_arg{
  unordered_xOpayz_g5_user_arg(
			       const OLattice< TVec >& x_,
			       const OLattice< TVec >& y_,
			       OLattice< TVec >& d_,
			       REAL* aptr_,
			       int Ns_,
			       const int* tab_,
			       void (*func_)(REAL*, REAL*, REAL*, REAL*, int)) :
  x(x_),y(y_),d(d_),aptr(aptr_),Ns(Ns_),tab(tab_),func(func_) {}

  const OLattice< TVec >& x;
  const OLattice< TVec >& y;
  OLattice< TVec >& d;
  REAL* aptr;
  int Ns;
  const int* tab;
  void (*func)(REAL*, REAL*, REAL*, REAL*, int);
};

// user func for xOpayz_g5 of NOT having order
inline
void unordered_xOpayz_g5_evaluate_function (int lo, int hi, int myId, unordered_xOpayz_g5_user_arg* a){

  const OLattice< TVec >& x = a->x;
  const OLattice< TVec >& y = a->y;
  OLattice< TVec >& d = a->d;
  REAL* aptr = a->aptr;
  int Ns = a->Ns;
  const int* tab = a->tab;
  void (*func)(REAL*, REAL*, REAL*, REAL*, int) = a->func;
  
  for(int j=lo; j < hi; j++) { 
    int i=tab[j];
    REAL *xptr = (REAL *) &(x.elem(i).elem(0).elem(0).real());
    REAL *yptr = (REAL *) &(y.elem(i).elem(0).elem(0).real());
    REAL* zptr =  &(d.elem(i).elem(0).elem(0).real());
      
    func(zptr, aptr, xptr, yptr, Ns);
  }
}




//
// for axOpbyz_g5: axpbyz_g5, g5_axmbyz,  axpbyz_ig5,  axmbyz_ig5
//

// structure for axOpbyz_g5  of having order
struct ordered_axOpbyz_g5_user_arg{
  REAL* Out;
  REAL* scalep;
  REAL* InScale;
  REAL* scalep2;
  REAL* Add;
  void (*func)(REAL*, REAL*, REAL*, REAL*, REAL*, int);
};

// user func for axOpbyz_g5 of having order
inline
void ordered_axOpbyz_g5_evaluate_function (int lo, int hi, int myId, ordered_axOpbyz_g5_user_arg* a){

  REAL* Out = a->Out;
  REAL* scalep = a->scalep;
  REAL* InScale = a->InScale;
  REAL* scalep2 = a->scalep2;
  REAL* Add = a->Add;
  void (*func)(REAL*, REAL*, REAL*, REAL*, REAL*, int) = a->func;

  int n_4vec = hi - lo;

  int index = lo * 24;
  Out = &Out[index];
  InScale = &InScale[index];
  Add = &Add[index];

  func(Out, scalep, InScale, scalep2, Add, n_4vec);

}

// structure for axOpbyz_g5 of NOT having order
struct unordered_axOpbyz_g5_user_arg{
  unordered_axOpbyz_g5_user_arg(  const OLattice< TVec >& x_,
				  const OLattice< TVec >& y_,
				  OLattice< TVec >& d_,
				  REAL* aptr_,
				  REAL* bptr_,
				  int Ns_,
				  const int* tab_,
				  void (*func_)(REAL*, REAL*, REAL*, REAL*, REAL*, int)) : x(x_),y(y_),d(d_),aptr(aptr_),bptr(bptr_),Ns(Ns_),tab(tab_),func(func_) {}

  const OLattice< TVec >& x;
  const OLattice< TVec >& y;
  OLattice< TVec >& d;
  REAL* aptr; 
  REAL* bptr;
  int Ns;
  const int* tab;
  void (*func)(REAL*, REAL*, REAL*, REAL*, REAL*, int);
};

// user func for axOpbyz_g5 of NOT having order
inline
void unordered_axOpbyz_g5_evaluate_function (int lo, int hi, int myId, unordered_axOpbyz_g5_user_arg* a){

  const OLattice< TVec >& x = a->x;
  const OLattice< TVec >& y = a->y;
  OLattice< TVec >& d = a->d;
  REAL* aptr = a->aptr;
  REAL* bptr = a->bptr;
  int Ns = a->Ns;
  const int* tab = a->tab;
  void (*func)(REAL*, REAL*, REAL*, REAL*, REAL*, int) = a->func;
  
  for(int j=lo; j < hi; j++) { 
    int i=tab[j];
    REAL *xptr = (REAL *) &(x.elem(i).elem(0).elem(0).real());
    REAL *yptr = (REAL *) &(y.elem(i).elem(0).elem(0).real());
    REAL* zptr =  &(d.elem(i).elem(0).elem(0).real());
      
    func(zptr, aptr, xptr, bptr, yptr, Ns);
  }
}




//
// for xOpayz_ig5:  xpayz_ig5, xmayz_ig5
//

// structure for x0payz_ig5 of having order
struct ordered_xOpayz_ig5_user_arg{
  REAL* Out;
  REAL* scalep;
  REAL* InScale;
  REAL* Add;
  void (*func)(REAL*, REAL*, REAL*, REAL*, int);
};

// user func for xOpayz_ig5 of having order
inline
void ordered_xOpayz_ig5_evaluate_function (int lo, int hi, int myId, ordered_xOpayz_ig5_user_arg* a){

  REAL* Out = a->Out;
  REAL* scalep = a->scalep;
  REAL* InScale = a->InScale;
  REAL* Add = a->Add;
  void (*func)(REAL*, REAL*, REAL*, REAL*, int) = a->func;

  int n_4vec = hi - lo;

  int index = lo * 24;
  Out = &Out[index];
  InScale = &InScale[index];
  Add = &Add[index];

  func(Out, scalep, InScale, Add, n_4vec);

}

// structure for xOpayz_ig5 of NOT having order (with xyz)
struct unordered_xOpayz_ig5_y_user_arg{
  unordered_xOpayz_ig5_y_user_arg(
				  const OLattice< TVec >& x_,
				  const OLattice< TVec >& y_,
				  OLattice< TVec >& d_,
				  REAL* aptr_,
				  int Ns_,
				  const int* tab_,
				  void (*func_)(REAL*, REAL*, REAL*, REAL*, int)) : x(x_), y(y_),d(d_),aptr(aptr_),Ns(Ns_),tab(tab_),func(func_) {}

  const OLattice< TVec >& x;
  const OLattice< TVec >& y;
  OLattice< TVec >& d;
  REAL* aptr;
  int Ns;
  const int* tab;
  void (*func)(REAL*, REAL*, REAL*, REAL*, int);
};

// user func for xOpayz_ig5 of NOT having order (with xz)
inline
void unordered_xOpayz_ig5_y_evaluate_function (int lo, int hi, int myId, unordered_xOpayz_ig5_y_user_arg* a){

  const OLattice< TVec >& x = a->x;
  const OLattice< TVec >& y = a->y;
  OLattice< TVec >& d = a->d;
  REAL* aptr = a->aptr;
  int Ns = a->Ns;
  const int* tab = a->tab;
  void (*func)(REAL*, REAL*, REAL*, REAL*, int) = a->func;
  
  for(int j=lo; j < hi; j++) { 
    int i=tab[j];
    REAL *xptr = (REAL *) &(x.elem(i).elem(0).elem(0).real());
    REAL *yptr = (REAL *) &(y.elem(i).elem(0).elem(0).real());
    REAL* zptr =  &(d.elem(i).elem(0).elem(0).real());
      
    func(zptr, aptr, xptr, yptr, Ns);
  }
}




// structure for xOpayz_ig5 of NOT having order (with xz)
struct unordered_xOpayz_ig5_z_user_arg{
  unordered_xOpayz_ig5_z_user_arg(
				  const OLattice< TVec >& x_,
				  OLattice< TVec >& d_,
				  REAL* aptr_,
				  int Ns_,
				  const int* tab_,
				  void (*func_)(REAL*, REAL*, REAL*, REAL*, int)) : x(x_),d(d_),aptr(aptr_),Ns(Ns_),tab(tab_),func(func_) {}
  const OLattice< TVec >& x;
  OLattice< TVec >& d;
  REAL* aptr;
  int Ns;
  const int* tab;
  void (*func)(REAL*, REAL*, REAL*, REAL*, int);
};

// user func for xOpayz_ig5 of NOT having order (with xz)
inline
void unordered_xOpayz_ig5_z_evaluate_function (int lo, int hi, int myId, unordered_xOpayz_ig5_z_user_arg* a){

  const OLattice< TVec >& x = a->x;
  OLattice< TVec >& d = a->d;
  REAL* aptr = a->aptr;
  int Ns = a->Ns;
  const int* tab = a->tab;
  void (*func)(REAL*, REAL*, REAL*, REAL*, int) = a->func;
  
  for(int j=lo; j < hi; j++) { 
    int i=tab[j];
    REAL *xptr = (REAL *) &(x.elem(i).elem(0).elem(0).real());
    REAL* zptr =  &(d.elem(i).elem(0).elem(0).real());
      
    func(zptr, aptr, zptr, xptr, Ns);
  }
}






#endif
