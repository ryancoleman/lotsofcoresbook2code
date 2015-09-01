#ifndef QDP_SCALARSITE_SSE_BLAS_WRAPPER_H
#define QDP_SCALARSITE_SSE_BLAS_WRAPPER_H

////////////////////////////////
// Threading evaluates wrappers
//
// by Xu Guo, EPCC, 6 October, 2008
////////////////////////////////

//
// vaxOpy3: for vaxpy3 and vaxmy3
//

// structure for vaxOpy3 of having order
struct ordered_sse_vaxOpy3_user_arg{
  REAL32* Out;
  REAL32* scalep;
  REAL32* InScale;
  REAL32* Opt;
  void (*func)(REAL32*, REAL32*, REAL32*, REAL32*, int);
};

// user func for vaxOpy3 of having order
inline
void ordered_sse_vaxOpy3_evaluate_function (int lo, int hi, int myId, ordered_sse_vaxOpy3_user_arg* a){

  REAL32* Out = a->Out;
  REAL32* scalep = a->scalep;
  REAL32* InScale = a->InScale;
  REAL32* Opt = a->Opt;
  void (*func)(REAL32*, REAL32*, REAL32*, REAL32*, int) = a->func;

  int n_3vec = (hi - lo); // no of numbers

  int index = lo; // Rescale to real32

  InScale = &InScale[24*index];
  Opt = &Opt[24*index];
  Out = &Out[24*index];

  func(Out, scalep, InScale, Opt, n_3vec);
}

// structure for vaxOpy3 (with yptr only) of NOT having order
struct unordered_sse_vaxOpy3_y_user_arg{
  unordered_sse_vaxOpy3_y_user_arg(  const OLattice< TVec >& x_,
				     OLattice< TVec >& d_,
				     REAL32* scalep_,
				     const int* tab_,
				     int xy_order_,
				     void (*func_)(REAL32*, REAL32*, REAL32*, REAL32*, int)) 
  :x(x_),d(d_),scalep(scalep_),tab(tab_),xy_order(xy_order_),func(func_) {}

  const OLattice< TVec >& x;
  OLattice< TVec >& d;
  REAL32* scalep;
  const int* tab;
  int xy_order;
  void (*func)(REAL32*, REAL32*, REAL32*, REAL32*, int);
};

// user func for vaxOpy3 (with yptr only) of NOT having order
inline
void unordered_sse_vaxOpy3_y_evaluate_function (int lo, int hi, int myId, unordered_sse_vaxOpy3_y_user_arg* a){

  const OLattice< TVec >& x = a->x;
  OLattice< TVec >& d = a->d;
  REAL32* scalep = a->scalep;
  const int* tab = a->tab;
  int xy_order = a->xy_order;
  void (*func)(REAL32*, REAL32*, REAL32*, REAL32*, int) = a->func;

  if (xy_order){
     for(int j=lo; j < hi; j++) { 
       int i=tab[j];
       REAL32* xptr = (REAL32 *)&(x.elem(i).elem(0).elem(0).real());
       REAL32* yptr = &(d.elem(i).elem(0).elem(0).real());
       func(yptr, scalep, xptr, yptr, 1);
     }
  }
  else {
    for(int j=lo; j < hi; j++) { 
     int i=tab[j];
     REAL32* xptr = (REAL32 *)&(x.elem(i).elem(0).elem(0).real());
     REAL32* yptr = &(d.elem(i).elem(0).elem(0).real());
     func(yptr, scalep, yptr, xptr, 1);
   }
  }
 

}

// structure for vaxOpy3 (with zptr) of NOT having order
struct unordered_sse_vaxOpy3_z_user_arg{
  unordered_sse_vaxOpy3_z_user_arg(
				   const OLattice< TVec >& x_,
				   const OLattice< TVec >& y_,
				   OLattice< TVec >& d_,
				   REAL32* scalep_,
				   const int* tab_,
				   void (*func_)(REAL32*, REAL32*, REAL32*, REAL32*, int)) 
  : x(x_),y(y_),d(d_),scalep(scalep_),tab(tab_),func(func_) {}


  const OLattice< TVec >& x;
  const OLattice< TVec >& y;
  OLattice< TVec >& d;
  REAL32* scalep;
  const int* tab;
  void (*func)(REAL32*, REAL32*, REAL32*, REAL32*, int);
};

// user func for vaxpy3 (with zptr only) of NOT having order
inline
void unordered_sse_vaxOpy3_z_evaluate_function (int lo, int hi, int myId, unordered_sse_vaxOpy3_z_user_arg* a){

  const OLattice< TVec >& x = a->x;
  const OLattice< TVec >& y = a->y;
  OLattice< TVec >& d = a->d;
  REAL32* scalep = a->scalep;
  const int* tab = a->tab;
  void (*func)(REAL32*, REAL32*, REAL32*, REAL32*, int) = a->func;

  for(int j=lo; j < hi; j++) { 
     int i=tab[j];
     REAL32* xptr = (REAL32 *)&(x.elem(i).elem(0).elem(0).real());
     REAL32* yptr = (REAL32 *)&(y.elem(i).elem(0).elem(0).real());
     REAL32* zptr =  &(d.elem(i).elem(0).elem(0).real());
     func(zptr, scalep, xptr, yptr, 1);
   }

}




//
// vOp: for vadd and vsub
//

// structure for vOp of having order
struct ordered_sse_vOp_user_arg{
  REAL32* Out;
  REAL32* In1;
  REAL32* In2;
  void (*func)(REAL32*, REAL32*, REAL32*, int);
};

// user func for vOp of having order
inline
void ordered_sse_vOp_evaluate_function (int lo, int hi, int myId, ordered_sse_vOp_user_arg* a){

  REAL32* Out = a->Out;
  REAL32* In1 = a->In1;
  REAL32* In2 = a->In2;
  void (*func)(REAL32*, REAL32*, REAL32*, int) = a->func;

  int n_3vec = hi - lo;

  int index = lo;
  In1 = &In1[24*index];
  In2 = &In2[24*index];
  Out = &Out[24*index];

  func(Out, In1, In2, n_3vec);
}

// structure for vOp (with yptr only) of NOT having order
struct unordered_sse_vOp_y_user_arg{
  unordered_sse_vOp_y_user_arg(
			       const OLattice< TVec >& x_,
			       OLattice< TVec >& d_,
			       const int* tab_,
			       void (*func_)(REAL32*, REAL32*, REAL32*, int)) 
  : x(x_),d(d_),tab(tab_),func(func_) {}

  const OLattice< TVec >& x;
  OLattice< TVec >& d;
  const int* tab;
  void (*func)(REAL32*, REAL32*, REAL32*, int);
};

// user func for vOp (with yptr only) of NOT having order
inline
void unordered_sse_vOp_y_evaluate_function (int lo, int hi, int myId, unordered_sse_vOp_y_user_arg* a){

  const OLattice< TVec >& x = a->x;
  OLattice< TVec >& d = a->d;
  const int* tab = a->tab;
  void (*func)(REAL32*, REAL32*, REAL32*, int) = a->func;


  for(int j=lo; j < hi; j++) { 
      int i = tab[j];
      REAL32 *xptr = (REAL32 *)(&x.elem(i).elem(0).elem(0).real());
      REAL32 *yptr = (REAL32 *)(&d.elem(i).elem(0).elem(0).real());
      func(yptr, yptr, xptr,1);
  }

}

// structure for vOp (with zptr) of NOT having order
struct unordered_sse_vOp_z_user_arg{
  unordered_sse_vOp_z_user_arg(
			       const OLattice< TVec >& x_,
			       const OLattice< TVec >& y_,
			       OLattice< TVec >& d_,
			       const int* tab_,
			       void (*func_)(REAL32*, REAL32*, REAL32*, int))
  : x(x_),y(y_),d(d_),tab(tab_),func(func_) {}
  const OLattice< TVec >& x;
  const OLattice< TVec >& y;
  OLattice< TVec >& d;
  const int* tab;
  void (*func)(REAL32*, REAL32*, REAL32*, int);
};

// user func for vOp (with zptr only) of NOT having order
inline
void unordered_vOp_z_evaluate_function (int lo, int hi, int myId, unordered_sse_vOp_z_user_arg* a){

  const OLattice< TVec >& x = a->x;
  const OLattice< TVec >& y = a->y;
  OLattice< TVec >& d = a->d;
  const int* tab = a->tab;
  void (*func)(REAL32*, REAL32*, REAL32*, int) = a->func;

  for(int j=lo; j < hi; j++) { 
      int i = tab[j];
      REAL32 *xptr = (REAL32 *) &(x.elem(i).elem(0).elem(0).real());
      REAL32 *yptr = (REAL32 *) &(y.elem(i).elem(0).elem(0).real());
      REAL32* zptr =  &(d.elem(i).elem(0).elem(0).real());
            
      // Get the no of 3vecs. s.start() and s.end() are inclusive so add +1
      func(zptr, xptr, yptr, 1);

   }

}



//
// vscal
//

// structure for vscal of having order
struct ordered_sse_vscal_user_arg{
  REAL32* Out;
  REAL32* scalep;
  REAL32* In;
};

// user func for vscal of having order
inline
void ordered_sse_vscal_evaluate_function (int lo, int hi, int myId, ordered_sse_vscal_user_arg* a){

  REAL32* Out = a->Out;
  REAL32* scalep = a->scalep;
  REAL32* In = a->In;

  int n_3vec = hi - lo;
  int index = lo;
  Out = &Out[24*index];
  In = &In[24*index];
  vscal(Out, scalep, In, n_3vec);
}

// structure for vscal of NOT having order
struct unordered_sse_vscal_user_arg{
  unordered_sse_vscal_user_arg(
			       const OLattice< TVec >& x_,
			       OLattice< TVec >& d_,
			       REAL32* scalep_,
			       const int* tab_) : x(x_),d(d_),scalep(scalep_),tab(tab_) {}
  const OLattice< TVec >& x;
  OLattice< TVec >& d;
  REAL32* scalep;
  const int* tab;
};

// user func for vscal of NOT having order
inline
void unordered_sse_vscal_evaluate_function (int lo, int hi, int myId, unordered_sse_vscal_user_arg* a){

  const OLattice< TVec >& x = a->x;
  OLattice< TVec >& d = a->d;
  REAL32* scalep = a->scalep;
  const int* tab = a->tab;

  for(int j=lo; j < hi; j++) { 
     int i=tab[j];
     REAL32* xptr = (REAL32 *) &(x.elem(i).elem(0).elem(0).real());
     REAL32* zptr =  &(d.elem(i).elem(0).elem(0).real()); 
     vscal(zptr, scalep, xptr, 1);
   }

}



//
// vaxOpby3: for vaxpby3 and vaxmby3
//

// structure for vaxOpby3 of having order
struct ordered_sse_vaxOpby3_user_arg{
  REAL32* zptr;
  REAL32* aptr;
  REAL32* xptr;
  REAL32* bptr;
  REAL32* yptr;
  void (*func)(REAL32*, REAL32*, REAL32*, REAL32*, REAL32*, int);
};

// user func for vaxOpby3 of having order
inline
void ordered_sse_vaxOpby3_evaluate_function (int lo, int hi, int myId, ordered_sse_vaxOpby3_user_arg* arg){

  REAL32* zptr = arg->zptr;
  REAL32* aptr = arg->aptr;
  REAL32* xptr = arg->xptr;
  REAL32* bptr = arg->bptr;
  REAL32* yptr = arg->yptr;
  void (*func)(REAL32*, REAL32*, REAL32*, REAL32*, REAL32*, int) = arg->func;

  int n_3vec = (hi - lo);

  int index = 24*lo;
  zptr = &zptr[index];
  xptr = &xptr[index];
  yptr = &yptr[index];

  func(zptr, aptr, xptr, bptr, yptr, n_3vec);
}

// structure for vaxOpby3 of NOT having order
struct unordered_sse_vaxOpby3_user_arg{
  unordered_sse_vaxOpby3_user_arg(
				  REAL32* aptr_,
				  const OLattice< TVec >& x_,
				  REAL32* bptr_,
				  const OLattice< TVec >& y_,
				  OLattice< TVec >& d_,
				  const int* tab_,
				  void (*func_)(REAL32*, REAL32*, REAL32*, REAL32*, REAL32*, int)) 
  : aptr(aptr_),x(x_),bptr(bptr_),y(y_),d(d_),tab(tab_),func(func_) {}


  REAL32* aptr;
  const OLattice< TVec >& x;
  REAL32* bptr;
  const OLattice< TVec >& y;
  OLattice< TVec >& d;
  const int* tab;
  void (*func)(REAL32*, REAL32*, REAL32*, REAL32*, REAL32*, int);
};

// user func for vaxOpby3 of NOT having order
inline
void unordered_sse_vaxOpby3_evaluate_function (int lo, int hi, int myId, unordered_sse_vaxOpby3_user_arg* arg){
  REAL32* aptr = arg->aptr;
  const OLattice< TVec >& x = arg->x;
  REAL32* bptr = arg->bptr;
  const OLattice< TVec >& y = arg->y;
  OLattice< TVec >& d = arg->d;
  const int* tab = arg->tab;
  void (*func)(REAL32*, REAL32*, REAL32*, REAL32*, REAL32*, int) = arg->func;


  for(int j=lo; j < hi; j++) { 
      int i = tab[j];

      REAL32 *xptr = (REAL32 *) &(x.elem(i).elem(0).elem(0).real());
      REAL32 *yptr = (REAL32 *) &(y.elem(i).elem(0).elem(0).real());
      REAL32 * zptr =  &(d.elem(i).elem(0).elem(0).real());
      
      func(zptr, aptr, xptr, bptr, yptr, 1);
    
  }

}


struct ordered_sse_norm_single_user_arg{
  REAL32* vptr;    // Single Prec vector
  REAL64* results; // Always store sums in doubles
  void (*func)(REAL64*, REAL32*, int);
};

inline void ordered_norm_single_func(int lo, int hi, int myId, ordered_sse_norm_single_user_arg* a) 
  {

    int nvec = hi - lo;
    int index = 24*lo;
    REAL32* vptr = &(a->vptr[index]);
    void (*func)(REAL64*, REAL32*, int) = a->func;
    func( &(a->results[myId]), vptr, nvec);
  }

#endif
