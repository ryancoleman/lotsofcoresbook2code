#ifndef QDP_SCALARSITE_GENERIC_BLAS_DOUBLE_WRAPPER_H
#define QDP_SCALARSITE_GENERIC_BLAS_DOUBLE_WRAPPER_H



////////////////////////////////
// Threading evaluates wrappers
//
// by Xu Guo, EPCC, 6 October, 2008
////////////////////////////////

//
// vaxOpy4_double: for vaxpy4 and vaxmy4 (double precision
//

// structure for vaxOpy4 of having order
struct ordered_sse_vaxOpy4_double_user_arg{
  REAL64* Out;
  REAL64* scalep;
  REAL64* InScale;
  void (*func)(REAL64*, REAL64*, REAL64*, int);
};

// user func for vaxOpy4 of having order
inline
void ordered_sse_vaxOpy4_double_evaluate_function (int lo, int hi, int myId, ordered_sse_vaxOpy4_double_user_arg* a){

  REAL64* Out = a->Out;
  REAL64* scalep = a->scalep;
  REAL64* InScale = a->InScale;
  void (*func)(REAL64*, REAL64*, REAL64*, int) = a->func;

  int n_4vec = hi - lo;

  int index = lo * 24;
  
  InScale = &InScale[index];
  Out = &Out[index];

  func(Out, scalep, InScale, n_4vec);
}


// structure for vaxOpy4_double of NOT having order
struct unordered_sse_vaxOpy4_double_user_arg{
  unordered_sse_vaxOpy4_double_user_arg(  const OLattice< DVec >& x_,
					  OLattice< DVec >& d_,
					  REAL64* scalep_,
					  int Ns_,
					  const int* tab_,
					  void (*func_)(REAL64*, REAL64*, REAL64*, int)) : 
  x(x_),d(d_),scalep(scalep_),Ns(Ns_),tab(tab_),func(func_) {}

  const OLattice< DVec >& x;
  OLattice< DVec >& d;
  REAL64* scalep;
  int Ns;
  const int* tab;
  void (*func)(REAL64*, REAL64*, REAL64*, int);
};

// user func for vaxOpy4_double of NOT having order
inline
void unordered_sse_vaxOpy4_double_evaluate_function (int lo, int hi, int myId, unordered_sse_vaxOpy4_double_user_arg* a){

  const OLattice< DVec >& x = a->x;
  OLattice< DVec >& d = a->d;
  REAL64* scalep = a->scalep;
  int Ns = a->Ns;
  const int* tab = a->tab;
  void (*func)(REAL64*, REAL64*, REAL64*, int) = a->func;

   for(int j=lo; j < hi; j++) { 
     int i=tab[j];
     REAL64* xptr = (REAL64 *)&(x.elem(i).elem(0).elem(0).real());
     REAL64* yptr = &(d.elem(i).elem(0).elem(0).real());
     func(yptr, scalep, xptr, Ns);
  }

}



//
// vaxOpyz4_double: for vaxpyz4 and vaxmyz4 (double precision)
//

// structure for vaxOpyz4 of having order
struct ordered_sse_vaxOpyz4_double_user_arg{
  REAL64* Out;
  REAL64* scalep;
  REAL64* InScale;
  REAL64* Opt;
  void (*func)(REAL64*, REAL64*, REAL64*, REAL64*, int);
};

// user func for vaxOpyz4 of having order
inline
void ordered_sse_vaxOpyz4_double_evaluate_function (int lo, int hi, int myId, ordered_sse_vaxOpyz4_double_user_arg* a){

  REAL64* Out = a->Out;
  REAL64* scalep = a->scalep;
  REAL64* InScale = a->InScale;
  REAL64* Opt = a->Opt;
  void (*func)(REAL64*, REAL64*, REAL64*, REAL64*, int) = a->func;

  int n_4vec = hi - lo;

  int index = lo * 24;

  InScale = &InScale[index];
  Out = &Out[index];
  Opt = &Opt[index];

  func(Out, scalep, InScale, Opt, n_4vec);
}


// structure for vaxOpyz4_double of NOT having order
struct unordered_sse_vaxOpyz4_double_user_arg{
  unordered_sse_vaxOpyz4_double_user_arg(
					 const OLattice< DVec >& x_,
					 const OLattice< DVec >& y_,
					 OLattice< DVec >& d_,
					 REAL64* scalep_,
					 int Ns_,
					 const int* tab_,
					 void (*func_)(REAL64*, REAL64*, REAL64*, REAL64*, int))
  :x(x_),y(y_),d(d_),scalep(scalep_),Ns(Ns_),tab(tab_),func(func_) {}
  const OLattice< DVec >& x;
  const OLattice< DVec >& y;
  OLattice< DVec >& d;
  REAL64* scalep;
  int Ns;
  const int* tab;
  void (*func)(REAL64*, REAL64*, REAL64*, REAL64*, int);
};

// user func for vaxOpyz4_double of NOT having order
inline
void unordered_sse_vaxOpyz4_double_evaluate_function (int lo, int hi, int myId, unordered_sse_vaxOpyz4_double_user_arg* a){

  const OLattice< DVec >& x = a->x;
  const OLattice< DVec >& y = a->y;
  OLattice< DVec >& d = a->d;
  REAL64* scalep = a->scalep;
  int Ns = a->Ns;
  const int* tab = a->tab;
  void (*func)(REAL64*, REAL64*, REAL64*, REAL64*, int) = a->func;

 for(int j=lo; j < hi; j++) { 
      int i=tab[j];
      REAL64* xptr = (REAL64 *)&(x.elem(i).elem(0).elem(0).real());
      REAL64* yptr = (REAL64 *)&(y.elem(i).elem(0).elem(0).real());
      REAL64* zptr =  &(d.elem(i).elem(0).elem(0).real());

      func(zptr, scalep, xptr, yptr, Ns);
 }


}



//
// vscal4_double: for vscal4 (double precision)
//

// structure for vscal4 of having order
struct ordered_sse_vscal4_double_user_arg{
  REAL64* zptr;
  REAL64* aptr;
  REAL64* xptr;
};

// user func for vscal4 of having order
inline
void ordered_sse_vscal4_double_evaluate_function (int lo, int hi, int myId, ordered_sse_vscal4_double_user_arg* a){
  REAL64* zptr = a->zptr;
  REAL64* aptr = a->aptr;
  REAL64* xptr = a->xptr;

  int n_4vec = hi - lo;

  int index = lo * 24;

  zptr = &zptr[index];
  xptr = &xptr[index];

  vscal4(zptr,aptr,xptr,n_4vec);
}


// structure for vscal4_double of NOT having order
struct unordered_sse_vscal4_double_user_arg{
  unordered_sse_vscal4_double_user_arg(
				       const OLattice< DVec >& x_,
				       OLattice< DVec >& d_,
				       REAL64* aptr_,
				       int Ns_,
				       const int* tab_)
  : x(x_),d(d_),aptr(aptr_),Ns(Ns_),tab(tab_) {}
  const OLattice< DVec >& x;
  OLattice< DVec >& d;
  REAL64* aptr;
  int Ns;
  const int* tab;
};

// user func for vscal4_double of NOT having order
inline
void unordered_sse_vscal4_double_evaluate_function (int lo, int hi, int myId, unordered_sse_vscal4_double_user_arg* a){

  const OLattice< DVec >& x = a->x;
  OLattice< DVec >& d = a->d;
  REAL64* aptr = a->aptr;
  int Ns = a->Ns;
  const int* tab = a->tab;

   for(int j=lo; j < hi; j++) { 
     int i=tab[j];
     REAL64 *xptr = (REAL64 *) &(x.elem(i).elem(0).elem(0).real());
     REAL64 *zptr =  &(d.elem(i).elem(0).elem(0).real());
 
     vscal4(zptr, aptr, xptr, Ns);
  }

}



//
// vaxOpby4_double: for vaxpby4 and vaxmby4 (double precision)
//

// structure for vaxOpby4 of having order
struct ordered_sse_vaxOpby4_double_user_arg{
  REAL64* yptr;
  REAL64* aptr;
  REAL64* xptr;
  REAL64* bptr;
  void (*func)(REAL64*, REAL64*, REAL64*, REAL64*, int);
};

// user func for vaxOpby4 of having order
inline
void ordered_sse_vaxOpby4_double_evaluate_function (int lo, int hi, int myId, ordered_sse_vaxOpby4_double_user_arg* a){

  REAL64* yptr = a->yptr;
  REAL64* aptr = a->aptr;
  REAL64* xptr = a->xptr;
  REAL64* bptr = a->bptr;
  void (*func)(REAL64*, REAL64*, REAL64*, REAL64*, int) = a->func;

  int n_4vec = hi - lo;

  int index = lo * 24;

  yptr = &yptr[index];
  xptr = &xptr[index];

  func(yptr, aptr, xptr, bptr, n_4vec);
}



//
// vaxOpbyz4_double: for vaxpbyz4 and vaxmbyz4 (double precision)
//

// structure for vaxOpbyz4 of having order
struct ordered_sse_vaxOpbyz4_double_user_arg{
  REAL64* zptr;
  REAL64* aptr;
  REAL64* xptr;
  REAL64* bptr;
  REAL64* yptr;
  void (*func)(REAL64*, REAL64*, REAL64*, REAL64*, REAL64*, int);
};

// user func for vaxOpbyz4 of having order
inline
void ordered_sse_vaxOpbyz4_double_evaluate_function (int lo, int hi, int myId, ordered_sse_vaxOpbyz4_double_user_arg* a){

  REAL64* zptr = a->zptr;
  REAL64* aptr = a->aptr;
  REAL64* xptr = a->xptr;
  REAL64* bptr = a->bptr;
  REAL64* yptr = a->yptr;
  void (*func)(REAL64*, REAL64*, REAL64*, REAL64*, REAL64*, int) = a->func;
  
  int n_4vec = hi - lo;

  int index = lo * 24;

  zptr = &zptr[index];
  xptr = &xptr[index];
  yptr = &yptr[index];

  func(zptr,aptr,xptr,bptr,yptr,n_4vec);
}


// structure for vaxOpbyz4_double of NOT having order
struct unordered_sse_vaxOpbyz4_double_user_arg{
  unordered_sse_vaxOpbyz4_double_user_arg(  const OLattice< DVec >& x_,
					    const OLattice< DVec >& y_,
					    OLattice< DVec >& d_,
					    REAL64* aptr_,
					    REAL64* bptr_,
					    int Ns_,
					    const int* tab_,
					    void (*func_)(REAL64*, REAL64*, REAL64*, REAL64*, REAL64*, int))
  :x(x_),y(y_),d(d_),aptr(aptr_),bptr(bptr_),Ns(Ns_),tab(tab_),func(func_) {}


  const OLattice< DVec >& x;
  const OLattice< DVec >& y;
  OLattice< DVec >& d;
  REAL64* aptr;
  REAL64* bptr;
  int Ns;
  const int* tab;
  void (*func)(REAL64*, REAL64*, REAL64*, REAL64*, REAL64*, int);
};

// user func for vaxOpbyz4_double of NOT having order
inline
void unordered_sse_vaxOpbyz4_double_evaluate_function (int lo, int hi, int myId, unordered_sse_vaxOpbyz4_double_user_arg* a){

  const OLattice< DVec >& x = a->x;
  const OLattice< DVec >& y = a->y;
  OLattice< DVec >& d = a->d;
  REAL64* aptr = a->aptr;
  REAL64* bptr = a->bptr;
  int Ns = a->Ns;
  const int* tab = a->tab;
  void (*func)(REAL64*, REAL64*, REAL64*, REAL64*, REAL64*, int) = a->func;

  for(int j=lo; j < hi; j++) { 
      int i=tab[j];
   
      REAL64 *xptr = (REAL64 *) &(x.elem(i).elem(0).elem(0).real());
      REAL64 *yptr = (REAL64 *) &(y.elem(i).elem(0).elem(0).real());
      REAL64* zptr =  &(d.elem(i).elem(0).elem(0).real());
      
      func(zptr, aptr, xptr, bptr, yptr, Ns);
 }


  
}

struct ordered_norm_double_user_arg {
  REAL64* vptr;
  REAL64* results;
  void (*func)(REAL64*,  REAL64*, int);
};


inline void ordered_norm_double_func(int lo, int hi, int myId, ordered_norm_double_user_arg* a) 
  {

    int nvec = hi - lo;
    int index = lo*24;
    REAL64* vptr = &(a->vptr[index]);
    void (*func)(REAL64*, REAL64*, int) = a->func;
    func( &(a->results[myId]), vptr, nvec);

  }


struct ordered_inner_product_double_user_arg {
  REAL64* xptr;
  REAL64* yptr;
  REAL64* results;
  void (*func)(REAL64*,  REAL64*, REAL64*, int);
};


inline void ordered_inner_product_double_func(int lo, int hi, int myId, ordered_inner_product_double_user_arg* a) 
  {

    int nvec = hi - lo;
    int index = lo*24;
    REAL64* xptr = &(a->xptr[index]);
    REAL64* yptr = &(a->yptr[index]);

    void (*func)(REAL64*, REAL64*, REAL64*, int) = a->func;
    func( &(a->results[2*myId]), xptr, yptr, nvec);
  }



#endif
