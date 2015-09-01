#ifndef QDP_SCALARSITE_GENERIC_CBLAS_WRAPPER_H
#define QDP_SCALARSITE_GENERIC_CBLAS_WRAPPER_H

////////////////////////////////
// Threading evaluates wrappers
//
// by Xu Guo, EPCC, 26 August, 2008
////////////////////////////////

//
// vcscal
//

// structure for vcscal of having order
struct ordered_vcscal_user_arg{
  REAL* Out;
  REAL* scalep;
  REAL* In;
};

// user func for vcscal of having order
inline
void ordered_vcscal_evaluate_function (int lo, int hi, int myId, ordered_vcscal_user_arg* a){

  REAL* Out = a->Out;
  REAL* scalep = a->scalep;
  REAL* In = a->In;

  int n_3vec = hi - lo;
  int index = lo * 24;
  Out = &Out[index];
  In = &In[index];
  vcscal(Out, scalep, In, n_3vec);
}

// structure for vcscal of NOT having order
struct unordered_vcscal_user_arg{
  unordered_vcscal_user_arg(
			    const OLattice< CTVec > &x_,
			    OLattice< CTVec >& d_,
			    REAL* scalep_,
			    const int* tab_) : x(x_), d(d_), scalep(scalep_), tab(tab_) {}

  const OLattice< CTVec > &x;
  OLattice< CTVec >& d;
  REAL* scalep;
  const int* tab;
};

// user func for vcscal of NOT having order
inline
void unordered_vcscal_evaluate_function (int lo, int hi, int myId, unordered_vcscal_user_arg* a){
  const OLattice< CTVec > &x = a->x;
  OLattice< CTVec >& d = a->d;
  REAL* scalep = a->scalep;
  const int* tab = a->tab;

  for(int j=lo; j < hi; j++) { 
     int i=tab[j];
     REAL *d_start = &(d.elem(i).elem(0).elem(0).real());
     REAL *x_start = (REAL *) &(x.elem(i).elem(0).elem(0).real());
     vcscal(d_start, scalep, x_start, 1);
   }

}




//
// vcaxpy3
//

// structure for vcaxpy3 of having order
struct ordered_vcaxpy3_user_arg{
  REAL* Out;
  REAL* scalep;
  REAL* InScale;
  REAL* Add;
};

// user func for vcaxpy3 of having order
inline
void ordered_vcaxpy3_evaluate_function (int lo, int hi, int myId, ordered_vcaxpy3_user_arg* a){

  REAL* Out = a->Out;
  REAL* scalep = a->scalep;
  REAL* InScale = a->InScale;
  REAL* Add = a->Add;

  int n_3vec = hi - lo;
  int index = lo * 24;
  Out = &Out[index];
  InScale = &InScale[index];
  Add = &Add[index];
  vcaxpy3(Out, scalep, InScale, Add, n_3vec);
}

// structure for vcaxpy3 of NOT having order (with y pointer only)
struct unordered_vcaxpy3_y_user_arg{
  unordered_vcaxpy3_y_user_arg(
			       const OLattice< CTVec > &x_,
			       OLattice< CTVec >& d_,
			       REAL* scalep_,
			       const int* tab_) : x(x_), d(d_), scalep(scalep_), tab(tab_) {}

  const OLattice< CTVec > &x;
  OLattice< CTVec >& d;
  REAL* scalep;
  const int* tab;
};

// user func for vcaxpy3 of NOT having order (with y pointer only)
inline
void unordered_vcaxpy3_y_evaluate_function (int lo, int hi, int myId, unordered_vcaxpy3_y_user_arg* a){
  const OLattice< CTVec > &x = a->x;
  OLattice< CTVec >& d = a->d;
  REAL* scalep = a->scalep;
  const int* tab = a->tab;

  for(int j=lo; j < hi; j++) { 
     int i=tab[j];
   
     REAL* xptr = (REAL *)&(x.elem(i).elem(0).elem(0).real());
     REAL* yptr = &(d.elem(i).elem(0).elem(0).real());

     vcaxpy3(yptr, scalep, xptr, yptr, 1);
    
   }

}

// structure for vcaxpy3 of NOT having order (with z pointer)
struct unordered_vcaxpy3_z_user_arg{
  unordered_vcaxpy3_z_user_arg(
			       const OLattice< CTVec > &x_,
			       const OLattice< CTVec >& y_,
			       OLattice< CTVec >& d_,
			       REAL* scalep_,
			       const int* tab_) : x(x_), y(y_), d(d_), scalep(scalep_),  tab(tab_) {}

  const OLattice< CTVec > &x;
  const OLattice< CTVec >& y;
  OLattice< CTVec >& d;
  REAL* scalep;
  const int* tab;
};

// user func for vcaxpy3 of NOT having order (with z pointer)
inline
void unordered_vcaxpy3_z_evaluate_function (int lo, int hi, int myId, unordered_vcaxpy3_z_user_arg* a){
  const OLattice< CTVec > &x = a->x;
  const OLattice< CTVec >& y = a->y;
  OLattice< CTVec >& d = a->d;
  REAL* scalep = a->scalep;
  const int* tab = a->tab;

  for(int j=lo; j < hi; j++) { 
     int i=tab[j];
   
     REAL *xptr = (REAL *) &(x.elem(i).elem(0).elem(0).real());
     REAL *yptr = (REAL *) &(y.elem(i).elem(0).elem(0).real());
     REAL *zptr =          &(d.elem(i).elem(0).elem(0).real());

     // Get the no of 3vecs. s.start() and s.end() are inclusive so add +1
     vcaxpy3(zptr, scalep, xptr, yptr, 1);
   }

}





//
// vcaxmy3
//

// structure for vcaxmy3 of having order
struct ordered_vcaxmy3_user_arg{
  REAL* Out;
  REAL* scalep;
  REAL* InScale;
  REAL* Sub;
};

// user func for vcaxmy3 of having order
inline
void ordered_vcaxmy3_evaluate_function (int lo, int hi, int myId, ordered_vcaxmy3_user_arg* a){

  REAL* Out = a->Out;
  REAL* scalep = a->scalep;
  REAL* InScale = a->InScale;
  REAL* Sub = a->Sub;

  int n_3vec = hi - lo;
  int index = lo * 24;
  Out = &Out[index];
  InScale = &InScale[index];
  Sub = &Sub[index];
  vcaxmy3(Out, scalep, InScale, Sub, n_3vec);
}

// structure for vcaxmy3 of NOT having order
struct unordered_vcaxmy3_user_arg{
  unordered_vcaxmy3_user_arg(
			     const OLattice< CTVec > &x_,
			     const OLattice< CTVec >& y_,
			     OLattice< CTVec >& d_,
			     REAL* scalep_,
			     const int* tab_) : x(x_),y(y_),d(d_),scalep(scalep_),  tab(tab_) {}

  const OLattice< CTVec > &x;
  const OLattice< CTVec >& y;
  OLattice< CTVec >& d;
  REAL* scalep;
  const int* tab;
};

// user func for vcaxmy3 of NOT having order 
inline
void unordered_vcaxmy3_evaluate_function (int lo, int hi, int myId, unordered_vcaxmy3_user_arg* a){
  const OLattice< CTVec > &x = a->x;
  const OLattice< CTVec >& y = a->y;
  OLattice< CTVec >& d = a->d;
  REAL* scalep = a->scalep;
  const int* tab = a->tab;

  for(int j=lo; j < hi; j++) { 
     int i=tab[j];
   
     REAL *xptr = (REAL *) &(x.elem(i).elem(0).elem(0).real());
     REAL *yptr = (REAL *) &(y.elem(i).elem(0).elem(0).real());
     REAL *zptr =          &(d.elem(i).elem(0).elem(0).real());

     // Get the no of 3vecs. s.start() and s.end() are inclusive so add +1
     vcaxmy3(zptr, scalep, xptr, yptr, 1);
   }

}



//
// vcaxpby3
//

// structure for vcaxpby3 of having order
struct ordered_vcaxpby3_user_arg{
  REAL* Out;
  REAL* ap;
  REAL* xp;
  REAL* bp;
  REAL* yp;
};

// user func for vcaxpby3 of having order
inline
void ordered_vcaxpby3_evaluate_function (int lo, int hi, int myId, ordered_vcaxpby3_user_arg* a){

  REAL* Out = a->Out;
  REAL* ap = a->ap;
  REAL* xp = a->xp;
  REAL* bp = a->bp;
  REAL* yp = a->yp;

  int n_3vec = hi - lo;
  int index = lo * 24;
  Out = &Out[index];
  xp = &xp[index];
  yp = &yp[index];
  vcaxpby3(Out, ap, xp, bp, yp, n_3vec);
}

// structure for vcaxpby3 of NOT having order
struct unordered_vcaxpby3_user_arg{
  unordered_vcaxpby3_user_arg(
			      const OLattice< CTVec > &x_,
			      const OLattice< CTVec >& y_,
			      OLattice< CTVec >& d_,
			      REAL* aptr_,
			      REAL* bptr_,
			      const int* tab_) : x(x_), y(y_),d(d_),aptr(aptr_), bptr(bptr_), tab(tab_) {}

  const OLattice< CTVec > &x;
  const OLattice< CTVec >& y;
  OLattice< CTVec >& d;
  REAL* aptr;
  REAL* bptr;
  const int* tab;
};

// user func for vcaxpby3 of NOT having order 
inline
void unordered_vcaxpby3_evaluate_function (int lo, int hi, int myId, unordered_vcaxpby3_user_arg* a){
  const OLattice< CTVec > &x = a->x;
  const OLattice< CTVec >& y = a->y;
  OLattice< CTVec >& d = a->d;
  REAL* aptr = a->aptr;
  REAL* bptr = a->bptr;
  const int* tab = a->tab;

  for(int j=lo; j < hi; j++) { 
     int i=tab[j];
   
     REAL *xptr = (REAL *) &(x.elem(i).elem(0).elem(0).real());
     REAL *yptr = (REAL *) &(y.elem(i).elem(0).elem(0).real());
     REAL *zptr =          &(d.elem(i).elem(0).elem(0).real());

     // Get the no of 3vecs. s.start() and s.end() are inclusive so add +1
     vcaxpby3(zptr, aptr, xptr, bptr, yptr, 1);  
   }

}



//
// vcaxmby3
//

// structure for vcaxmby3 of having order
struct ordered_vcaxmby3_user_arg{
  REAL* Out;
  REAL* ap;
  REAL* xp;
  REAL* bp;
  REAL* yp;
};

// user func for vcaxmby3 of having order
inline
void ordered_vcaxmby3_evaluate_function (int lo, int hi, int myId, ordered_vcaxmby3_user_arg* a){

  REAL* Out = a->Out;
  REAL* ap = a->ap;
  REAL* xp = a->xp;
  REAL* bp = a->bp;
  REAL* yp = a->yp;

  int n_3vec = hi - lo;
  int index = lo * 24;
  Out = &Out[index];
  xp = &xp[index];
  yp = &yp[index];
  vcaxmby3(Out, ap, xp, bp, yp, n_3vec);
}

// structure for vcaxmby3 of NOT having order
struct unordered_vcaxmby3_user_arg{
  unordered_vcaxmby3_user_arg(
			      const OLattice< CTVec > &x_,
			      const OLattice< CTVec >& y_,
			      OLattice< CTVec >& d_,
			      REAL* aptr_,
			      REAL* bptr_,
			      const int* tab_): x(x_),y(y_),d(d_),aptr(aptr_),bptr(bptr_),tab(tab_) {}

  const OLattice< CTVec > &x;
  const OLattice< CTVec >& y;
  OLattice< CTVec >& d;
  REAL* aptr;
  REAL* bptr;
  const int* tab;
};

// user func for vcaxmby3 of NOT having order 
inline
void unordered_vcaxmby3_evaluate_function (int lo, int hi, int myId, unordered_vcaxmby3_user_arg* a){
  const OLattice< CTVec > &x = a->x;
  const OLattice< CTVec >& y = a->y;
  OLattice< CTVec >& d = a->d;
  REAL* aptr = a->aptr;
  REAL* bptr = a->bptr;
  const int* tab = a->tab;

  for(int j=lo; j < hi; j++) { 
     int i=tab[j];
   
     REAL *xptr = (REAL *) &(x.elem(i).elem(0).elem(0).real());
     REAL *yptr = (REAL *) &(y.elem(i).elem(0).elem(0).real());
     REAL *zptr =          &(d.elem(i).elem(0).elem(0).real());

     // Get the no of 3vecs. s.start() and s.end() are inclusive so add +1
     vcaxmby3(zptr, aptr, xptr, bptr, yptr, 1);  
   }

}



#endif
