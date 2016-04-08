#ifndef QDP_SCALARSITE_GENERIC_BLAS_WRAPPER_H
#define QDP_SCALARSITE_GENERIC_BLAS_WRAPPER_H

////////////////////////////////
// Threading evaluates wrappers
//
// by Xu Guo, EPCC, 12 August, 2008
////////////////////////////////

//
// for vaxpy3
//

// structure for vaxpy3 of having order
struct ordered_vaxpy3_user_arg{
  REAL* Out;
  REAL* scalep;
  REAL* InScale;
  REAL* Add;
};

// user func for vaxpy3 of having order
inline
void ordered_vaxpy3_evaluate_function (int lo, int hi, int myId, ordered_vaxpy3_user_arg* a){

  REAL* Out = a->Out;
  REAL* scalep = a->scalep;
  REAL* InScale = a->InScale;
  REAL* Add = a->Add;

  int n_4vec = hi - lo;

  int index = lo * 24;
  InScale = &InScale[index];
  Add = &Add[index];
  Out = &Out[index];

  vaxpy3(Out, scalep, InScale, Add, n_4vec);
}

// structure for vaxpy3 (with yptr only) of NOT having order
struct unordered_vaxpy3_y_user_arg{
  unordered_vaxpy3_y_user_arg(  const OLattice< TVec >& x_,
				OLattice< TVec >& d_,
				REAL* scalep_,
				const int* tab_,
				int xy_order_) : x(x_), d(d_), scalep(scalep_), tab(tab_), xy_order(xy_order_) {}

  const OLattice< TVec >& x;
  OLattice< TVec >& d;
  REAL* scalep;
  const int* tab;
  int xy_order;
};

// user func for vaxpy3 (with yptr only) of NOT having order
inline
void unordered_vaxpy3_y_evaluate_function (int lo, int hi, int myId, unordered_vaxpy3_y_user_arg* a){

  const OLattice< TVec >& x = a->x;
  OLattice< TVec >& d = a->d;
  REAL* scalep = a->scalep;
  const int* tab = a->tab;
  int xy_order = a->xy_order;

  if (xy_order){
     for(int j=lo; j < hi; j++) { 
       int i=tab[j];
       REAL* xptr = (REAL *)&(x.elem(i).elem(0).elem(0).real());
       REAL* yptr = &(d.elem(i).elem(0).elem(0).real());
       vaxpy3(yptr, scalep, xptr, yptr, 1);
     }
  }
  else {
    for(int j=lo; j < hi; j++) { 
     int i=tab[j];
     REAL* xptr = (REAL *)&(x.elem(i).elem(0).elem(0).real());
     REAL* yptr = &(d.elem(i).elem(0).elem(0).real());
     vaxpy3(yptr, scalep, yptr, xptr, 1);
   }
  }
 

}

// structure for vaxpy3 (with zptr) of NOT having order
struct unordered_vaxpy3_z_user_arg{
  unordered_vaxpy3_z_user_arg(  const OLattice< TVec >& x_,
				const OLattice< TVec >& y_,
				OLattice< TVec >& d_,
				REAL* scalep_,
				const int* tab_) : x(x_), y(y_), d(d_), scalep(scalep_), tab(tab_) {}

  const OLattice< TVec >& x;
  const OLattice< TVec >& y;
  OLattice< TVec >& d;
  REAL* scalep;
  const int* tab;
};

// user func for vaxpy3 (with zptr only) of NOT having order
inline
void unordered_vaxpy3_z_evaluate_function (int lo, int hi, int myId, unordered_vaxpy3_z_user_arg* a){

  const OLattice< TVec >& x = a->x;
  const OLattice< TVec >& y = a->y;
  OLattice< TVec >& d = a->d;
  REAL* scalep = a->scalep;
  const int* tab = a->tab;

  for(int j=lo; j < hi; j++) { 
     int i=tab[j];
     REAL* xptr = (REAL *)&(x.elem(i).elem(0).elem(0).real());
     REAL* yptr = (REAL *)&(y.elem(i).elem(0).elem(0).real());
     REAL* zptr =  &(d.elem(i).elem(0).elem(0).real());
     vaxpy3(zptr, scalep, xptr, yptr, 1);
   }

}


//
// for vaxmy3
//

// structure for vaxmy3 of having order
struct ordered_vaxmy3_user_arg{
  REAL* Out;
  REAL* scalep;
  REAL* InScale;
  REAL* Sub;
};

// user func for vaxmy3 of having order
inline
void ordered_vaxmy3_evaluate_function (int lo, int hi, int myId, ordered_vaxmy3_user_arg* a){

  REAL* Out = a->Out;
  REAL* scalep = a->scalep;
  REAL* InScale = a->InScale;
  REAL* Sub = a->Sub;

  int n_3vec = hi - lo;
  int index = lo * 24;
  InScale = &InScale[index];
  Sub = &Sub[index];
  Out = &Out[index];

  vaxmy3(Out, scalep, InScale, Sub, n_3vec);

}

// structure for vaxmy3 (with yptr only) of NOT having order
struct unordered_vaxmy3_y_user_arg {
  unordered_vaxmy3_y_user_arg(  const OLattice< TVec >& x_,
				OLattice< TVec >& d_,
				REAL* scalep_,
				const int* tab_) : x(x_), d(d_), scalep(scalep_), tab(tab_) {}

  const OLattice< TVec >& x;
  OLattice< TVec >& d;
  REAL* scalep;
  const int* tab;
};

// user func for vampy3 (with yptr only) of NOT having order
inline
void unordered_vaxmy3_y_evaluate_function (int lo, int hi, int myId, unordered_vaxmy3_y_user_arg* a){

  const OLattice< TVec >& x = a->x;
  OLattice< TVec >& d = a->d;
  REAL* scalep = a->scalep;
  const int* tab = a->tab;

  for(int j=lo; j < hi; j++) { 
    int i=tab[j];
    REAL* xptr = (REAL *)&(x.elem(i).elem(0).elem(0).real());
    REAL* yptr = &(d.elem(i).elem(0).elem(0).real());
    vaxmy3(yptr, scalep, yptr, xptr, 1);
  }
}

// structure for vaxmy3 (with zptr) of NOT having order
struct unordered_vaxmy3_z_user_arg{
  unordered_vaxmy3_z_user_arg(
			      const OLattice< TVec >& x_,
			      const OLattice< TVec >& y_,
			      OLattice< TVec >& d_,
			      REAL* scalep_,
			      const int* tab_) : x(x_), y(y_), d(d_), scalep(scalep_), tab(tab_) {}

  const OLattice< TVec >& x;
  const OLattice< TVec >& y;
  OLattice< TVec >& d;
  REAL* scalep;
  const int* tab;
};

// user func for vaxmy3 (with zptr only) of NOT having order
inline
void unordered_vaxmy3_z_evaluate_function (int lo, int hi, int myId, unordered_vaxmy3_z_user_arg* a){

  const OLattice< TVec >& x = a->x;
  const OLattice< TVec >& y = a->y;
  OLattice< TVec >& d = a->d;
  REAL* scalep = a->scalep;
  const int* tab = a->tab;

  for(int j=lo; j < hi; j++) { 
     int i=tab[j];
     REAL* xptr = (REAL *)&(x.elem(i).elem(0).elem(0).real());
     REAL* yptr = (REAL *)&(y.elem(i).elem(0).elem(0).real());
     REAL* zptr =  &(d.elem(i).elem(0).elem(0).real());
     vaxmy3(zptr, scalep, xptr, yptr, 1);
   }

}


//
// vscal
//

// structure for vscal of having order
struct ordered_vscal_user_arg{
  REAL* Out;
  REAL* scalep;
  REAL* In;
};

// user func for vscal of having order
inline
void ordered_vscal_evaluate_function (int lo, int hi, int myId, ordered_vscal_user_arg* a){

  REAL* Out = a->Out;
  REAL* scalep = a->scalep;
  REAL* In = a->In;

  int n_3vec = hi - lo;
  int index = lo * 24;
  Out = &Out[index];
  In = &In[index];
  vscal(Out, scalep, In, n_3vec);
}

// structure for vscal of NOT having order
struct unordered_vscal_user_arg {
  unordered_vscal_user_arg(
  
			   const OLattice< TVec >& x_,
			   OLattice< TVec >& d_,
			   REAL* scalep_,
			   const int* tab_
			   ) : x(x_), d(d_), scalep(scalep_), tab(tab_) {}
  const OLattice< TVec >& x;
  OLattice< TVec >& d;
  REAL* scalep;
  const int* tab;
};

// user func for vscal of NOT having order
inline
void unordered_vscal_evaluate_function (int lo, int hi, int myId, unordered_vscal_user_arg* a){

  const OLattice< TVec >& x = a->x;
  OLattice< TVec >& d = a->d;
  REAL* scalep = a->scalep;
  const int* tab = a->tab;

  for(int j=lo; j < hi; j++) { 
     int i=tab[j];
     REAL *xptr = (REAL *) &(x.elem(i).elem(0).elem(0).real());
     REAL *zptr =  &(d.elem(i).elem(0).elem(0).real()); 
     vscal(zptr, scalep, xptr, 1);
   }

}



//
// vaxpby3
//

// structure for vaxpby3 of having order
struct ordered_vaxpby3_user_arg{
  REAL* Out;
  REAL* ap;
  REAL* xp;
  REAL* bp;
  REAL* yp;
};

// user func for vaxpby3 of having order
inline
void ordered_vaxpby3_evaluate_function (int lo, int hi, int myId, ordered_vaxpby3_user_arg* a){

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
  vaxpby3(Out, ap, xp, bp, yp, n_3vec);
}

// structure for vaxpby3 of NOT having order
struct unordered_vaxpby3_user_arg{
unordered_vaxpby3_user_arg(  const OLattice< TVec >& x_,
			     const OLattice< TVec >& y_,
			     OLattice< TVec >& d_,
			     REAL* aptr_,
			     REAL* bptr_,
			     const int* tab_): x(x_), y(y_), d(d_), aptr(aptr_), bptr(bptr_), tab(tab_) {}
  const OLattice< TVec >& x;
  const OLattice< TVec >& y;
  OLattice< TVec >& d;
  REAL* aptr;
  REAL* bptr;
  const int* tab;
};

// user func for vaxpby3 of NOT having order
inline
void unordered_vaxpby3_evaluate_function (int lo, int hi, int myId, unordered_vaxpby3_user_arg* a){

  const OLattice< TVec >& x = a->x;
  const OLattice< TVec >& y = a->y;
  OLattice< TVec >& d = a->d;
  REAL* aptr = a->aptr;
  REAL* bptr = a->bptr;
  const int* tab = a->tab;

  for(int j=lo; j < hi; j++) { 
     int i=tab[j];
     REAL *xptr = (REAL *) &(x.elem(i).elem(0).elem(0).real());
     REAL *yptr = (REAL *) &(y.elem(i).elem(0).elem(0).real());
     REAL* zptr =  &(d.elem(i).elem(0).elem(0).real());
      
     // Get the no of 3vecs. s.start() and s.end() are inclusive so add +1
     vaxpby3(zptr, aptr, xptr, bptr, yptr, 1);
  
   }

}



//
// vaxmby3
//

// structure for vaxmby3 of having order
struct ordered_vaxmby3_user_arg{
  REAL* Out;
  REAL* ap;
  REAL* xp;
  REAL* bp;
  REAL* yp;
};

// user func for vaxpby3 of having order
inline
void ordered_vaxmby3_evaluate_function (int lo, int hi, int myId, ordered_vaxmby3_user_arg* a){

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
  vaxmby3(Out, ap, xp, bp, yp, n_3vec);
}

// structure for vaxmby3 of NOT having order
struct unordered_vaxmby3_user_arg{
  unordered_vaxmby3_user_arg(  const OLattice< TVec >& x_,
			       const OLattice< TVec >& y_,
			       OLattice< TVec >& d_,
			       REAL* aptr_,
			       REAL* bptr_,
			       const int* tab_): x(x_), y(y_), d(d_), aptr(aptr_), bptr(bptr_), tab(tab_) {}

  const OLattice< TVec >& x;
  const OLattice< TVec >& y;
  OLattice< TVec >& d;
  REAL* aptr;
  REAL* bptr;
  const int* tab;
};

// user func for vaxmby3 of NOT having order
inline
void unordered_vaxmby3_evaluate_function (int lo, int hi, int myId, unordered_vaxmby3_user_arg* a){

  const OLattice< TVec >& x = a->x;
  const OLattice< TVec >& y = a->y;
  OLattice< TVec >& d = a->d;
  REAL* aptr = a->aptr;
  REAL* bptr = a->bptr;
  const int* tab = a->tab;

  for(int j=lo; j < hi; j++) { 
     int i=tab[j];
     REAL *xptr = (REAL *) &(x.elem(i).elem(0).elem(0).real());
     REAL *yptr = (REAL *) &(y.elem(i).elem(0).elem(0).real());
     REAL* zptr =  &(d.elem(i).elem(0).elem(0).real());
      
     // Get the no of 3vecs. s.start() and s.end() are inclusive so add +1
     vaxmby3(zptr, aptr, xptr, bptr, yptr,1);
  
   }

}







#endif
