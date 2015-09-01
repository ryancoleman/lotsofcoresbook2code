#ifndef GENERIC_MAT_VEC_H
#define GENERIC_MAT_VEC_H
/* 
 * Generic routines for multiplying SU(3) matrix with Vector
 * and SU(3) matrix dagger with color vector
 * 
 * Author: $Id: generic_mat_vec.h,v 1.3 2007-02-06 15:01:58 bjoo Exp $
*/

/* SU(3)*color vector:  aa * bb -> cc  */
/* aa = QDP++ PColorMatrix
   bb = QDP++ PColorVector
   cc = QDP++ PColorVector */

/* Ideally want to registerize b, have a in cache, and store gather c? */
#define _inline_generic_mult_su3_mat_vec(aa,bb,cc) \
{\
  cc.elem(0).real()  = aa.elem(0,0).real()*bb.elem(0).real() ; \
  cc.elem(0).real() -= aa.elem(0,0).imag()*bb.elem(0).imag() ; \
\
  cc.elem(0).real() += aa.elem(0,1).real()*bb.elem(1).real() ; \
  cc.elem(0).real() -= aa.elem(0,1).imag()*bb.elem(1).imag() ; \
\
  cc.elem(0).real() += aa.elem(0,2).real()*bb.elem(2).real() ; \
  cc.elem(0).real() -= aa.elem(0,2).imag()*bb.elem(2).imag() ; \
\
  cc.elem(0).imag()  = aa.elem(0,0).real()*bb.elem(0).imag() ; \
  cc.elem(0).imag() += aa.elem(0,0).imag()*bb.elem(0).real() ; \
\
  cc.elem(0).imag() += aa.elem(0,1).real()*bb.elem(1).imag() ; \
  cc.elem(0).imag() += aa.elem(0,1).imag()*bb.elem(1).real() ; \
\
  cc.elem(0).imag() += aa.elem(0,2).real()*bb.elem(2).imag() ; \
  cc.elem(0).imag() += aa.elem(0,2).imag()*bb.elem(2).real() ; \
\
\
  cc.elem(1).real()  = aa.elem(1,0).real()*bb.elem(0).real() ; \
  cc.elem(1).real() -= aa.elem(1,0).imag()*bb.elem(0).imag() ; \
\
  cc.elem(1).real() += aa.elem(1,1).real()*bb.elem(1).real() ; \
  cc.elem(1).real() -= aa.elem(1,1).imag()*bb.elem(1).imag() ; \
\
  cc.elem(1).real() += aa.elem(1,2).real()*bb.elem(2).real() ; \
  cc.elem(1).real() -= aa.elem(1,2).imag()*bb.elem(2).imag() ; \
\
  cc.elem(1).imag()  = aa.elem(1,0).real()*bb.elem(0).imag() ; \
  cc.elem(1).imag() += aa.elem(1,0).imag()*bb.elem(0).real() ; \
\
  cc.elem(1).imag() += aa.elem(1,1).real()*bb.elem(1).imag() ; \
  cc.elem(1).imag() += aa.elem(1,1).imag()*bb.elem(1).real() ; \
\
  cc.elem(1).imag() += aa.elem(1,2).real()*bb.elem(2).imag() ; \
  cc.elem(1).imag() += aa.elem(1,2).imag()*bb.elem(2).real() ; \
\
\
  cc.elem(2).real()  = aa.elem(2,0).real()*bb.elem(0).real() ; \
  cc.elem(2).real() -= aa.elem(2,0).imag()*bb.elem(0).imag() ; \
\
  cc.elem(2).real() += aa.elem(2,1).real()*bb.elem(1).real() ; \
  cc.elem(2).real() -= aa.elem(2,1).imag()*bb.elem(1).imag() ; \
\
  cc.elem(2).real() += aa.elem(2,2).real()*bb.elem(2).real() ; \
  cc.elem(2).real() -= aa.elem(2,2).imag()*bb.elem(2).imag() ; \
\
  cc.elem(2).imag()  = aa.elem(2,0).real()*bb.elem(0).imag() ; \
  cc.elem(2).imag() += aa.elem(2,0).imag()*bb.elem(0).real() ; \
\
  cc.elem(2).imag() += aa.elem(2,1).real()*bb.elem(1).imag() ; \
  cc.elem(2).imag() += aa.elem(2,1).imag()*bb.elem(1).real() ; \
\
  cc.elem(2).imag() += aa.elem(2,2).real()*bb.elem(2).imag() ; \
  cc.elem(2).imag() += aa.elem(2,2).imag()*bb.elem(2).real() ; \
}

/* Unrolled adj SU(3)*color vector: aa^\dagger * bb -> cc 
 *  aa = QDP++ PColorMatrix
 *  bb = QDP++ PColorVector
 *  cc = QDP++ PColorVector
 */
/* Ideally want to registerize b, have a in cache, and store gather c? */
#define _inline_generic_mult_adj_su3_mat_vec(aa,bb,cc) \
{\
  cc.elem(0).real()  = aa.elem(0,0).real()*bb.elem(0).real() ; \
  cc.elem(0).real() += aa.elem(0,0).imag()*bb.elem(0).imag() ; \
\
  cc.elem(0).real() += aa.elem(1,0).real()*bb.elem(1).real() ; \
  cc.elem(0).real() += aa.elem(1,0).imag()*bb.elem(1).imag() ; \
\
  cc.elem(0).real() += aa.elem(2,0).real()*bb.elem(2).real() ; \
  cc.elem(0).real() += aa.elem(2,0).imag()*bb.elem(2).imag() ; \
\
  cc.elem(0).imag()  = aa.elem(0,0).real()*bb.elem(0).imag() ; \
  cc.elem(0).imag() -= aa.elem(0,0).imag()*bb.elem(0).real() ; \
\
  cc.elem(0).imag() += aa.elem(1,0).real()*bb.elem(1).imag() ; \
  cc.elem(0).imag() -= aa.elem(1,0).imag()*bb.elem(1).real() ; \
\
  cc.elem(0).imag() += aa.elem(2,0).real()*bb.elem(2).imag() ; \
  cc.elem(0).imag() -= aa.elem(2,0).imag()*bb.elem(2).real() ; \
\
\
  cc.elem(1).real()  = aa.elem(0,1).real()*bb.elem(0).real() ; \
  cc.elem(1).real() += aa.elem(0,1).imag()*bb.elem(0).imag() ; \
\
  cc.elem(1).real() += aa.elem(1,1).real()*bb.elem(1).real() ; \
  cc.elem(1).real() += aa.elem(1,1).imag()*bb.elem(1).imag() ; \
\
  cc.elem(1).real() += aa.elem(2,1).real()*bb.elem(2).real() ; \
  cc.elem(1).real() += aa.elem(2,1).imag()*bb.elem(2).imag() ; \
\
  cc.elem(1).imag()  = aa.elem(0,1).real()*bb.elem(0).imag() ; \
  cc.elem(1).imag() -= aa.elem(0,1).imag()*bb.elem(0).real() ; \
\
  cc.elem(1).imag() += aa.elem(1,1).real()*bb.elem(1).imag() ; \
  cc.elem(1).imag() -= aa.elem(1,1).imag()*bb.elem(1).real() ; \
\
  cc.elem(1).imag() += aa.elem(2,1).real()*bb.elem(2).imag() ; \
  cc.elem(1).imag() -= aa.elem(2,1).imag()*bb.elem(2).real() ; \
\
\
  cc.elem(2).real()  = aa.elem(0,2).real()*bb.elem(0).real() ; \
  cc.elem(2).real() += aa.elem(0,2).imag()*bb.elem(0).imag() ; \
\
  cc.elem(2).real() += aa.elem(1,2).real()*bb.elem(1).real() ; \
  cc.elem(2).real() += aa.elem(1,2).imag()*bb.elem(1).imag() ; \
\
  cc.elem(2).real() += aa.elem(2,2).real()*bb.elem(2).real() ; \
  cc.elem(2).real() += aa.elem(2,2).imag()*bb.elem(2).imag() ; \
\
  cc.elem(2).imag()  = aa.elem(0,2).real()*bb.elem(0).imag() ; \
  cc.elem(2).imag() -= aa.elem(0,2).imag()*bb.elem(0).real() ; \
\
  cc.elem(2).imag() += aa.elem(1,2).real()*bb.elem(1).imag() ; \
  cc.elem(2).imag() -= aa.elem(1,2).imag()*bb.elem(1).real() ; \
\
  cc.elem(2).imag() += aa.elem(2,2).real()*bb.elem(2).imag() ; \
  cc.elem(2).imag() -= aa.elem(2,2).imag()*bb.elem(2).real() ; \
}

#endif
