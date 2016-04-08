#ifndef SSE_SPIN_PROJ_INLINES_H
#define SSE_SPIN_PROJ_INLINES_H

#include "qdp_sse_intrin.h"

/* File: generic_spin_proj_inlines.h
   Purpose: Supply inline functions to do spin projection
   Author: $Id: sse_spin_proj_inlines.h,v 1.6 2009-02-11 20:50:45 bjoo Exp $
*/
namespace QDP {
#include <stdio.h>

/** \brief Spin Project (1/2)(1+\gamma_0)
 *
 * \param src (pointer to 4 vector source)
 * \param dst (pointer to 2 vector dest)
 * \param n_vec (number of vectors to project)

 * It is assumeed that src points to an array of  4 spin
 * 3 colour vectors with complex components. 
 * It is assumed that dst points to an array of 2 spin
 * 3 colour vectors with complex components. In both cases
 * ordering is that spin is slowest, and complex components
 * are fastest. */
inline  
void inlineSpinProjDir0Plus(const REAL32* src, REAL32 *dst, unsigned int n_vec) 
{

#ifdef DEBUG_GENERIC_SPIN_PROJ_INLINES
  QDPIO::cout << "inlineSpinProjDir0Plus" << endl;
#endif


  SSEVec v0, v1, v2, v3, v4, v5, v6, v7;

  /* 1 + \gamma_0 =  1  0  0  i 
                     0  1  i  0
                     0 -i  1  0
                    -i  0  0  1 
 
   *      ( d0r + i d0i )  =  ( {x0r - x3i} + i{x0i + x3r} )
   *      ( d1r + i d1i )     ( {x1r - x2i} + i{x1i + x2r} )
   */
  const REAL32* src_shadow = src;
  REAL32* dst_shadow = dst;


  // Store in the spinor - top half
  v0.vector = _mm_load_ps(src_shadow);
  v1.vector = _mm_load_ps(src_shadow+4);
  v2.vector = _mm_load_ps(src_shadow+8);
  v3.vector = _mm_load_ps(src_shadow+12);
  v4.vector = _mm_load_ps(src_shadow+16);
  v5.vector = _mm_load_ps(src_shadow+20);
  
  v7.floats[0] = -1;
  v7.floats[1] = +1;
  v7.floats[2] = -1;
  v7.floats[3] = +1;

  for(unsigned int site=0; site < n_vec-1; site++) {
    src_shadow += 24;

    _mm_prefetch((const char *)src_shadow, _MM_HINT_T0);

    v6.vector = v4.vector;         // V6 is dest so we can move its 
    // I want to shuffle so that:
    //  v6[0] <- v6[3]=v4[3]  [1:0]=3 = x11 <- tmp[3][col=0][im]
    //  v6[1] <- v6[2]=v4[2]  [3:2]=2 = x10 <- tmp[3][col=0][re]
    //  v6[2] <-       v5[1]  [5:4]=1 = x01 <- tmp[3][col=1][im]
    //  v6[3] <-       v5[0]  [7:6]=0 = x00 <- tmp[3][col=1][re]
    //  So Immediate for shufps is: 0x0001 1011 = x1B = 27
    v6.vector = _mm_shuffle_ps(v6.vector, v5.vector, 0x1B);
    
    // Does a v0 + v6*v7  
    //  v0[0] <- v0[0] + v6[0]*v7[0] =  tmp[0][col=0][re] - tmp[3][col=0][im]
    //  v0[1] <- v0[1] + v6[1]*v7[1] =  tmp[0][col=0][im] + tmp[3][col=0][re]
    //  v0[2] <- v0[2] + v6[2]*v7[2] =  tmp[0][col=1][re] - tmp[3][col=1][im]
    //  v0[3] <- v0[3] + v6[3]*v7[3] =  tmp[0][col=1][im] + tmp[3][col=1][re]
    v6.vector = _mm_mul_ps(v6.vector, v7.vector);
    v0.vector = _mm_add_ps(v0.vector, v6.vector);
    _mm_store_ps(dst_shadow, v0.vector);

    v6.vector = v5.vector;
    // Now setup v7 so that 
    //  v6[0] <- v5[3]  [1:0]=3 = x11 <- tmp[3][col=2][im]
    //  v6[1] <- v5[2]  [3:2]=2 = x10 <- tmp[3][col=2][re]
    //  v6[2] <- v3[1]  [5:4]=1 = x01 <- tmp[2][col=0][im]
    //  v6[3] <- v3[0]  [7:6]=0 = x00 <- tmp[2][col=0][re]
    //  So Immediate for shufps is: 0x0001 1011 = x1B = 27
    v6.vector = _mm_shuffle_ps(v6.vector, v3.vector, 0x1B);

    // Does a v1 + v6*v7  
    //  v1[0] <- v1[0] + v6[0]*v7[0] =  tmp[0][col=2][re] - tmp[3][col=2][im]
    //  v1[1] <- v1[1] + v6[1]*v7[1] =  tmp[0][col=2][im] + tmp[3][col=2][re]
    //  v1[2] <- v1[2] + v6[2]*v7[2] =  tmp[1][col=0][re] - tmp[2][col=0][im]
    //  v1[3] <- v1[3] + v6[3]*v7[3] =  tmp[1][col=0][im] + tmp[2][col=0][re]
    v6.vector  = _mm_mul_ps(v6.vector, v7.vector);
    v1.vector  = _mm_add_ps(v1.vector, v6.vector);
    _mm_store_ps(dst_shadow+4, v1.vector);

    v6.vector = v3.vector;
    //  v6[0] <- v6[3]=v3[3]  [1:0]=3 = x11 <- tmp[2][col=1][im]
    //  v6[1] <- v6[3]=v3[2]  [3:2]=2 = x10 <- tmp[2][col=1][re]
    //  v6[2] <-       v4[1]  [5:4]=1 = x01 <- tmp[2][col=2][im]
    //  v6[3] <-       v4[0]  [7:6]=0 = x00 <- tmp[2][col=2][re]
    //  So Immediate for shufps is: 0x0001 1011 = x1B = 27
    v6.vector = _mm_shuffle_ps(v6.vector, v4.vector, 0x1B);

    // Does a v2 + v6*v7  
    //  v2[0] <- v2[0] + v6[0]*v7[0] =  tmp[1][col=1][re] - tmp[2][col=1][im]
    //  v2[1] <- v2[1] + v6[1]*v7[1] =  tmp[1][col=1][im] + tmp[2][col=1][re]
    //  v2[2] <- v2[2] + v6[2]*v7[2] =  tmp[1][col=2][re] - tmp[2][col=2][im]
    //  v2[3] <- v2[3] + v6[3]*v7[3] =  tmp[1][col=2][im] + tmp[2][col=2][re]
    v6.vector = _mm_mul_ps(v6.vector, v7.vector);
    v2.vector = _mm_add_ps(v2.vector, v6.vector);
    _mm_store_ps(dst_shadow+8, v2.vector);

    // Push the next one
    dst_shadow+=12;

    // Store in the next spinor -- should be prefetched if all has gone well
    v0.vector = _mm_load_ps(src_shadow);
    v1.vector = _mm_load_ps(src_shadow+4);
    v2.vector = _mm_load_ps(src_shadow+8);
    v3.vector = _mm_load_ps(src_shadow+12);
    v4.vector = _mm_load_ps(src_shadow+16);
    v5.vector = _mm_load_ps(src_shadow+20);

  }
  
  // Last bit
  
  v6.vector = v4.vector;         // V6 is dest so we can move its 
  // I want to shuffle so that:
  //  v6[0] <- v6[3]=v4[3]  [1:0]=3 = x11 <- tmp[3][col=0][im]
  //  v6[1] <- v6[2]=v4[2]  [3:2]=2 = x10 <- tmp[3][col=0][re]
  //  v6[2] <-       v5[1]  [5:4]=1 = x01 <- tmp[3][col=1][im]
  //  v6[3] <-       v5[0]  [7:6]=0 = x00 <- tmp[3][col=1][re]
  //  So Immediate for shufps is: 0x0001 1011 = x1B = 27
  v6.vector = _mm_shuffle_ps(v6.vector, v5.vector, 0x1B);
  
  // Does a v0 + v6*v7  
  //  v0[0] <- v0[0] + v6[0]*v7[0] =  tmp[0][col=0][re] - tmp[3][col=0][im]
  //  v0[1] <- v0[1] + v6[1]*v7[1] =  tmp[0][col=0][im] + tmp[3][col=0][re]
  //  v0[2] <- v0[2] + v6[2]*v7[2] =  tmp[0][col=1][re] - tmp[3][col=1][im]
  //  v0[3] <- v0[3] + v6[3]*v7[3] =  tmp[0][col=1][im] + tmp[3][col=1][re]
  v6.vector = _mm_mul_ps(v6.vector, v7.vector);
  v0.vector = _mm_add_ps(v0.vector, v6.vector);
  _mm_store_ps(dst_shadow, v0.vector);
  
  v6.vector = v5.vector;
  // Now setup v7 so that 
  //  v6[0] <- v5[3]  [1:0]=3 = x11 <- tmp[3][col=2][im]
  //  v6[1] <- v5[2]  [3:2]=2 = x10 <- tmp[3][col=2][re]
  //  v6[2] <- v3[1]  [5:4]=1 = x01 <- tmp[2][col=0][im]
  //  v6[3] <- v3[0]  [7:6]=0 = x00 <- tmp[2][col=0][re]
  //  So Immediate for shufps is: 0x0001 1011 = x1B = 27
  v6.vector = _mm_shuffle_ps(v6.vector, v3.vector, 0x1B);
  
  // Does a v1 + v6*v7  
  //  v1[0] <- v1[0] + v6[0]*v7[0] =  tmp[0][col=2][re] - tmp[3][col=2][im]
  //  v1[1] <- v1[1] + v6[1]*v7[1] =  tmp[0][col=2][im] + tmp[3][col=2][re]
  //  v1[2] <- v1[2] + v6[2]*v7[2] =  tmp[1][col=0][re] - tmp[2][col=0][im]
  //  v1[3] <- v1[3] + v6[3]*v7[3] =  tmp[1][col=0][im] + tmp[2][col=0][re]
  v6.vector  = _mm_mul_ps(v6.vector, v7.vector);
  v1.vector  = _mm_add_ps(v1.vector, v6.vector);
  _mm_store_ps(dst_shadow+4, v1.vector);
  
  v6.vector = v3.vector;
  //  v6[0] <- v6[3]=v3[3]  [1:0]=3 = x11 <- tmp[2][col=1][im]
  //  v6[1] <- v6[3]=v3[2]  [3:2]=2 = x10 <- tmp[2][col=1][re]
  //  v6[2] <-       v4[1]  [5:4]=1 = x01 <- tmp[2][col=2][im]
  //  v6[3] <-       v4[0]  [7:6]=0 = x00 <- tmp[2][col=2][re]
  //  So Immediate for shufps is: 0x0001 1011 = x1B = 27
  v6.vector = _mm_shuffle_ps(v6.vector, v4.vector, 0x1B);
  
  // Does a v2 + v6*v7  
  //  v2[0] <- v2[0] + v6[0]*v7[0] =  tmp[1][col=1][re] - tmp[2][col=1][im]
  //  v2[1] <- v2[1] + v6[1]*v7[1] =  tmp[1][col=1][im] + tmp[2][col=1][re]
  //  v2[2] <- v2[2] + v6[2]*v7[2] =  tmp[1][col=2][re] - tmp[2][col=2][im]
  //  v2[3] <- v2[3] + v6[3]*v7[3] =  tmp[1][col=2][im] + tmp[2][col=2][re]
  v6.vector = _mm_mul_ps(v6.vector, v7.vector);
  v2.vector = _mm_add_ps(v2.vector, v6.vector);
  _mm_store_ps(dst_shadow+8, v2.vector);
  
}

/** \brief Spin Project (1/2)(1-\gamma_0)
 *
 * \param src (pointer to 4 vector source)
 * \param dst (pointer to 2 vector dest)
 *
 * It is assumeed that src points to an array of  4 spin
 * 3 colour vectors with complex components. 
 * It is assumed that dst points to an array of 2 spin
 * 3 colour vectors with complex components. In both cases
 * ordering is that spin is slowest, and complex components
 * are fastest. */
inline 
void inlineSpinProjDir0Minus(const REAL32* src, REAL32 *dst, unsigned int n_vec) 
{

#ifdef DEBUG_GENERIC_SPIN_PROJ_INLINES
  QDPIO::cout << "inlineSpinProjDir0Minus" << endl;
#endif


  /*                              ( 1  0  0 -i)  ( a0 )    ( a0 - i a3 )
   *  B  :=  ( 1 - Gamma  ) A  =  ( 0  1 -i  0)  ( a1 )  = ( a1 - i a2 )
   *                    0         ( 0  i  1  0)  ( a2 )    ( a2 + i a1 )
   *                              ( i  0  0  1)  ( a3 )    ( a3 + i a0 )

   * Therefore the top components are

   *      ( b0r + i b0i )  =  ( {a0r + a3i} + i{a0i - a3r} )
   *      ( b1r + i b1i )     ( {a1r + a2i} + i{a1i - a2r} )
   */


  const REAL32* src_shadow = src;
  REAL32* dst_shadow = dst;
  SSEVec v0, v1, v2, v3, v4, v5, v6, v7;

  v7.floats[0] = -1;
  v7.floats[1] = +1;
  v7.floats[2] = -1;
  v7.floats[3] = +1;
  
  // Store in the spinor 
  v0.vector = _mm_load_ps(src_shadow);
  v1.vector = _mm_load_ps(src_shadow+4);
  v2.vector = _mm_load_ps(src_shadow+8);
  v3.vector = _mm_load_ps(src_shadow+12);
  v4.vector = _mm_load_ps(src_shadow+16);
  v5.vector = _mm_load_ps(src_shadow+20);

  for(unsigned int site=0; site < n_vec-1; site++) {
    src_shadow += 24;

    _mm_prefetch((const char *)src_shadow, _MM_HINT_T0);
    
    v6.vector = v4.vector;         // V6 is dest so we can move its 
    // I want to shuffle so that:
    //  v6[0] <- v6[3]=v4[3]  [1:0]=3 = x11 <- tmp[3][col=0][im]
    //  v6[1] <- v6[2]=v4[2]  [3:2]=2 = x10 <- tmp[3][col=0][re]
    //  v6[2] <-       v5[1]  [5:4]=1 = x01 <- tmp[3][col=1][im]
    //  v6[3] <-       v5[0]  [7:6]=0 = x00 <- tmp[3][col=1][re]
    //  So Immediate for shufps is: 0x0001 1011 = x1B = 27
    v6.vector = _mm_shuffle_ps(v6.vector, v5.vector, 0x1B);

    // Does a v0 - v6*v7
    //  v0[0] <- v0[0] - v6[0]*v7[0] =  tmp[0][col=0][re] + tmp[3][col=0][im]
    //  v0[1] <- v0[1] - v6[1]*v7[1] =  tmp[0][col=0][im] - tmp[3][col=0][re]
    //  v0[2] <- v0[2] - v6[2]*v7[2] =  tmp[0][col=1][re] + tmp[3][col=1][im]
    //  v0[3] <- v0[3] - v6[3]*v7[3] =  tmp[0][col=1][im] - tmp[3][col=1][re]    
    v6.vector = _mm_mul_ps(v6.vector, v7.vector);
    v0.vector = _mm_sub_ps(v0.vector, v6.vector);
    _mm_store_ps(dst_shadow, v0.vector);

    v6.vector = v5.vector;
    // I want to shuffle so that:
    //  v6[0] <- v6[3]=v5[3]  [1:0]=3 = x11 <- tmp[3][col=0][im]
    //  v6[1] <- v6[2]=v5[2]  [3:2]=2 = x10 <- tmp[3][col=0][re]
    //  v6[2] <-       v3[1]  [5:4]=1 = x01 <- tmp[3][col=1][im]
    //  v6[3] <-       v3[0]  [7:6]=0 = x00 <- tmp[3][col=1][re]
    //  So Immediate for shufps is: 0x0001 1011 = x1B = 27
    v6.vector = _mm_shuffle_ps(v6.vector, v3.vector, 0x1B);

    // Does a v1 - v6*v7
    //  v1[0] <- v1[0] - v6[0]*v7[0] =  tmp[0][col=2][re] + tmp[3][col=2][im]
    //  v1[1] <- v1[1] - v6[1]*v7[1] =  tmp[0][col=2][im] - tmp[3][col=2][re]
    //  v1[2] <- v1[2] - v6[2]*v7[2] =  tmp[1][col=0][re] + tmp[2][col=0][im]
    //  v1[3] <- v1[3] - v6[3]*v7[3] =  tmp[1][col=0][im] - tmp[2][col=0][re]
    v6.vector  = _mm_mul_ps(v6.vector, v7.vector);
    v1.vector  = _mm_sub_ps(v1.vector, v6.vector);
    _mm_store_ps(dst_shadow+4, v1.vector);

    v6.vector = v3.vector;
    //  v6[0] <- v6[3]=v3[3]  [1:0]=3 = x11 <- tmp[2][col=1][im]
    //  v6[1] <- v6[3]=v3[2]  [3:2]=2 = x10 <- tmp[2][col=1][re]
    //  v6[2] <-       v4[1]  [5:4]=1 = x01 <- tmp[2][col=2][im]
    //  v6[3] <-       v4[0]  [7:6]=0 = x00 <- tmp[2][col=2][re]
    //  So Immediate for shufps is: 0x0001 1011 = x1B = 27
    v6.vector = _mm_shuffle_ps(v6.vector, v4.vector, 0x1B);

    // Does a v2 - v6*v7
    //  v2[0] <- v2[0] - v6[0]*v7[0] =  tmp[1][col=1][re] + tmp[2][col=1][im]
    //  v2[1] <- v2[1] - v6[1]*v7[1] =  tmp[1][col=1][im] - tmp[2][col=1][re]
    //  v2[2] <- v2[2] - v6[2]*v7[2] =  tmp[1][col=2][re] + tmp[2][col=2][im]
    //  v2[3] <- v2[3] - v6[3]*v7[3] =  tmp[1][col=2][im] - tmp[2][col=2][re]
    v6.vector = _mm_mul_ps(v6.vector, v7.vector);
    v2.vector = _mm_sub_ps(v2.vector, v6.vector);
    _mm_store_ps(dst_shadow+8, v2.vector);

    dst_shadow+=12;

    // Store in the next spinor -- should be prefetched if all has gone well
    v0.vector = _mm_load_ps(src_shadow);
    v1.vector = _mm_load_ps(src_shadow+4);
    v2.vector = _mm_load_ps(src_shadow+8);
    v3.vector = _mm_load_ps(src_shadow+12);
    v4.vector = _mm_load_ps(src_shadow+16);
    v5.vector = _mm_load_ps(src_shadow+20);

  }

  
  v6.vector = v4.vector;         // V6 is dest so we can move its 
  // I want to shuffle so that:
  //  v6[0] <- v6[3]=v4[3]  [1:0]=3 = x11 <- tmp[3][col=0][im]
  //  v6[1] <- v6[2]=v4[2]  [3:2]=2 = x10 <- tmp[3][col=0][re]
  //  v6[2] <-       v5[1]  [5:4]=1 = x01 <- tmp[3][col=1][im]
  //  v6[3] <-       v5[0]  [7:6]=0 = x00 <- tmp[3][col=1][re]
  //  So Immediate for shufps is: 0x0001 1011 = x1B = 27
  v6.vector = _mm_shuffle_ps(v6.vector, v5.vector, 0x1B);
  
  // Does a v0 - v6*v7
  //  v0[0] <- v0[0] - v6[0]*v7[0] =  tmp[0][col=0][re] + tmp[3][col=0][im]
  //  v0[1] <- v0[1] - v6[1]*v7[1] =  tmp[0][col=0][im] - tmp[3][col=0][re]
  //  v0[2] <- v0[2] - v6[2]*v7[2] =  tmp[0][col=1][re] + tmp[3][col=1][im]
  //  v0[3] <- v0[3] - v6[3]*v7[3] =  tmp[0][col=1][im] - tmp[3][col=1][re]    
  v6.vector = _mm_mul_ps(v6.vector, v7.vector);
  v0.vector = _mm_sub_ps(v0.vector, v6.vector);
  _mm_store_ps(dst_shadow, v0.vector);
  
  v6.vector = v5.vector;
  // I want to shuffle so that:
  //  v6[0] <- v6[3]=v5[3]  [1:0]=3 = x11 <- tmp[3][col=0][im]
  //  v6[1] <- v6[2]=v5[2]  [3:2]=2 = x10 <- tmp[3][col=0][re]
  //  v6[2] <-       v3[1]  [5:4]=1 = x01 <- tmp[3][col=1][im]
  //  v6[3] <-       v3[0]  [7:6]=0 = x00 <- tmp[3][col=1][re]
  //  So Immediate for shufps is: 0x0001 1011 = x1B = 27
  v6.vector = _mm_shuffle_ps(v6.vector, v3.vector, 0x1B);
  
  // Does a v1 - v6*v7
  //  v1[0] <- v1[0] - v6[0]*v7[0] =  tmp[0][col=2][re] + tmp[3][col=2][im]
  //  v1[1] <- v1[1] - v6[1]*v7[1] =  tmp[0][col=2][im] - tmp[3][col=2][re]
  //  v1[2] <- v1[2] - v6[2]*v7[2] =  tmp[1][col=0][re] + tmp[2][col=0][im]
  //  v1[3] <- v1[3] - v6[3]*v7[3] =  tmp[1][col=0][im] - tmp[2][col=0][re]
  v6.vector  = _mm_mul_ps(v6.vector, v7.vector);
  v1.vector  = _mm_sub_ps(v1.vector, v6.vector);
  _mm_store_ps(dst_shadow+4, v1.vector);
  
  v6.vector = v3.vector;
  //  v6[0] <- v6[3]=v3[3]  [1:0]=3 = x11 <- tmp[2][col=1][im]
  //  v6[1] <- v6[3]=v3[2]  [3:2]=2 = x10 <- tmp[2][col=1][re]
  //  v6[2] <-       v4[1]  [5:4]=1 = x01 <- tmp[2][col=2][im]
  //  v6[3] <-       v4[0]  [7:6]=0 = x00 <- tmp[2][col=2][re]
  //  So Immediate for shufps is: 0x0001 1011 = x1B = 27
  v6.vector = _mm_shuffle_ps(v6.vector, v4.vector, 0x1B);
  
  // Does a v2 - v6*v7
  //  v2[0] <- v2[0] - v6[0]*v7[0] =  tmp[1][col=1][re] + tmp[2][col=1][im]
  //  v2[1] <- v2[1] - v6[1]*v7[1] =  tmp[1][col=1][im] - tmp[2][col=1][re]
  //  v2[2] <- v2[2] - v6[2]*v7[2] =  tmp[1][col=2][re] + tmp[2][col=2][im]
  //  v2[3] <- v2[3] - v6[3]*v7[3] =  tmp[1][col=2][im] - tmp[2][col=2][re]
  v6.vector = _mm_mul_ps(v6.vector, v7.vector);
  v2.vector = _mm_sub_ps(v2.vector, v6.vector);
  _mm_store_ps(dst_shadow+8, v2.vector);

}



/** \brief Spin Project (1/2)(1+\gamma_1)
 *
 * \param src (pointer to 4 vector source)
 * \param dst (pointer to 2 vector dest)
 * \param n_vec (number of vectors to project)

 * It is assumeed that src points to an array of  4 spin
 * 3 colour vectors with complex components. 
 * It is assumed that dst points to an array of 2 spin
 * 3 colour vectors with complex components. In both cases
 * ordering is that spin is slowest, and complex components
 * are fastest. */
inline  
void inlineSpinProjDir1Plus(const REAL32* src, REAL32 *dst, unsigned int n_vec) 
{

#ifdef DEBUG_GENERIC_SPIN_PROJ_INLINES
  QDPIO::cout << "inlineSpinProjDir1Plus" << endl;
#endif

 /* 1 + \gamma_1 =  1  0  0 -1 
                     0  1  1  0
                     0  1  1  0
                    -1  0  0  1 
 
   *      ( b0r + i b0i )  =  ( {a0r - a3r} + i{a0i - a3i} )
   *      ( b1r + i b1i )     ( {a1r + a2r} + i{a1i + a2i} )
   */

  const REAL32* src_shadow = src;
  REAL32* dst_shadow = dst;
  SSEVec v0, v1, v2, v3, v4, v5, v6, v7;

  v7.floats[0] = -1;
  v7.floats[1] = -1;
  v7.floats[2] = +1;
  v7.floats[3] = +1;
  
  // Store in the spinor - top half
  v0.vector = _mm_load_ps(src_shadow);
  v1.vector = _mm_load_ps(src_shadow+4);
  v2.vector = _mm_load_ps(src_shadow+8);
  v3.vector = _mm_load_ps(src_shadow+12);
  v4.vector = _mm_load_ps(src_shadow+16);
  v5.vector = _mm_load_ps(src_shadow+20);

  for(unsigned int site=0; site < n_vec-1; site++) {
    src_shadow += 24;
    
    _mm_prefetch((const char *)src_shadow, _MM_HINT_T0);
    
    v6.vector = v4.vector;         // V6 is dest so we can move its 
    // I want to shuffle so that:
    //  v6[0] <- v6[2]=v4[2]  [1:0]=2 = x10 <- tmp[3][col=0][re]
    //  v6[1] <- v6[3]=v4[3]  [3:2]=3 = x11 <- tmp[3][col=0][im]
    //  v6[2] <-       v5[0]  [5:4]=0 = x00 <- tmp[3][col=1][re]
    //  v6[3] <-       v5[1]  [7:6]=1 = x01 <- tmp[3][col=1][im]
    //  So Immediate for shufps is: 0x0100 1110 = x4E
    v6.vector = _mm_shuffle_ps(v6.vector, v5.vector, 0x4E);
    
    // Does a v0 - v6 
    //  v0[0] <- v0[0] - v6[0] =  tmp[0][col=0][re] - tmp[3][col=0][re]
    //  v0[1] <- v0[1] - v6[1] =  tmp[0][col=0][im] - tmp[3][col=0][im]
    //  v0[2] <- v0[2] - v6[2] =  tmp[0][col=1][re] - tmp[3][col=1][re]
    //  v0[3] <- v0[3] - v6[3] =  tmp[0][col=1][im] - tmp[3][col=1][im]
    v0.vector = _mm_sub_ps(v0.vector, v6.vector);
    _mm_store_ps(dst_shadow, v0.vector);

    v6.vector = v5.vector;
    // Now setup v6 so that 
    //  v6[0] <- v6[2]=v5[2]  [1:0]=2 = x11 <- tmp[3][col=2][re]
    //  v6[1] <- v6[3]=v5[3]  [3:2]=3 = x10 <- tmp[3][col=2][im]
    //  v6[2] <-       v3[0]  [5:4]=0 = x01 <- tmp[2][col=0][re]
    //  v6[3] <-       v3[1]  [7:6]=1 = x00 <- tmp[2][col=0][im]
    //  So Immediate for shufps is: 0x0100 1110 = x4E
    v6.vector = _mm_shuffle_ps(v6.vector, v3.vector, 0x4E);

    // Does a v1 + v7*v6  
    //  v1[0] <- v1[0] + v7[0]*v6[0] =  tmp[0][col=2][re] - tmp[3][col=2][re]
    //  v1[1] <- v1[1] + v7[1]*v6[1] =  tmp[0][col=2][im] - tmp[3][col=2][im]
    //  v1[2] <- v1[2] + v7[2]*v6[2] =  tmp[1][col=0][re] + tmp[2][col=0][re]
    //  v1[3] <- v1[3] + v7[3]*v6[3] =  tmp[1][col=0][im] + tmp[2][col=0][im]
    v6.vector  = _mm_mul_ps(v6.vector, v7.vector);
    v1.vector  = _mm_add_ps(v1.vector, v6.vector);
    _mm_store_ps(dst_shadow+4, v1.vector);

    v6.vector = v3.vector;
    //  v6[0] <- v6[2]=v3[2]  [1:0]=2 = x10 <- tmp[2][col=1][re]
    //  v6[1] <- v6[3]=v3[3]  [3:2]=3 = x11 <- tmp[2][col=1][im]
    //  v6[2] <-       v4[0]  [5:4]=0 = x00 <- tmp[2][col=2][re]
    //  v6[3] <-       v4[1]  [7:6]=1 = x01 <- tmp[2][col=2][im]
    //  So Immediate for shufps is: 0x0100 1110 = x4E

    v6.vector = _mm_shuffle_ps(v6.vector, v4.vector, 0x4E);

    // Does a v2 + v6 
    //  v2[0] <- v2[0] + v6[0] =  tmp[1][col=1][re] + tmp[2][col=1][re]
    //  v2[1] <- v2[1] + v6[1] =  tmp[1][col=1][im] + tmp[2][col=1][im]
    //  v2[2] <- v2[2] + v6[2] =  tmp[1][col=2][re] + tmp[2][col=2][re]
    //  v2[3] <- v2[3] + v6[3] =  tmp[1][col=2][im] + tmp[2][col=2][im]
    v2.vector = _mm_add_ps(v2.vector, v6.vector);
    _mm_store_ps(dst_shadow+8, v2.vector);

    dst_shadow+=12;

    // Store in the next spinor -- should be prefetched if all has gone well
    v0.vector = _mm_load_ps(src_shadow);
    v1.vector = _mm_load_ps(src_shadow+4);
    v2.vector = _mm_load_ps(src_shadow+8);
    v3.vector = _mm_load_ps(src_shadow+12);
    v4.vector = _mm_load_ps(src_shadow+16);
    v5.vector = _mm_load_ps(src_shadow+20);

  }
  
  v6.vector = v4.vector;         // V6 is dest so we can move its 
  // I want to shuffle so that:
  //  v6[0] <- v6[2]=v4[2]  [1:0]=2 = x10 <- tmp[3][col=0][re]
  //  v6[1] <- v6[3]=v4[3]  [3:2]=3 = x11 <- tmp[3][col=0][im]
  //  v6[2] <-       v5[0]  [5:4]=0 = x00 <- tmp[3][col=1][re]
  //  v6[3] <-       v5[1]  [7:6]=1 = x01 <- tmp[3][col=1][im]
  //  So Immediate for shufps is: 0x0100 1110 = x4E
  v6.vector = _mm_shuffle_ps(v6.vector, v5.vector, 0x4E);
  
  // Does a v0 - v6 
  //  v0[0] <- v0[0] - v6[0] =  tmp[0][col=0][re] - tmp[3][col=0][re]
  //  v0[1] <- v0[1] - v6[1] =  tmp[0][col=0][im] - tmp[3][col=0][im]
  //  v0[2] <- v0[2] - v6[2] =  tmp[0][col=1][re] - tmp[3][col=1][re]
  //  v0[3] <- v0[3] - v6[3] =  tmp[0][col=1][im] - tmp[3][col=1][im]
  v0.vector = _mm_sub_ps(v0.vector, v6.vector);
  _mm_store_ps(dst_shadow, v0.vector);
  
  v6.vector = v5.vector;
  // Now setup v6 so that 
  //  v6[0] <- v6[2]=v5[2]  [1:0]=2 = x11 <- tmp[3][col=2][re]
  //  v6[1] <- v6[3]=v5[3]  [3:2]=3 = x10 <- tmp[3][col=2][im]
  //  v6[2] <-       v3[0]  [5:4]=0 = x01 <- tmp[2][col=0][re]
  //  v6[3] <-       v3[1]  [7:6]=1 = x00 <- tmp[2][col=0][im]
  //  So Immediate for shufps is: 0x0100 1110 = x4E
  v6.vector = _mm_shuffle_ps(v6.vector, v3.vector, 0x4E);
  
  // Does a v1 + v7*v6  
  //  v1[0] <- v1[0] + v7[0]*v6[0] =  tmp[0][col=2][re] - tmp[3][col=2][re]
  //  v1[1] <- v1[1] + v7[1]*v6[1] =  tmp[0][col=2][im] - tmp[3][col=2][im]
  //  v1[2] <- v1[2] + v7[2]*v6[2] =  tmp[1][col=0][re] + tmp[2][col=0][re]
  //  v1[3] <- v1[3] + v7[3]*v6[3] =  tmp[1][col=0][im] + tmp[2][col=0][im]
  v6.vector  = _mm_mul_ps(v6.vector, v7.vector);
  v1.vector  = _mm_add_ps(v1.vector, v6.vector);
  _mm_store_ps(dst_shadow+4, v1.vector);
  
  v6.vector = v3.vector;
  //  v6[0] <- v6[2]=v3[2]  [1:0]=2 = x10 <- tmp[2][col=1][re]
  //  v6[1] <- v6[3]=v3[3]  [3:2]=3 = x11 <- tmp[2][col=1][im]
  //  v6[2] <-       v4[0]  [5:4]=0 = x00 <- tmp[2][col=2][re]
  //  v6[3] <-       v4[1]  [7:6]=1 = x01 <- tmp[2][col=2][im]
  //  So Immediate for shufps is: 0x0100 1110 = x4E
  
  v6.vector = _mm_shuffle_ps(v6.vector, v4.vector, 0x4E);
  
  // Does a v2 + v6 
  //  v2[0] <- v2[0] + v6[0] =  tmp[1][col=1][re] + tmp[2][col=1][re]
  //  v2[1] <- v2[1] + v6[1] =  tmp[1][col=1][im] + tmp[2][col=1][im]
  //  v2[2] <- v2[2] + v6[2] =  tmp[1][col=2][re] + tmp[2][col=2][re]
  //  v2[3] <- v2[3] + v6[3] =  tmp[1][col=2][im] + tmp[2][col=2][im]
  v2.vector = _mm_add_ps(v2.vector, v6.vector);
  _mm_store_ps(dst_shadow+8, v2.vector);

}

/** \brief Spin Project (1/2)(1-\gamma_1)
 *
 * \param src (pointer to 4 vector source)
 * \param dst (pointer to 2 vector dest)
 *
 * It is assumeed that src points to an array of  4 spin
 * 3 colour vectors with complex components. 
 * It is assumed that dst points to an array of 2 spin
 * 3 colour vectors with complex components. In both cases
 * ordering is that spin is slowest, and complex components
 * are fastest. */
inline 
void inlineSpinProjDir1Minus(const REAL32* src, REAL32 *dst, unsigned int n_vec) 
{

#ifdef DEBUG_GENERIC_SPIN_PROJ_INLINES
  QDPIO::cout << "inlineSpinProjDir1Minus" << endl;
#endif

  /* 1 - \gamma_1 =  1  0  0 +1 
                     0  1 -1  0
                     0 -1  1  0
                    +1  0  0  1 
 
   *      ( b0r + i b0i )  =  ( {a0r + a3r} + i{a0i + a3i} )
   *      ( b1r + i b1i )     ( {a1r - a2r} + i{a1i - a2i} )
   */
  const REAL32* src_shadow = src;
  REAL32* dst_shadow = dst;
  SSEVec v0, v1, v2, v3, v4, v5, v6, v7;

  v7.floats[0] = -1;
  v7.floats[1] = -1;
  v7.floats[2] = +1;
  v7.floats[3] = +1;

  // Store in the spinor - top half
  v0.vector = _mm_load_ps(src_shadow);
  v1.vector = _mm_load_ps(src_shadow+4);
  v2.vector = _mm_load_ps(src_shadow+8);
  v3.vector = _mm_load_ps(src_shadow+12);
  v4.vector = _mm_load_ps(src_shadow+16);
  v5.vector = _mm_load_ps(src_shadow+20);
  
  for(unsigned int site=0; site < n_vec-1; site++) {
    src_shadow += 24;
    _mm_prefetch((const char *)src_shadow, _MM_HINT_T0);
    
    v6.vector = v4.vector;         // V6 is dest so we can move its 
    // I want to shuffle so that:
    //  v6[0] <- v6[2]=v4[2]  [1:0]=2 = x10 <- tmp[3][col=0][re]
    //  v6[1] <- v6[3]=v4[3]  [3:2]=3 = x11 <- tmp[3][col=0][im]
    //  v6[2] <-       v5[0]  [5:4]=0 = x00 <- tmp[3][col=1][re]
    //  v6[3] <-       v5[1]  [7:6]=1 = x01 <- tmp[3][col=1][im]
    //  So Immediate for shufps is: 0x0100 1110 = x4E
    v6.vector = _mm_shuffle_ps(v6.vector, v5.vector, 0x4E);
    
    // Does a v0 + v6
    //  v0[0] <- v0[0] + v6[0] =  tmp[0][col=0][re] + tmp[3][col=0][re]
    //  v0[1] <- v0[1] + v6[1] =  tmp[0][col=0][im] + tmp[3][col=0][im]
    //  v0[2] <- v0[2] + v6[2] =  tmp[0][col=1][re] + tmp[3][col=1][re]
    //  v0[3] <- v0[3] + v6[3] =  tmp[0][col=1][im] + tmp[3][col=1][im]
    v0.vector = _mm_add_ps(v0.vector, v6.vector);
    _mm_store_ps(dst_shadow, v0.vector);

    v6.vector = v5.vector;
    // Now setup v6 so that 
    //  v6[0] <- v6[2]=v5[2]  [1:0]=2 = x11 <- tmp[3][col=2][re]
    //  v6[1] <- v6[3]=v5[3]  [3:2]=3 = x10 <- tmp[3][col=2][im]
    //  v6[2] <-       v3[0]  [5:4]=0 = x01 <- tmp[2][col=0][re]
    //  v6[3] <-       v3[1]  [7:6]=1 = x00 <- tmp[2][col=0][im]
    //  So Immediate for shufps is: 0x0100 1110 = x4E
    v6.vector = _mm_shuffle_ps(v6.vector, v3.vector, 0x4E);


    // Does a v1 - v7*v6  
    //  v1[0] <- v1[0] - v7[0]*v6[0] =  tmp[0][col=2][re] + tmp[3][col=2][re]
    //  v1[1] <- v1[1] - v7[1]*v6[1] =  tmp[0][col=2][im] + tmp[3][col=2][im]
    //  v1[2] <- v1[2] - v7[2]*v6[2] =  tmp[1][col=0][re] - tmp[2][col=0][re]
    //  v1[3] <- v1[3] - v7[3]*v6[3] =  tmp[1][col=0][im] - tmp[2][col=0][im]
    v6.vector  = _mm_mul_ps(v6.vector, v7.vector);
    v1.vector  = _mm_sub_ps(v1.vector, v6.vector);
    _mm_store_ps(dst_shadow+4, v1.vector);

    v6.vector = v3.vector;
    //  v6[0] <- v6[2]=v3[2]  [1:0]=2 = x10 <- tmp[2][col=1][re]
    //  v6[1] <- v6[3]=v3[3]  [3:2]=3 = x11 <- tmp[2][col=1][im]
    //  v6[2] <-       v4[0]  [5:4]=0 = x00 <- tmp[2][col=2][re]
    //  v6[3] <-       v4[1]  [7:6]=1 = x01 <- tmp[2][col=2][im]
    //  So Immediate for shufps is: 0x0100 1110 = x4E

    v6.vector = _mm_shuffle_ps(v6.vector, v4.vector, 0x4E);

    // Does a v2 - v6
    //  v2[0] <- v2[0] - v6[0] =  tmp[1][col=1][re] - tmp[2][col=1][re]
    //  v2[1] <- v2[1] - v6[1] =  tmp[1][col=1][im] - tmp[2][col=1][im]
    //  v2[2] <- v2[2] - v6[2] =  tmp[1][col=2][re] - tmp[2][col=2][re]
    //  v2[3] <- v2[3] - v6[3] =  tmp[1][col=2][im] - tmp[2][col=2][im]
    v2.vector = _mm_sub_ps(v2.vector, v6.vector);
    _mm_store_ps(dst_shadow+8, v2.vector);

    dst_shadow+=12;

    v0.vector = _mm_load_ps(src_shadow);
    v1.vector = _mm_load_ps(src_shadow+4);
    v2.vector = _mm_load_ps(src_shadow+8);
    v3.vector = _mm_load_ps(src_shadow+12);
    v4.vector = _mm_load_ps(src_shadow+16);
    v5.vector = _mm_load_ps(src_shadow+20);
   
  }
  
  
  v6.vector = v4.vector;         // V6 is dest so we can move its 
  // I want to shuffle so that:
  //  v6[0] <- v6[2]=v4[2]  [1:0]=2 = x10 <- tmp[3][col=0][re]
  //  v6[1] <- v6[3]=v4[3]  [3:2]=3 = x11 <- tmp[3][col=0][im]
  //  v6[2] <-       v5[0]  [5:4]=0 = x00 <- tmp[3][col=1][re]
  //  v6[3] <-       v5[1]  [7:6]=1 = x01 <- tmp[3][col=1][im]
  //  So Immediate for shufps is: 0x0100 1110 = x4E
  v6.vector = _mm_shuffle_ps(v6.vector, v5.vector, 0x4E);
    
  // Does a v0 + v6
  //  v0[0] <- v0[0] + v6[0] =  tmp[0][col=0][re] + tmp[3][col=0][re]
  //  v0[1] <- v0[1] + v6[1] =  tmp[0][col=0][im] + tmp[3][col=0][im]
  //  v0[2] <- v0[2] + v6[2] =  tmp[0][col=1][re] + tmp[3][col=1][re]
  //  v0[3] <- v0[3] + v6[3] =  tmp[0][col=1][im] + tmp[3][col=1][im]
  v0.vector = _mm_add_ps(v0.vector, v6.vector);
  _mm_store_ps(dst_shadow, v0.vector);
  
  v6.vector = v5.vector;
  // Now setup v6 so that 
  //  v6[0] <- v6[2]=v5[2]  [1:0]=2 = x11 <- tmp[3][col=2][re]
  //  v6[1] <- v6[3]=v5[3]  [3:2]=3 = x10 <- tmp[3][col=2][im]
  //  v6[2] <-       v3[0]  [5:4]=0 = x01 <- tmp[2][col=0][re]
  //  v6[3] <-       v3[1]  [7:6]=1 = x00 <- tmp[2][col=0][im]
  //  So Immediate for shufps is: 0x0100 1110 = x4E
  v6.vector = _mm_shuffle_ps(v6.vector, v3.vector, 0x4E);
  
  
  // Does a v1 - v7*v6  
  //  v1[0] <- v1[0] - v7[0]*v6[0] =  tmp[0][col=2][re] + tmp[3][col=2][re]
  //  v1[1] <- v1[1] - v7[1]*v6[1] =  tmp[0][col=2][im] + tmp[3][col=2][im]
  //  v1[2] <- v1[2] - v7[2]*v6[2] =  tmp[1][col=0][re] - tmp[2][col=0][re]
  //  v1[3] <- v1[3] - v7[3]*v6[3] =  tmp[1][col=0][im] - tmp[2][col=0][im]
  v6.vector  = _mm_mul_ps(v6.vector, v7.vector);
  v1.vector  = _mm_sub_ps(v1.vector, v6.vector);
  _mm_store_ps(dst_shadow+4, v1.vector);
  
  v6.vector = v3.vector;
  //  v6[0] <- v6[2]=v3[2]  [1:0]=2 = x10 <- tmp[2][col=1][re]
  //  v6[1] <- v6[3]=v3[3]  [3:2]=3 = x11 <- tmp[2][col=1][im]
  //  v6[2] <-       v4[0]  [5:4]=0 = x00 <- tmp[2][col=2][re]
  //  v6[3] <-       v4[1]  [7:6]=1 = x01 <- tmp[2][col=2][im]
  //  So Immediate for shufps is: 0x0100 1110 = x4E
  
  v6.vector = _mm_shuffle_ps(v6.vector, v4.vector, 0x4E);
  
  // Does a v2 - v6
  //  v2[0] <- v2[0] - v6[0] =  tmp[1][col=1][re] - tmp[2][col=1][re]
  //  v2[1] <- v2[1] - v6[1] =  tmp[1][col=1][im] - tmp[2][col=1][im]
  //  v2[2] <- v2[2] - v6[2] =  tmp[1][col=2][re] - tmp[2][col=2][re]
  //  v2[3] <- v2[3] - v6[3] =  tmp[1][col=2][im] - tmp[2][col=2][im]
  v2.vector = _mm_sub_ps(v2.vector, v6.vector);
  _mm_store_ps(dst_shadow+8, v2.vector);
  
}


/** \brief Spin Project (1/2)(1+\gamma_2)
 *
 * \param src (pointer to 4 vector source)
 * \param dst (pointer to 2 vector dest)
 * \param n_vec (number of vectors to project)

 * It is assumeed that src points to an array of  4 spin
 * 3 colour vectors with complex components. 
 * It is assumed that dst points to an array of 2 spin
 * 3 colour vectors with complex components. In both cases
 * ordering is that spin is slowest, and complex components
 * are fastest. */
inline  
void inlineSpinProjDir2Plus(const REAL32* src, REAL32 *dst, unsigned int n_vec) 
{

#ifdef DEBUG_GENERIC_SPIN_PROJ_INLINES
  QDPIO::cout << "inlineSpinProjDir2Plus" << endl;
#endif
  /* 1 + \gamma_2 =  1  0  i  0 
                     0  1  0 -i
                    -i  0  1  0
                     0  i  0  1 


   *      ( b0r + i b0i )  =  ( {a0r - a2i} + i{a0i + a2r} )
   *      ( b1r + i b1i )     ( {a1r + a3i} + i{a1i - a3r} )
   */

  const REAL32* src_shadow = src;
  REAL32* dst_shadow = dst;
  SSEVec v0, v1, v2, v3, v4, v5, v6, v7;

  v7.floats[0] = -1;
  v7.floats[1] = +1;
  v7.floats[2] = -1;
  v7.floats[3] = +1;
  
  // The shuffling will not need extra registers 
  // so I can reuse v6 for a set of signs.
  v6.floats[0] = -1;
  v6.floats[1] = +1;
  v6.floats[2] = +1;
  v6.floats[3] = -1;
  
  // Store in the spinor - top half
  v0.vector = _mm_load_ps(src_shadow);
  v1.vector = _mm_load_ps(src_shadow+4);
  v2.vector = _mm_load_ps(src_shadow+8);
  v3.vector = _mm_load_ps(src_shadow+12);
  v4.vector = _mm_load_ps(src_shadow+16);
  v5.vector = _mm_load_ps(src_shadow+20);

  for(unsigned int site=0; site < n_vec-1; site++) {
    src_shadow += 24;
    _mm_prefetch((const char *)src_shadow, _MM_HINT_T0);
    
    // I want to shuffle so that
    //  v3[0] <- v3[1]  [1:0]=1 = x01 <- tmp[2][col=0][im]
    //  v3[1] <- v3[0]  [3:2]=0 = x00 <- tmp[2][col=0][re]
    //  v3[2] <- v3[3]  [5:4]=3 = x11 <- tmp[2][col=1][im]
    //  v3[3] <- v3[2]  [7:6]=2 = x10 <- tmp[2][col=1][re]
    //
    // Note I don't need an extra  here. Shuf from v3 to 
    // v3 directly.
    //
    //  So Immediate for shufps is: 0x1011 0001 = xB1
    v3.vector = _mm_shuffle_ps(v3.vector, v3.vector, 0xB1);
    
    // Does a v0 + v7*v3 
    //  v0[0] <- v0[0] + v7[0]*v3[0] =  tmp[0][col=0][re] - tmp[2][col=0][im]
    //  v0[1] <- v0[1] + v7[1]*v3[1] =  tmp[0][col=0][im] + tmp[2][col=0][re]
    //  v0[2] <- v0[2] + v7[2]*v3[2] =  tmp[0][col=1][re] - tmp[2][col=1][im]
    //  v0[3] <- v0[3] + v7[3]*v3[3] =  tmp[0][col=1][im] + tmp[2][col=1][re]
    v3.vector = _mm_mul_ps(v3.vector, v7.vector);
    v0.vector = _mm_add_ps(v0.vector, v3.vector);
    _mm_store_ps(dst_shadow, v0.vector);

    // Now setup v5 so that 
    //  v4[0] <- v4[1]  [1:0]=1 = x01 <- tmp[2][col=2][im]
    //  v4[1] <- v4[0]  [3:2]=0 = x00 <- tmp[2][col=2][re]
    //  v4[2] <- v4[3]  [5:4]=3 = x11 <- tmp[3][col=0][im]
    //  v4[3] <- v4[2]  [7:6]=2 = x10 <- tmp[3][col=0][re]
    //  So Immediate for shufps is: 0x1011 0001 = xB1
    v4.vector = _mm_shuffle_ps(v4.vector, v4.vector, 0xB1);
    
    // Does a v1 + v6*v4 
    //  v1[0] <- v1[0] + v6[0]*v4[0] =  tmp[0][col=2][re] - tmp[2][col=2][im]
    //  v1[1] <- v1[1] + v6[1]*v4[1] =  tmp[0][col=2][im] + tmp[2][col=2][re]
    //  v1[2] <- v1[2] + v6[2]*v4[2] =  tmp[1][col=0][re] + tmp[3][col=0][im]
    //  v1[3] <- v1[3] + v6[3]*v4[3] =  tmp[1][col=0][im] - tmp[3][col=0][re]
    v4.vector  = _mm_mul_ps(v4.vector, v6.vector);
    v1.vector  = _mm_add_ps(v1.vector, v4.vector);
    _mm_store_ps(dst_shadow+4, v1.vector);

    //  v5[0] <- v5[1]  [1:0]=1 = x01 <- tmp[3][col=1][im]
    //  v5[1] <- v5[0]  [3:2]=0 = x00 <- tmp[3][col=1][re]
    //  v5[2] <- v5[3]  [5:4]=3 = x11 <- tmp[3][col=2][im]
    //  v5[3] <- v5[2]  [7:6]=2 = x10 <- tmp[3][col=2][re]
    //  So Immediate for shufps is: 0x1011 0001 = xB1
    v5.vector = _mm_shuffle_ps(v5.vector, v5.vector, 0xB1);

    // Does a v0 - v7*v5 
    //  v2[0] <- v2[0] - v7[0]*v5[0] =  tmp[1][col=1][re] + tmp[2][col=1][im]
    //  v2[1] <- v2[1] - v7[1]*v5[1] =  tmp[1][col=1][im] - tmp[2][col=1][re]
    //  v2[2] <- v2[2] - v7[2]*v5[2] =  tmp[1][col=2][re] + tmp[2][col=2][im]
    //  v2[3] <- v2[3] - v7[3]*v5[3] =  tmp[1][col=2][im] - tmp[2][col=2][re]
    v5.vector = _mm_mul_ps(v5.vector, v7.vector );
    v2.vector = _mm_sub_ps(v2.vector, v5.vector);
    _mm_store_ps(dst_shadow+8, v2.vector);

    dst_shadow+=12;

    // Store in the next spinor - top half
    v0.vector = _mm_load_ps(src_shadow);
    v1.vector = _mm_load_ps(src_shadow+4);
    v2.vector = _mm_load_ps(src_shadow+8);
    v3.vector = _mm_load_ps(src_shadow+12);
    v4.vector = _mm_load_ps(src_shadow+16);
    v5.vector = _mm_load_ps(src_shadow+20);
   
  }
  
  // I want to shuffle so that
  //  v3[0] <- v3[1]  [1:0]=1 = x01 <- tmp[2][col=0][im]
  //  v3[1] <- v3[0]  [3:2]=0 = x00 <- tmp[2][col=0][re]
  //  v3[2] <- v3[3]  [5:4]=3 = x11 <- tmp[2][col=1][im]
  //  v3[3] <- v3[2]  [7:6]=2 = x10 <- tmp[2][col=1][re]
  //
  // Note I don't need an extra  here. Shuf from v3 to 
  // v3 directly.
  //
  //  So Immediate for shufps is: 0x1011 0001 = xB1
  v3.vector = _mm_shuffle_ps(v3.vector, v3.vector, 0xB1);
  
  // Does a v0 + v7*v3 
  //  v0[0] <- v0[0] + v7[0]*v3[0] =  tmp[0][col=0][re] - tmp[2][col=0][im]
  //  v0[1] <- v0[1] + v7[1]*v3[1] =  tmp[0][col=0][im] + tmp[2][col=0][re]
  //  v0[2] <- v0[2] + v7[2]*v3[2] =  tmp[0][col=1][re] - tmp[2][col=1][im]
  //  v0[3] <- v0[3] + v7[3]*v3[3] =  tmp[0][col=1][im] + tmp[2][col=1][re]
  v3.vector = _mm_mul_ps(v3.vector, v7.vector);
  v0.vector = _mm_add_ps(v0.vector, v3.vector);
  _mm_store_ps(dst_shadow, v0.vector);
  
  // Now setup v5 so that 
  //  v4[0] <- v4[1]  [1:0]=1 = x01 <- tmp[2][col=2][im]
  //  v4[1] <- v4[0]  [3:2]=0 = x00 <- tmp[2][col=2][re]
  //  v4[2] <- v4[3]  [5:4]=3 = x11 <- tmp[3][col=0][im]
  //  v4[3] <- v4[2]  [7:6]=2 = x10 <- tmp[3][col=0][re]
  //  So Immediate for shufps is: 0x1011 0001 = xB1
  v4.vector = _mm_shuffle_ps(v4.vector, v4.vector, 0xB1);
  
  // Does a v1 + v6*v4 
  //  v1[0] <- v1[0] + v6[0]*v4[0] =  tmp[0][col=2][re] - tmp[2][col=2][im]
  //  v1[1] <- v1[1] + v6[1]*v4[1] =  tmp[0][col=2][im] + tmp[2][col=2][re]
  //  v1[2] <- v1[2] + v6[2]*v4[2] =  tmp[1][col=0][re] + tmp[3][col=0][im]
  //  v1[3] <- v1[3] + v6[3]*v4[3] =  tmp[1][col=0][im] - tmp[3][col=0][re]
  v4.vector  = _mm_mul_ps(v4.vector, v6.vector);
  v1.vector  = _mm_add_ps(v1.vector, v4.vector);
  _mm_store_ps(dst_shadow+4, v1.vector);
  
  //  v5[0] <- v5[1]  [1:0]=1 = x01 <- tmp[3][col=1][im]
  //  v5[1] <- v5[0]  [3:2]=0 = x00 <- tmp[3][col=1][re]
  //  v5[2] <- v5[3]  [5:4]=3 = x11 <- tmp[3][col=2][im]
  //  v5[3] <- v5[2]  [7:6]=2 = x10 <- tmp[3][col=2][re]
  //  So Immediate for shufps is: 0x1011 0001 = xB1
  v5.vector = _mm_shuffle_ps(v5.vector, v5.vector, 0xB1);
  
  // Does a v0 - v7*v5 
  //  v2[0] <- v2[0] - v7[0]*v5[0] =  tmp[1][col=1][re] + tmp[2][col=1][im]
  //  v2[1] <- v2[1] - v7[1]*v5[1] =  tmp[1][col=1][im] - tmp[2][col=1][re]
  //  v2[2] <- v2[2] - v7[2]*v5[2] =  tmp[1][col=2][re] + tmp[2][col=2][im]
  //  v2[3] <- v2[3] - v7[3]*v5[3] =  tmp[1][col=2][im] - tmp[2][col=2][re]
  v5.vector = _mm_mul_ps(v5.vector, v7.vector );
  v2.vector = _mm_sub_ps(v2.vector, v5.vector);
  _mm_store_ps(dst_shadow+8, v2.vector);
  
}

/** \brief Spin Project (1/2)(1-\gamma_2)
 *
 * \param src (pointer to 4 vector source)
 * \param dst (pointer to 2 vector dest)
 *
 * It is assumeed that src points to an array of  4 spin
 * 3 colour vectors with complex components. 
 * It is assumed that dst points to an array of 2 spin
 * 3 colour vectors with complex components. In both cases
 * ordering is that spin is slowest, and complex components
 * are fastest. */
inline 
void inlineSpinProjDir2Minus(const REAL32* src, REAL32 *dst, unsigned int n_vec) 
{

#ifdef DEBUG_GENERIC_SPIN_PROJ_INLINES
  QDPIO::cout << "inlineSpinProjDir2Minus" << endl;
#endif

   /* 1 - \gamma_2 =  1  0  -i  0 
                     0  1  0  +i
                    +i  0  1   0
                     0 -i  0   1 


   *      ( b0r + i b0i )  =  ( {a0r + a2i} + i{a0i - a2r} )
   *      ( b1r + i b1i )     ( {a1r - a3i} + i{a1i + a3r} )
   */

  const REAL32* src_shadow = src;
  REAL32* dst_shadow = dst;
  SSEVec v0, v1, v2, v3, v4, v5, v6, v7;

  v7.floats[0] = -1;
  v7.floats[1] = +1;
  v7.floats[2] = -1;
  v7.floats[3] = +1;
  
  // The shuffling will not need extra registers 
  // so I can reuse v6 for a set of signs.
  v6.floats[0] = -1;
  v6.floats[1] = +1;
  v6.floats[2] = +1;
  v6.floats[3] = -1;

  // Store in the spinor - top half
  v0.vector = _mm_load_ps(src_shadow);
  v1.vector = _mm_load_ps(src_shadow+4);
  v2.vector = _mm_load_ps(src_shadow+8);
  v3.vector = _mm_load_ps(src_shadow+12);
  v4.vector = _mm_load_ps(src_shadow+16);
  v5.vector = _mm_load_ps(src_shadow+20);

  for(unsigned int site=0; site < n_vec-1; site++) {
    src_shadow += 24;
    _mm_prefetch((const char *)src_shadow, _MM_HINT_T0);
    
    // I want to shuffle so that
    //  v3[0] <- v3[1]  [1:0]=1 = x01 <- tmp[2][col=0][im]
    //  v3[1] <- v3[0]  [3:2]=0 = x00 <- tmp[2][col=0][re]
    //  v3[2] <- v3[3]  [5:4]=3 = x11 <- tmp[2][col=1][im]
    //  v3[3] <- v3[2]  [7:6]=2 = x10 <- tmp[2][col=1][re]
    //
    // Note I don't need an extra  here. Shuf from v3 to 
    // v3 directly.
    //
    //  So Immediate for shufps is: 0x1011 0001 = xB1
    v3.vector = _mm_shuffle_ps(v3.vector, v3.vector, 0xB1);
    
    // Does a v0 - v7*v3 
    //  v0[0] <- v0[0] - v7[0]*v3[0] =  tmp[0][col=0][re] + tmp[2][col=0][im]
    //  v0[1] <- v0[1] - v7[1]*v3[1] =  tmp[0][col=0][im] - tmp[2][col=0][re]
    //  v0[2] <- v0[2] - v7[2]*v3[2] =  tmp[0][col=1][re] + tmp[2][col=1][im]
    //  v0[3] <- v0[3] - v7[3]*v3[3] =  tmp[0][col=1][im] - tmp[2][col=1][re]
    v3.vector = _mm_mul_ps(v3.vector, v7.vector);
    v0.vector = _mm_sub_ps(v0.vector, v3.vector);
    _mm_store_ps(dst_shadow, v0.vector);

    // Now setup v4 so that 
    //  v4[0] <- v4[1]  [1:0]=1 = x01 <- tmp[2][col=2][im]
    //  v4[1] <- v4[0]  [3:2]=0 = x00 <- tmp[2][col=2][re]
    //  v4[2] <- v4[3]  [5:4]=3 = x11 <- tmp[3][col=0][im]
    //  v4[3] <- v4[2]  [7:6]=2 = x10 <- tmp[3][col=0][re]
    //  So Immediate for shufps is: 0x1011 0001 = xB1
    v4.vector = _mm_shuffle_ps(v4.vector, v4.vector, 0xB1);
    
    // Does a v1 - v6*v4 
    //  v1[0] <- v1[0] + v6[0]*v4[0] =  tmp[0][col=2][re] + tmp[2][col=2][im]
    //  v1[1] <- v1[1] + v6[1]*v4[1] =  tmp[0][col=2][im] - tmp[2][col=2][re]
    //  v1[2] <- v1[2] + v6[2]*v4[2] =  tmp[1][col=0][re] - tmp[3][col=0][im]
    //  v1[3] <- v1[3] + v6[3]*v4[3] =  tmp[1][col=0][im] + tmp[3][col=0][re]
    v4.vector  = _mm_mul_ps(v4.vector, v6.vector);
    v1.vector  = _mm_sub_ps(v1.vector, v4.vector);
    _mm_store_ps(dst_shadow+4, v1.vector);

    //  v5[0] <- v5[1]  [1:0]=1 = x01 <- tmp[3][col=1][im]
    //  v5[1] <- v5[0]  [3:2]=0 = x00 <- tmp[3][col=1][re]
    //  v5[2] <- v5[3]  [5:4]=3 = x11 <- tmp[3][col=2][im]
    //  v5[3] <- v5[2]  [7:6]=2 = x10 <- tmp[3][col=2][re]
    //  So Immediate for shufps is: 0x1011 0001 = xB1
    v5.vector = _mm_shuffle_ps(v5.vector, v5.vector, 0xB1);

    // Does a v0 + v7*v5 
    //  v2[0] <- v2[0] + v7[0]*v5[0] =  tmp[1][col=1][re] - tmp[2][col=1][im]
    //  v2[1] <- v2[1] + v7[1]*v5[1] =  tmp[1][col=1][im] + tmp[2][col=1][re]
    //  v2[2] <- v2[2] + v7[2]*v5[2] =  tmp[1][col=2][re] - tmp[2][col=2][im]
    //  v2[3] <- v2[3] + v7[3]*v5[3] =  tmp[1][col=2][im] + tmp[2][col=2][re]
    v5.vector = _mm_mul_ps(v5.vector, v7.vector );
    v2.vector = _mm_add_ps(v2.vector, v5.vector);
    _mm_store_ps(dst_shadow+8, v2.vector);

    dst_shadow+=12;

    // Store in the spinor - top half
    v0.vector = _mm_load_ps(src_shadow);
    v1.vector = _mm_load_ps(src_shadow+4);
    v2.vector = _mm_load_ps(src_shadow+8);
    v3.vector = _mm_load_ps(src_shadow+12);
    v4.vector = _mm_load_ps(src_shadow+16);
    v5.vector = _mm_load_ps(src_shadow+20);
   
  }
  
  // I want to shuffle so that
  //  v3[0] <- v3[1]  [1:0]=1 = x01 <- tmp[2][col=0][im]
  //  v3[1] <- v3[0]  [3:2]=0 = x00 <- tmp[2][col=0][re]
  //  v3[2] <- v3[3]  [5:4]=3 = x11 <- tmp[2][col=1][im]
  //  v3[3] <- v3[2]  [7:6]=2 = x10 <- tmp[2][col=1][re]
  //
  // Note I don't need an extra  here. Shuf from v3 to 
  // v3 directly.
  //
  //  So Immediate for shufps is: 0x1011 0001 = xB1
  v3.vector = _mm_shuffle_ps(v3.vector, v3.vector, 0xB1);
  
  // Does a v0 - v7*v3 
  //  v0[0] <- v0[0] - v7[0]*v3[0] =  tmp[0][col=0][re] + tmp[2][col=0][im]
  //  v0[1] <- v0[1] - v7[1]*v3[1] =  tmp[0][col=0][im] - tmp[2][col=0][re]
  //  v0[2] <- v0[2] - v7[2]*v3[2] =  tmp[0][col=1][re] + tmp[2][col=1][im]
  //  v0[3] <- v0[3] - v7[3]*v3[3] =  tmp[0][col=1][im] - tmp[2][col=1][re]
  v3.vector = _mm_mul_ps(v3.vector, v7.vector);
  v0.vector = _mm_sub_ps(v0.vector, v3.vector);
  _mm_store_ps(dst_shadow, v0.vector);
  
  // Now setup v4 so that 
  //  v4[0] <- v4[1]  [1:0]=1 = x01 <- tmp[2][col=2][im]
  //  v4[1] <- v4[0]  [3:2]=0 = x00 <- tmp[2][col=2][re]
  //  v4[2] <- v4[3]  [5:4]=3 = x11 <- tmp[3][col=0][im]
  //  v4[3] <- v4[2]  [7:6]=2 = x10 <- tmp[3][col=0][re]
  //  So Immediate for shufps is: 0x1011 0001 = xB1
  v4.vector = _mm_shuffle_ps(v4.vector, v4.vector, 0xB1);
  
  // Does a v1 - v6*v4 
  //  v1[0] <- v1[0] + v6[0]*v4[0] =  tmp[0][col=2][re] + tmp[2][col=2][im]
  //  v1[1] <- v1[1] + v6[1]*v4[1] =  tmp[0][col=2][im] - tmp[2][col=2][re]
  //  v1[2] <- v1[2] + v6[2]*v4[2] =  tmp[1][col=0][re] - tmp[3][col=0][im]
  //  v1[3] <- v1[3] + v6[3]*v4[3] =  tmp[1][col=0][im] + tmp[3][col=0][re]
  v4.vector  = _mm_mul_ps(v4.vector, v6.vector);
  v1.vector  = _mm_sub_ps(v1.vector, v4.vector);
  _mm_store_ps(dst_shadow+4, v1.vector);
  
  //  v5[0] <- v5[1]  [1:0]=1 = x01 <- tmp[3][col=1][im]
  //  v5[1] <- v5[0]  [3:2]=0 = x00 <- tmp[3][col=1][re]
  //  v5[2] <- v5[3]  [5:4]=3 = x11 <- tmp[3][col=2][im]
  //  v5[3] <- v5[2]  [7:6]=2 = x10 <- tmp[3][col=2][re]
  //  So Immediate for shufps is: 0x1011 0001 = xB1
  v5.vector = _mm_shuffle_ps(v5.vector, v5.vector, 0xB1);
  
  // Does a v0 + v7*v5 
  //  v2[0] <- v2[0] + v7[0]*v5[0] =  tmp[1][col=1][re] - tmp[2][col=1][im]
  //  v2[1] <- v2[1] + v7[1]*v5[1] =  tmp[1][col=1][im] + tmp[2][col=1][re]
  //  v2[2] <- v2[2] + v7[2]*v5[2] =  tmp[1][col=2][re] - tmp[2][col=2][im]
  //  v2[3] <- v2[3] + v7[3]*v5[3] =  tmp[1][col=2][im] + tmp[2][col=2][re]
  v5.vector = _mm_mul_ps(v5.vector, v7.vector );
  v2.vector = _mm_add_ps(v2.vector, v5.vector);
  _mm_store_ps(dst_shadow+8, v2.vector);
 
}

/** \brief Spin Project (1/2)(1+\gamma_3)
 *
 * \param src (pointer to 4 vector source)
 * \param dst (pointer to 2 vector dest)
 * \param n_vec (number of vectors to project)

 * It is assumeed that src points to an array of  4 spin
 * 3 colour vectors with complex components. 
 * It is assumed that dst points to an array of 2 spin
 * 3 colour vectors with complex components. In both cases
 * ordering is that spin is slowest, and complex components
 * are fastest. */
inline  
void inlineSpinProjDir3Plus(const REAL32* src, REAL32 *dst, unsigned int n_vec) 
{

#ifdef DEBUG_GENERIC_SPIN_PROJ_INLINES
  QDPIO::cout << "inlineSpinProjDir3Plus" << endl;
#endif
  /* 1 + \gamma_3 =  1  0  1  0 
                     0  1  0  1
                     1  0  1  0
                     0  1  0  1 

   *      ( b0r + i b0i )  =  ( {a0r + a2r} + i{a0i + a2i} )
   *      ( b1r + i b1i )     ( {a1r + a3r} + i{a1i + a3i} )
   */
  const REAL32* src_shadow = src;
  REAL32* dst_shadow = dst;
  SSEVec v0, v1, v2, v3, v4, v5;

  // Store in the spinor - top half
  v0.vector = _mm_load_ps(src_shadow);
  v1.vector = _mm_load_ps(src_shadow+4);
  v2.vector = _mm_load_ps(src_shadow+8);
  v3.vector = _mm_load_ps(src_shadow+12);
  v4.vector = _mm_load_ps(src_shadow+16);
  v5.vector = _mm_load_ps(src_shadow+20);

  for(unsigned int site=0; site < n_vec-1; site++) {
    src_shadow += 24;
    _mm_prefetch((const char *)src_shadow, _MM_HINT_T0);

    // This is easy - I don't need any shuffling
    v0.vector = _mm_add_ps(v0.vector, v3.vector);
    v1.vector = _mm_add_ps(v1.vector, v4.vector);
    v2.vector = _mm_add_ps(v2.vector, v5.vector);

    _mm_store_ps(dst_shadow, v0.vector);
    _mm_store_ps(dst_shadow+4, v1.vector);
    _mm_store_ps(dst_shadow+8, v2.vector);

    dst_shadow+=12;

    v0.vector = _mm_load_ps(src_shadow);
    v1.vector = _mm_load_ps(src_shadow+4);
    v2.vector = _mm_load_ps(src_shadow+8);
    v3.vector = _mm_load_ps(src_shadow+12);
    v4.vector = _mm_load_ps(src_shadow+16);
    v5.vector = _mm_load_ps(src_shadow+20);
   
  }
  
  
  // This is easy - I don't need any shuffling
  v0.vector = _mm_add_ps(v0.vector, v3.vector);
  v1.vector = _mm_add_ps(v1.vector, v4.vector);
  v2.vector = _mm_add_ps(v2.vector, v5.vector);
  
  _mm_store_ps(dst_shadow, v0.vector);
  _mm_store_ps(dst_shadow+4, v1.vector);
  _mm_store_ps(dst_shadow+8, v2.vector);
  
}

/** \brief Spin Project (1/2)(1-\gamma_3)
 *
 * \param src (pointer to 4 vector source)
 * \param dst (pointer to 2 vector dest)
 *
 * It is assumeed that src points to an array of  4 spin
 * 3 colour vectors with complex components. 
 * It is assumed that dst points to an array of 2 spin
 * 3 colour vectors with complex components. In both cases
 * ordering is that spin is slowest, and complex components
 * are fastest. */
inline 
void inlineSpinProjDir3Minus(const REAL32* src, REAL32 *dst, unsigned int n_vec) 
{

#ifdef DEBUG_GENERIC_SPIN_PROJ_INLINES
  QDPIO::cout << "inlineSpinProjDir3Minus" << endl;
#endif

  /* 1 - \gamma_3 =  1  0  -1  0 
                     0  1  0  -1
                    -1  0  1  0
                     0 -1  0  1 

   *      ( b0r + i b0i )  =  ( {a0r - a2r} + i{a0i - a2i} )
   *      ( b1r + i b1i )     ( {a1r - a3r} + i{a1i - a3i} )
   */

  const REAL32* src_shadow = src;
  REAL32* dst_shadow = dst;
  SSEVec v0, v1, v2, v3, v4, v5;

  // Store in the spinor - top half
  v0.vector = _mm_load_ps(src_shadow);
  v1.vector = _mm_load_ps(src_shadow+4);
  v2.vector = _mm_load_ps(src_shadow+8);
  v3.vector = _mm_load_ps(src_shadow+12);
  v4.vector = _mm_load_ps(src_shadow+16);
  v5.vector = _mm_load_ps(src_shadow+20);

  for(unsigned int site=0; site < n_vec-1; site++) {
    src_shadow += 24;
    _mm_prefetch((const char *)src_shadow, _MM_HINT_T0);

    // This is easy - I don't need any shuffling
    v0.vector = _mm_sub_ps(v0.vector, v3.vector);
    v1.vector = _mm_sub_ps(v1.vector, v4.vector);
    v2.vector = _mm_sub_ps(v2.vector, v5.vector);

    _mm_store_ps(dst_shadow, v0.vector);
    _mm_store_ps(dst_shadow+4, v1.vector);
    _mm_store_ps(dst_shadow+8, v2.vector);

    dst_shadow+=12;

    v0.vector = _mm_load_ps(src_shadow);
    v1.vector = _mm_load_ps(src_shadow+4);
    v2.vector = _mm_load_ps(src_shadow+8);
    v3.vector = _mm_load_ps(src_shadow+12);
    v4.vector = _mm_load_ps(src_shadow+16);
    v5.vector = _mm_load_ps(src_shadow+20);
  }

  
  // This is easy - I don't need any shuffling
  v0.vector = _mm_sub_ps(v0.vector, v3.vector);
  v1.vector = _mm_sub_ps(v1.vector, v4.vector);
  v2.vector = _mm_sub_ps(v2.vector, v5.vector);
  
  _mm_store_ps(dst_shadow, v0.vector);
  _mm_store_ps(dst_shadow+4, v1.vector);
  _mm_store_ps(dst_shadow+8, v2.vector);
  
}

} // namespace QDP;

#endif 
