#ifndef SSE_SPIN_RECON_INLINES_H
#define SSE_SPIN_RECON_INLINES_H

#include "qdp_sse_intrin.h"

/* File: generic_spin_recon_inlines.h
   Purpose: Supply inline functions to do spin reconstruction
   Author: $Id: sse_spin_recon_inlines.h,v 1.6 2007-08-31 14:41:18 bjoo Exp $
*/
namespace QDP {

/** \brief Spin recon (1/2)(1+\gamma_0)
 *
 * \param src (pointer to 2 vector source)
 * \param dst (pointer to 4 vector dest)
 * \param n_vec (number of vectors to reconstruct)

 * It is assumeed that src points to an array of  2 spin
 * 3 colour vectors with complex components. 
 * It is assumed that dst points to an array of 4 spin
 * 3 colour vectors with complex components. In both cases
 * ordering is that spin is slowest, and complex components
 * are fastest. */
inline  
void inlineSpinReconDir0Plus(const REAL32* src, REAL32 *dst, unsigned int n_vec) 
{

#ifdef DEBUG_GENERIC_SPIN_RECON_INLINES
  QDPIO::cout << "inlineSpinProjDir0Plus" << endl;
#endif

  /* 1 + \gamma_0 =  1  0  0  i 
                     0  1  i  0
                     0 -i  1  0
                    -i  0  0  1 
 
    *  ( b2r + i b2i )  =  ( {a2r + a1i} + i{a2i - a1r} )  =  ( b1i - i b1r )
    *  ( b3r + i b3i )     ( {a3r + a0i} + i{a3i - a0r} )     ( b0i - i b0r ) 
   */

  // 12 components in source, 24 in dest
  // 3 vectors for source

  const REAL32* src_shadow = src;
  REAL32* dst_shadow = dst;

  
  SSEVec v0, v1, v2, v3, v4, v5, v6;
  v6.floats[0] = +1;
  v6.floats[1] = -1;
  v6.floats[2] = +1;
  v6.floats[3] = -1;

  v0.vector = _mm_load_ps(src_shadow);
  v1.vector = _mm_load_ps(src_shadow+4);
  v2.vector = _mm_load_ps(src_shadow+8);
  
  for(unsigned int site=0; site < n_vec-1; site++) { 

    src_shadow += 12;
    _mm_prefetch((const char *)src_shadow, _MM_HINT_T0);

    // Now we have
    // v0.floats[0] = tmp[0][col=0][re]
    // v0.floats[1] = tmp[0][col=0][im]
    // v0.floats[2] = tmp[0][col=1][re]
    // v0.floats[3] = tmp[0][col=1][im]
    //
    // v1.floats[0] = tmp[0][col=2][re]
    // v1.floats[1] = tmp[0][col=2][im]
    // v1.floats[2] = tmp[1][col=0][re]
    // v1.floats[3] = tmp[1][col=0][im]
    //
    // v2.floats[0] = tmp[1][col=1][re]
    // v2.floats[1] = tmp[1][col=1][im]
    // v2.floats[2] = tmp[1][col=2][re]
    // v2.floats[3] = tmp[1][col=2][im]
    
    // These are the top 2 components of the result
    // so I can store them out already.
    _mm_store_ps(dst_shadow, v0.vector);
    _mm_store_ps(dst_shadow+4, v1.vector);
    _mm_store_ps(dst_shadow+8, v2.vector);

    // I want to set up 
    //
    //     v3[0] = v1[3] = tmp[1][col=0][im]
    //     v3[1] = v1[2] = tmp[1][col=0][re]
    //     v3[2] = v2[1] = tmp[1][col=1][im]
    //     v3[3] = v2[0] = tmp[1][col=1][re]
    // 
    // We can do this with v3 = v1,  v3 = shuf( v3, v1 ). Shuf code is 00 01 10 11 = x1B
    v3.vector = v1.vector;
    v3.vector = _mm_shuffle_ps(v3.vector, v2.vector, 0x1B);
    
    // Now I need to do multiply the mask (+1, -1, +1, -1)
    v3.vector = _mm_mul_ps(v3.vector, v6.vector);

    // I want to set up 
    //
    //     v4[0] = v2[3] = tmp[1][col=2][im]
    //     v4[1] = v2[2] = tmp[1][col=2][re]
    //     v4[2] = v0[1] = tmp[0][col=0][im]
    //     v4[3] = v0[0] = tmp[0][col=0][re]
    // 
    // We can do this with v4 = v2,  v4 = shuf( v4, v2 ). Shuf code is 00 01 10 11 = x1B
    v4.vector = v2.vector;
    v4.vector = _mm_shuffle_ps(v4.vector, v0.vector, 0x1B);
    
    // Now I need to do multiply the mask (+1, -1, +1, -1)
    v4.vector = _mm_mul_ps(v4.vector, v6.vector);

    // I want to set up 
    //
    //     v5[0] = v0[3] = tmp[0][col=1][im]
    //     v5[1] = v0[2] = tmp[0][col=1][re]
    //     v5[2] = v1[1] = tmp[0][col=2][im]
    //     v5[3] = v1[0] = tmp[0][col=2][re]
    // 
    // We can do this with v5 = v0,  v5 = shuf( v5, v1 ). Shuf code is 00 01 10 11 = x1B
    v5.vector = v0.vector;
    v5.vector = _mm_shuffle_ps(v5.vector, v1.vector, 0x1B);
    
    // Now I need to do multiply the mask (+1, -1, +1, -1)
    v5.vector = _mm_mul_ps(v5.vector, v6.vector);

    // so I can  them out already - bypass cache.
    _mm_store_ps(dst_shadow+12, v3.vector);
    _mm_store_ps(dst_shadow+16, v4.vector);
    _mm_store_ps(dst_shadow+20, v5.vector);

    dst_shadow+=24;

    v0.vector = _mm_load_ps(src_shadow);
    v1.vector = _mm_load_ps(src_shadow+4);
    v2.vector = _mm_load_ps(src_shadow+8);

  }

  // Now we have
  // v0.floats[0] = tmp[0][col=0][re]
  // v0.floats[1] = tmp[0][col=0][im]
  // v0.floats[2] = tmp[0][col=1][re]
  // v0.floats[3] = tmp[0][col=1][im]
  //
  // v1.floats[0] = tmp[0][col=2][re]
  // v1.floats[1] = tmp[0][col=2][im]
  // v1.floats[2] = tmp[1][col=0][re]
  // v1.floats[3] = tmp[1][col=0][im]
  //
  // v2.floats[0] = tmp[1][col=1][re]
  // v2.floats[1] = tmp[1][col=1][im]
  // v2.floats[2] = tmp[1][col=2][re]
  // v2.floats[3] = tmp[1][col=2][im]
  
  // These are the top 2 components of the result
  // so I can store them out already.
  _mm_store_ps(dst_shadow, v0.vector);
  _mm_store_ps(dst_shadow+4, v1.vector);
  _mm_store_ps(dst_shadow+8, v2.vector);

  // I want to set up 
  //
  //     v3[0] = v1[3] = tmp[1][col=0][im]
  //     v3[1] = v1[2] = tmp[1][col=0][re]
  //     v3[2] = v2[1] = tmp[1][col=1][im]
  //     v3[3] = v2[0] = tmp[1][col=1][re]
  // 
  // We can do this with v3 = v1,  v3 = shuf( v3, v1 ). Shuf code is 00 01 10 11 = x1B
  v3.vector = v1.vector;
  v3.vector = _mm_shuffle_ps(v3.vector, v2.vector, 0x1B);
  
  // Now I need to do multiply the mask (+1, -1, +1, -1)
  v3.vector = _mm_mul_ps(v3.vector, v6.vector);
  
  // I want to set up 
  //
  //     v4[0] = v2[3] = tmp[1][col=2][im]
  //     v4[1] = v2[2] = tmp[1][col=2][re]
  //     v4[2] = v0[1] = tmp[0][col=0][im]
  //     v4[3] = v0[0] = tmp[0][col=0][re]
  // 
  // We can do this with v4 = v2,  v4 = shuf( v4, v2 ). Shuf code is 00 01 10 11 = x1B
  v4.vector = v2.vector;
  v4.vector = _mm_shuffle_ps(v4.vector, v0.vector, 0x1B);
  
  // Now I need to do multiply the mask (+1, -1, +1, -1)
  v4.vector = _mm_mul_ps(v4.vector, v6.vector);
  
  // I want to set up 
  //
  //     v5[0] = v0[3] = tmp[0][col=1][im]
  //     v5[1] = v0[2] = tmp[0][col=1][re]
  //     v5[2] = v1[1] = tmp[0][col=2][im]
  //     v5[3] = v1[0] = tmp[0][col=2][re]
  // 
  // We can do this with v5 = v0,  v5 = shuf( v5, v1 ). Shuf code is 00 01 10 11 = x1B
  v5.vector = v0.vector;
  v5.vector = _mm_shuffle_ps(v5.vector, v1.vector, 0x1B);
  
  // Now I need to do multiply the mask (+1, -1, +1, -1)
  v5.vector = _mm_mul_ps(v5.vector, v6.vector);
  
  // so I can  them out already - bypass cache.
  _mm_store_ps(dst_shadow+12, v3.vector);
  _mm_store_ps(dst_shadow+16, v4.vector);
  _mm_store_ps(dst_shadow+20, v5.vector);

}

/** \brief Spin recon (1/2)(1-\gamma_0)
 *
 * \param src (pointer to 2 vector source)
 * \param dst (pointer to 4 vector dest)
 * \param n_vec (number of vectors to reconstruct)

 * It is assumeed that src points to an array of  2 spin
 * 3 colour vectors with complex components. 
 * It is assumed that dst points to an array of 4 spin
 * 3 colour vectors with complex components. In both cases
 * ordering is that spin is slowest, and complex components
 * are fastest. */
inline  
void inlineSpinReconDir0Minus(const REAL32* src, REAL32 *dst, unsigned int n_vec) 
{

#ifdef DEBUG_GENERIC_SPIN_RECON_INLINES
  QDPIO::cout << "inlineSpinReconDir0Minus" << endl;
#endif

   /*                              ( 1  0  0 -i)  ( a0 )    ( a0 - i a3 )
    *  B  :=  ( 1 - Gamma  ) A  =  ( 0  1 -i  0)  ( a1 )  = ( a1 - i a2 )
    *                    0         ( 0  i  1  0)  ( a2 )    ( a2 + i a1 )
    *                              ( i  0  0  1)  ( a3 )    ( a3 + i a0 )
    
    * The bottom components of be may be reconstructed using the formula
    *   ( b2r + i b2i )  =  ( {a2r - a1i} + i{a2i + a1r} )  =  ( - b1i + i b1r )
    *   ( b3r + i b3i )     ( {a3r - a0i} + i{a3i + a0r} )     ( - b0i + i b0r ) 
    */


  const REAL32* src_shadow = src;
  REAL32* dst_shadow = dst;

  
  SSEVec v0, v1, v2, v3, v4, v5, v6;
  v6.floats[0] = -1;
  v6.floats[1] = +1;
  v6.floats[2] = -1;
  v6.floats[3] = +1;

  v0.vector = _mm_load_ps(src_shadow);
  v1.vector = _mm_load_ps(src_shadow+4);
  v2.vector = _mm_load_ps(src_shadow+8);
  
  for(unsigned int site=0; site < n_vec-1; site++) { 
    src_shadow += 12;
    _mm_prefetch((const char *)src_shadow, _MM_HINT_T0);

    // Now we have
    // v0.floats[0] = tmp[0][col=0][re]
    // v0.floats[1] = tmp[0][col=0][im]
    // v0.floats[2] = tmp[0][col=1][re]
    // v0.floats[3] = tmp[0][col=1][im]
    //
    // v1.floats[0] = tmp[0][col=2][re]
    // v1.floats[1] = tmp[0][col=2][im]
    // v1.floats[2] = tmp[1][col=0][re]
    // v1.floats[3] = tmp[1][col=0][im]
    //
    // v2.floats[0] = tmp[1][col=1][re]
    // v2.floats[1] = tmp[1][col=1][im]
    // v2.floats[2] = tmp[1][col=2][re]
    // v2.floats[3] = tmp[1][col=2][im]
    
    // These are the top 2 components of the result
    // so I can store them out already.
    _mm_store_ps(dst_shadow, v0.vector);
    _mm_store_ps(dst_shadow+4, v1.vector);
    _mm_store_ps(dst_shadow+8, v2.vector);

    // I want to set up 
    //
    //     v3[0] = v1[3] = tmp[1][col=0][im]
    //     v3[1] = v1[2] = tmp[1][col=0][re]
    //     v3[2] = v2[1] = tmp[1][col=1][im]
    //     v3[3] = v2[0] = tmp[1][col=1][re]
    // 
    // We can do this with v3 = v1,  v3 = shuf( v3, v1 ). Shuf code is 00 01 10 11 = x1B
    v3.vector = v1.vector;
    v3.vector = _mm_shuffle_ps(v3.vector, v2.vector, 0x1B);
    
    // Now I need to do multiply the mask (-1,+1,-1,+1)
    v3.vector = _mm_mul_ps(v3.vector, v6.vector);

    // I want to set up 
    //
    //     v4[0] = v2[3] = tmp[1][col=2][im]
    //     v4[1] = v2[2] = tmp[1][col=2][re]
    //     v4[2] = v0[1] = tmp[0][col=0][im]
    //     v4[3] = v0[0] = tmp[0][col=0][re]
    // 
    // We can do this with v4 = v2,  v4 = shuf( v4, v2 ). Shuf code is 00 01 10 11 = x1B
    v4.vector = v2.vector;
    v4.vector = _mm_shuffle_ps(v4.vector, v0.vector, 0x1B);
    
    // Now I need to do multiply the mask (-1,+1,-1,+1)
    v4.vector = _mm_mul_ps(v4.vector, v6.vector);

    // I want to set up 
    //
    //     v5[0] = v0[3] = tmp[0][col=1][im]
    //     v5[1] = v0[2] = tmp[0][col=1][re]
    //     v5[2] = v1[1] = tmp[0][col=2][im]
    //     v5[3] = v1[0] = tmp[0][col=2][re]
    // 
    // We can do this with v5 = v0,  v5 = shuf( v5, v1 ). Shuf code is 00 01 10 11 = x1B
    v5.vector = v0.vector;
    v5.vector = _mm_shuffle_ps(v5.vector, v1.vector, 0x1B);
    
    // Now I need to do multiply the mask (-1,+1,-1,+1)
    v5.vector = _mm_mul_ps(v5.vector, v6.vector);

    // so I can  them out already - bypass cache.
    _mm_store_ps(dst_shadow+12, v3.vector);
    _mm_store_ps(dst_shadow+16, v4.vector);
    _mm_store_ps(dst_shadow+20, v5.vector);

    dst_shadow+=24;

    v0.vector = _mm_load_ps(src_shadow);
    v1.vector = _mm_load_ps(src_shadow+4);
    v2.vector = _mm_load_ps(src_shadow+8);

  }

  // Now we have
  // v0.floats[0] = tmp[0][col=0][re]
  // v0.floats[1] = tmp[0][col=0][im]
  // v0.floats[2] = tmp[0][col=1][re]
  // v0.floats[3] = tmp[0][col=1][im]
  //
  // v1.floats[0] = tmp[0][col=2][re]
  // v1.floats[1] = tmp[0][col=2][im]
  // v1.floats[2] = tmp[1][col=0][re]
  // v1.floats[3] = tmp[1][col=0][im]
  //
  // v2.floats[0] = tmp[1][col=1][re]
  // v2.floats[1] = tmp[1][col=1][im]
  // v2.floats[2] = tmp[1][col=2][re]
  // v2.floats[3] = tmp[1][col=2][im]
  
  // These are the top 2 components of the result
  // so I can store them out already.
  _mm_store_ps(dst_shadow, v0.vector);
  _mm_store_ps(dst_shadow+4, v1.vector);
  _mm_store_ps(dst_shadow+8, v2.vector);

  // I want to set up 
  //
  //     v3[0] = v1[3] = tmp[1][col=0][im]
  //     v3[1] = v1[2] = tmp[1][col=0][re]
  //     v3[2] = v2[1] = tmp[1][col=1][im]
  //     v3[3] = v2[0] = tmp[1][col=1][re]
  // 
  // We can do this with v3 = v1,  v3 = shuf( v3, v1 ). Shuf code is 00 01 10 11 = x1B
  v3.vector = v1.vector;
  v3.vector = _mm_shuffle_ps(v3.vector, v2.vector, 0x1B);
  
  // Now I need to do multiply the mask (-1,+1,-1,+1)
  v3.vector = _mm_mul_ps(v3.vector, v6.vector);
  
  // I want to set up 
  //
  //     v4[0] = v2[3] = tmp[1][col=2][im]
  //     v4[1] = v2[2] = tmp[1][col=2][re]
  //     v4[2] = v0[1] = tmp[0][col=0][im]
  //     v4[3] = v0[0] = tmp[0][col=0][re]
  // 
  // We can do this with v4 = v2,  v4 = shuf( v4, v2 ). Shuf code is 00 01 10 11 = x1B
  v4.vector = v2.vector;
  v4.vector = _mm_shuffle_ps(v4.vector, v0.vector, 0x1B);
  
  // Now I need to do multiply the mask (-1,+1,-1,+1)
  v4.vector = _mm_mul_ps(v4.vector, v6.vector);
  
  // I want to set up 
  //
  //     v5[0] = v0[3] = tmp[0][col=1][im]
  //     v5[1] = v0[2] = tmp[0][col=1][re]
  //     v5[2] = v1[1] = tmp[0][col=2][im]
  //     v5[3] = v1[0] = tmp[0][col=2][re]
  // 
  // We can do this with v5 = v0,  v5 = shuf( v5, v1 ). Shuf code is 00 01 10 11 = x1B
  v5.vector = v0.vector;
  v5.vector = _mm_shuffle_ps(v5.vector, v1.vector, 0x1B);
  
  // Now I need to do multiply the mask (-1,+1,-1,+1)
  v5.vector = _mm_mul_ps(v5.vector, v6.vector);
  
  // so I can  them out already - bypass cache.
  _mm_store_ps(dst_shadow+12, v3.vector);
  _mm_store_ps(dst_shadow+16, v4.vector);
  _mm_store_ps(dst_shadow+20, v5.vector);

}

/** \brief Spin recon (1/2)(1+\gamma_1)
 *
 * \param src (pointer to 2 vector source)
 * \param dst (pointer to 4 vector dest)
 * \param n_vec (number of vectors to reconstruct)

 * It is assumeed that src points to an array of  2 spin
 * 3 colour vectors with complex components. 
 * It is assumed that dst points to an array of 4 spin
 * 3 colour vectors with complex components. In both cases
 * ordering is that spin is slowest, and complex components
 * are fastest. */
inline  
void inlineSpinReconDir1Plus(const REAL32* src, REAL32 *dst, unsigned int n_vec) 
{

#ifdef DEBUG_GENERIC_SPIN_RECON_INLINES
  QDPIO::cout << "inlineSpinReconDir0Plus" << endl;
#endif

    /* 1 + \gamma_1 =  1  0  0 -1 
     *                 0  1  1  0
     *                 0  1  1  0
     *                -1  0  0  1 
 
     *   ( b2r + i b2i )  =  ( {a2r + a1r} + i{a2i + a1i} )  =  (   b1r + i b1i )
     *   ( b3r + i b3i )     ( {a3r - a0r} + i{a3i - a0i} )     ( - b0r - i b0i ) 
  
    */

  const REAL32* src_shadow = src;
  REAL32* dst_shadow = dst;

  
  SSEVec v0, v1, v2, v3, v4, v5, v6, v7;
  v6.floats[0] = +1;
  v6.floats[1] = +1;
  v6.floats[2] = -1;
  v6.floats[3] = -1;


  v7.floats[0] = -1;
  v7.floats[1] = -1;
  v7.floats[2] = -1;
  v7.floats[3] = -1;


  v0.vector = _mm_load_ps(src_shadow);
  v1.vector = _mm_load_ps(src_shadow+4);
  v2.vector = _mm_load_ps(src_shadow+8);
  
  for(unsigned int site=0; site < n_vec-1; site++) { 

    src_shadow += 12;
    _mm_prefetch((const char *)src_shadow, _MM_HINT_T0);

    // Now we have
    // v0.floats[0] = tmp[0][col=0][re]
    // v0.floats[1] = tmp[0][col=0][im]
    // v0.floats[2] = tmp[0][col=1][re]
    // v0.floats[3] = tmp[0][col=1][im]
    //
    // v1.floats[0] = tmp[0][col=2][re]
    // v1.floats[1] = tmp[0][col=2][im]
    // v1.floats[2] = tmp[1][col=0][re]
    // v1.floats[3] = tmp[1][col=0][im]
    //
    // v2.floats[0] = tmp[1][col=1][re]
    // v2.floats[1] = tmp[1][col=1][im]
    // v2.floats[2] = tmp[1][col=2][re]
    // v2.floats[3] = tmp[1][col=2][im]
    
    // These are the top 2 components of the result
    // so I can store them out already.
    _mm_store_ps(dst_shadow, v0.vector);
    _mm_store_ps(dst_shadow+4, v1.vector);
    _mm_store_ps(dst_shadow+8, v2.vector);

    // I want to set up 
    //
    //     v3[0] = v1[2] = tmp[1][col=0][re]
    //     v3[1] = v1[3] = tmp[1][col=0][im]
    //     v3[2] = v2[0] = tmp[1][col=1][re]
    //     v3[3] = v2[1] = tmp[1][col=1][im]
    // 
    // We can do this with v3 = v1,  v3 = shuf( v3, v2 ). Shuf code is 01 00 11 10 = x4E
    v3.vector = v1.vector;
    v3.vector = _mm_shuffle_ps(v3.vector, v2.vector, 0x4E);
    
    // I want to set up 
    //
    //     v4[0] = v2[2] = tmp[1][col=2][re]
    //     v4[1] = v2[3] = tmp[1][col=2][im]
    //     v4[2] = v0[0] = tmp[0][col=0][re]
    //     v4[3] = v0[1] = tmp[0][col=0][im]
    // 
    // We can do this with v4 = v2,  v4 = shuf( v4, v0 ). Shuf code is 01 00 11 10 = x4E
    v4.vector = v2.vector;
    v4.vector = _mm_shuffle_ps(v4.vector, v0.vector, 0x4E);

    // Need to multiply in mask (+1,+1,-1,-1) 
    v4.vector = _mm_mul_ps(v4.vector, v6.vector);

    // I want to set up 
    //
    //     v5[0] = v0[2] = tmp[0][col=1][re]
    //     v5[1] = v0[3] = tmp[0][col=1][im]
    //     v5[2] = v1[0] = tmp[0][col=2][re]
    //     v5[3] = v1[1] = tmp[0][col=2][im]
    // 
    // We can do this with v5 = v0,  v5 = shuf( v5, v1 ). Shuf code is 01 00 11 10 = x4E
    v5.vector = v0.vector;
    v5.vector = _mm_shuffle_ps(v5.vector, v1.vector, 0x4E);
    
    // Now I need to do multiply the mask (-1,-1,-1,-11)
    v5.vector = _mm_mul_ps(v5.vector, v7.vector);

    // so I can  them out already - bypass cache.
    _mm_store_ps(dst_shadow+12, v3.vector);
    _mm_store_ps(dst_shadow+16, v4.vector);
    _mm_store_ps(dst_shadow+20, v5.vector);

    dst_shadow+=24;

    v0.vector = _mm_load_ps(src_shadow);
    v1.vector = _mm_load_ps(src_shadow+4);
    v2.vector = _mm_load_ps(src_shadow+8);

  }

  // Now we have
  // v0.floats[0] = tmp[0][col=0][re]
  // v0.floats[1] = tmp[0][col=0][im]
  // v0.floats[2] = tmp[0][col=1][re]
  // v0.floats[3] = tmp[0][col=1][im]
  //
  // v1.floats[0] = tmp[0][col=2][re]
  // v1.floats[1] = tmp[0][col=2][im]
  // v1.floats[2] = tmp[1][col=0][re]
  // v1.floats[3] = tmp[1][col=0][im]
  //
  // v2.floats[0] = tmp[1][col=1][re]
  // v2.floats[1] = tmp[1][col=1][im]
  // v2.floats[2] = tmp[1][col=2][re]
  // v2.floats[3] = tmp[1][col=2][im]
  
  // These are the top 2 components of the result
  // so I can store them out already.
  _mm_store_ps(dst_shadow, v0.vector);
  _mm_store_ps(dst_shadow+4, v1.vector);
  _mm_store_ps(dst_shadow+8, v2.vector);
  
  // I want to set up 
  //
  //     v3[0] = v1[2] = tmp[1][col=0][re]
  //     v3[1] = v1[3] = tmp[1][col=0][im]
  //     v3[2] = v2[0] = tmp[1][col=1][re]
  //     v3[3] = v2[1] = tmp[1][col=1][im]
  // 
  // We can do this with v3 = v1,  v3 = shuf( v3, v2 ). Shuf code is 01 00 11 10 = x4E
  v3.vector = v1.vector;
  v3.vector = _mm_shuffle_ps(v3.vector, v2.vector, 0x4E);
  
  // I want to set up 
  //
  //     v4[0] = v2[2] = tmp[1][col=2][re]
  //     v4[1] = v2[3] = tmp[1][col=2][im]
  //     v4[2] = v0[0] = tmp[0][col=0][re]
  //     v4[3] = v0[1] = tmp[0][col=0][im]
  // 
  // We can do this with v4 = v2,  v4 = shuf( v4, v0 ). Shuf code is 01 00 11 10 = x4E
  v4.vector = v2.vector;
  v4.vector = _mm_shuffle_ps(v4.vector, v0.vector, 0x4E);
  
  // Need to multiply in mask (+1,+1,-1,-1) 
  v4.vector = _mm_mul_ps(v4.vector, v6.vector);
  
  // I want to set up 
  //
  //     v5[0] = v0[2] = tmp[0][col=1][re]
  //     v5[1] = v0[3] = tmp[0][col=1][im]
  //     v5[2] = v1[0] = tmp[0][col=2][re]
  //     v5[3] = v1[1] = tmp[0][col=2][im]
  // 
  // We can do this with v5 = v0,  v5 = shuf( v5, v1 ). Shuf code is 01 00 11 10 = x4E
  v5.vector = v0.vector;
  v5.vector = _mm_shuffle_ps(v5.vector, v1.vector, 0x4E);
  
  // Now I need to do multiply the mask (-1,-1,-1,-11)
  v5.vector = _mm_mul_ps(v5.vector, v7.vector);
  
  // so I can  them out already - bypass cache.
  _mm_store_ps(dst_shadow+12, v3.vector);
  _mm_store_ps(dst_shadow+16, v4.vector);
  _mm_store_ps(dst_shadow+20, v5.vector);
  
}

/** \brief Spin recon (1/2)(1-\gamma_1)
 *
 * \param src (pointer to 2 vector source)
 * \param dst (pointer to 4 vector dest)
 * \param n_vec (number of vectors to reconstruct)

 * It is assumeed that src points to an array of  2 spin
 * 3 colour vectors with complex components. 
 * It is assumed that dst points to an array of 4 spin
 * 3 colour vectors with complex components. In both cases
 * ordering is that spin is slowest, and complex components
 * are fastest. */
inline  
void inlineSpinReconDir1Minus(const REAL32* src, REAL32 *dst, unsigned int n_vec) 
{

#ifdef DEBUG_GENERIC_SPIN_RECON_INLINES
  QDPIO::cout << "inlineSpinReconDir0Minus" << endl;
#endif

  const REAL32* src_shadow = src;
  REAL32* dst_shadow = dst;

  
  SSEVec v0, v1, v2, v3, v4, v5, v6, v7;
  v6.floats[0] = -1;
  v6.floats[1] = -1;
  v6.floats[2] = +1;
  v6.floats[3] = +1;


  v7.floats[0] = -1;
  v7.floats[1] = -1;
  v7.floats[2] = -1;
  v7.floats[3] = -1;

  v0.vector = _mm_load_ps(src_shadow);
  v1.vector = _mm_load_ps(src_shadow+4);
  v2.vector = _mm_load_ps(src_shadow+8);

  for(unsigned int site=0; site < n_vec-1; site++) { 

    src_shadow += 12;
    _mm_prefetch((const char *)src_shadow, _MM_HINT_T0);

    // Now we have
    // v0.floats[0] = tmp[0][col=0][re]
    // v0.floats[1] = tmp[0][col=0][im]
    // v0.floats[2] = tmp[0][col=1][re]
    // v0.floats[3] = tmp[0][col=1][im]
    //
    // v1.floats[0] = tmp[0][col=2][re]
    // v1.floats[1] = tmp[0][col=2][im]
    // v1.floats[2] = tmp[1][col=0][re]
    // v1.floats[3] = tmp[1][col=0][im]
    //
    // v2.floats[0] = tmp[1][col=1][re]
    // v2.floats[1] = tmp[1][col=1][im]
    // v2.floats[2] = tmp[1][col=2][re]
    // v2.floats[3] = tmp[1][col=2][im]
    
    // These are the top 2 components of the result
    // so I can store them out already.
    _mm_store_ps(dst_shadow, v0.vector);
    _mm_store_ps(dst_shadow+4, v1.vector);
    _mm_store_ps(dst_shadow+8, v2.vector);

    // I want to set up 
    //
    //     v3[0] = v1[2] = tmp[1][col=0][re]
    //     v3[1] = v1[3] = tmp[1][col=0][im]
    //     v3[2] = v2[0] = tmp[1][col=1][re]
    //     v3[3] = v2[1] = tmp[1][col=1][im]
    // 
    // We can do this with v3 = v1,  v3 = shuf( v3, v2 ). Shuf code is 01 00 11 10 = x4E
    v3.vector = v1.vector;
    v3.vector = _mm_shuffle_ps(v3.vector, v2.vector, 0x4E);

    // Mask (-1,-1,-1,-1)
    v3.vector = _mm_mul_ps(v3.vector, v7.vector);

    // I want to set up 
    //
    //     v4[0] = v2[2] = tmp[1][col=2][re]
    //     v4[1] = v2[3] = tmp[1][col=2][im]
    //     v4[2] = v0[0] = tmp[0][col=0][re]
    //     v4[3] = v0[1] = tmp[0][col=0][im]
    // 
    // We can do this with v4 = v2,  v4 = shuf( v4, v0 ). Shuf code is 01 00 11 10 = x4E
    v4.vector = v2.vector;
    v4.vector = _mm_shuffle_ps(v4.vector, v0.vector, 0x4E);

    // Need to multiply in mask (-1,-1,+1,+1) 
    v4.vector = _mm_mul_ps(v4.vector, v6.vector);

    // I want to set up 
    //
    //     v5[0] = v0[2] = tmp[0][col=1][re]
    //     v5[1] = v0[3] = tmp[0][col=1][im]
    //     v5[2] = v1[0] = tmp[0][col=2][re]
    //     v5[3] = v1[1] = tmp[0][col=2][im]
    // 
    // We can do this with v5 = v0,  v5 = shuf( v5, v1 ). Shuf code is 01 00 11 10 = x4E
    v5.vector = v0.vector;
    v5.vector = _mm_shuffle_ps(v5.vector, v1.vector, 0x4E);
    
    // so I can  them out already - bypass cache.
    _mm_store_ps(dst_shadow+12, v3.vector);
    _mm_store_ps(dst_shadow+16, v4.vector);
    _mm_store_ps(dst_shadow+20, v5.vector);

    dst_shadow+=24;

    v0.vector = _mm_load_ps(src_shadow);
    v1.vector = _mm_load_ps(src_shadow+4);
    v2.vector = _mm_load_ps(src_shadow+8);


  }

  // Now we have
  // v0.floats[0] = tmp[0][col=0][re]
  // v0.floats[1] = tmp[0][col=0][im]
  // v0.floats[2] = tmp[0][col=1][re]
  // v0.floats[3] = tmp[0][col=1][im]
  //
  // v1.floats[0] = tmp[0][col=2][re]
  // v1.floats[1] = tmp[0][col=2][im]
  // v1.floats[2] = tmp[1][col=0][re]
  // v1.floats[3] = tmp[1][col=0][im]
  //
  // v2.floats[0] = tmp[1][col=1][re]
  // v2.floats[1] = tmp[1][col=1][im]
  // v2.floats[2] = tmp[1][col=2][re]
  // v2.floats[3] = tmp[1][col=2][im]
  
  // These are the top 2 components of the result
  // so I can store them out already.
  _mm_store_ps(dst_shadow, v0.vector);
  _mm_store_ps(dst_shadow+4, v1.vector);
  _mm_store_ps(dst_shadow+8, v2.vector);
  
  // I want to set up 
  //
  //     v3[0] = v1[2] = tmp[1][col=0][re]
  //     v3[1] = v1[3] = tmp[1][col=0][im]
  //     v3[2] = v2[0] = tmp[1][col=1][re]
  //     v3[3] = v2[1] = tmp[1][col=1][im]
  // 
  // We can do this with v3 = v1,  v3 = shuf( v3, v2 ). Shuf code is 01 00 11 10 = x4E
  v3.vector = v1.vector;
  v3.vector = _mm_shuffle_ps(v3.vector, v2.vector, 0x4E);
  
  // Mask (-1,-1,-1,-1)
  v3.vector = _mm_mul_ps(v3.vector, v7.vector);
  
  // I want to set up 
  //
  //     v4[0] = v2[2] = tmp[1][col=2][re]
  //     v4[1] = v2[3] = tmp[1][col=2][im]
  //     v4[2] = v0[0] = tmp[0][col=0][re]
  //     v4[3] = v0[1] = tmp[0][col=0][im]
  // 
  // We can do this with v4 = v2,  v4 = shuf( v4, v0 ). Shuf code is 01 00 11 10 = x4E
  v4.vector = v2.vector;
  v4.vector = _mm_shuffle_ps(v4.vector, v0.vector, 0x4E);
  
  // Need to multiply in mask (-1,-1,+1,+1) 
  v4.vector = _mm_mul_ps(v4.vector, v6.vector);
  
  // I want to set up 
  //
  //     v5[0] = v0[2] = tmp[0][col=1][re]
  //     v5[1] = v0[3] = tmp[0][col=1][im]
  //     v5[2] = v1[0] = tmp[0][col=2][re]
  //     v5[3] = v1[1] = tmp[0][col=2][im]
  // 
  // We can do this with v5 = v0,  v5 = shuf( v5, v1 ). Shuf code is 01 00 11 10 = x4E
  v5.vector = v0.vector;
  v5.vector = _mm_shuffle_ps(v5.vector, v1.vector, 0x4E);
  
  // so I can  them out already - bypass cache.
  _mm_store_ps(dst_shadow+12, v3.vector);
  _mm_store_ps(dst_shadow+16, v4.vector);
  _mm_store_ps(dst_shadow+20, v5.vector);

}


/** \brief Spin recon (1/2)(1+\gamma_2)
 *
 * \param src (pointer to 2 vector source)
 * \param dst (pointer to 4 vector dest)
 * \param n_vec (number of vectors to reconstruct)

 * It is assumeed that src points to an array of  2 spin
 * 3 colour vectors with complex components. 
 * It is assumed that dst points to an array of 4 spin
 * 3 colour vectors with complex components. In both cases
 * ordering is that spin is slowest, and complex components
 * are fastest. */
inline  
void inlineSpinReconDir2Plus(const REAL32* src, REAL32 *dst, unsigned int n_vec) 
{

#ifdef DEBUG_GENERIC_SPIN_RECON_INLINES
  QDPIO::cout << "inlineSpinReconDir0Plus" << endl;
#endif


  /* 1 + \gamma_2 =  1  0  i  0 
   *                 0  1  0 -i
   *                -i  0  1  0
   *                 0  i  0  1 
   *	     
   *  ( b2r + i b2i )  =  ( {a2r + a0i} + i{a2i - a0r} )  =  (   b0i - i b0r )
   *  ( b3r + i b3i )     ( {a3r - a1i} + i{a3i + a1r} )     ( - b1i + i b1r ) 
  */    
  
  SSEVec v0, v1, v2, v5, v6, v7;
  const REAL32* src_shadow = src;
  REAL32* dst_shadow = dst;

  v5.floats[0] = +1;
  v5.floats[1] = -1;
  v5.floats[2] = +1;
  v5.floats[3] = -1;

  v6.floats[0] = +1;
  v6.floats[1] = -1;
  v6.floats[2] = -1;
  v6.floats[3] = +1;

  v7.floats[0] = -1;
  v7.floats[1] = +1;
  v7.floats[2] = -1;
  v7.floats[3] = +1;

  v0.vector = _mm_load_ps(src_shadow);
  v1.vector = _mm_load_ps(src_shadow+4);
  v2.vector = _mm_load_ps(src_shadow+8);

  for(unsigned int site=0; site < n_vec-1; site++) { 

    src_shadow += 12;
    _mm_prefetch((const char *)src_shadow, _MM_HINT_T0);

    // Now we have
    // v0.floats[0] = tmp[0][col=0][re]
    // v0.floats[1] = tmp[0][col=0][im]
    // v0.floats[2] = tmp[0][col=1][re]
    // v0.floats[3] = tmp[0][col=1][im]
    //
    // v1.floats[0] = tmp[0][col=2][re]
    // v1.floats[1] = tmp[0][col=2][im]
    // v1.floats[2] = tmp[1][col=0][re]
    // v1.floats[3] = tmp[1][col=0][im]
    //
    // v2.floats[0] = tmp[1][col=1][re]
    // v2.floats[1] = tmp[1][col=1][im]
    // v2.floats[2] = tmp[1][col=2][re]
    // v2.floats[3] = tmp[1][col=2][im]
    
    // These are the top 2 components of the result
    // so I can store them out already.
    _mm_store_ps(dst_shadow, v0.vector);
    _mm_store_ps(dst_shadow+4, v1.vector);
    _mm_store_ps(dst_shadow+8, v2.vector);

    // I want to set up 
    //
    //     v0[0] = v0[1] = tmp[0][col=0][im]
    //     v0[1] = v0[0] = tmp[0][col=0][re]
    //     v0[2] = v0[3] = tmp[0][col=1][im]
    //     v0[3] = v0[2] = tmp[0][col=1][re]
    // 
    // We can do this with v0 = shuf( v0, v0 ). Shuf code is 10 11 00 01 = xB1
    v0.vector = _mm_shuffle_ps(v0.vector, v0.vector, 0xB1);

    // Mask (+1,-1,+1,-1)
    v0.vector = _mm_mul_ps(v0.vector, v5.vector);

    // I want to set up 
    //
    //     v1[0] = v1[1] = tmp[0][col=2][im]
    //     v1[1] = v1[0] = tmp[0][col=2][re]
    //     v1[2] = v1[3] = tmp[1][col=0][im]
    //     v1[3] = v1[2] = tmp[1][col=0][re]
    // 
    // We can do this with v1 = shuf( v1, v1 ). Shuf code is 10 11 00 01 = xB1
    v1.vector = _mm_shuffle_ps(v1.vector, v1.vector, 0xB1);

    // Need to multiply in mask (+1,-1,-1,+1) 
    v1.vector = _mm_mul_ps(v1.vector, v6.vector);

    // I want to set up 
    //
    //     v2[0] = v2[1] = tmp[1][col=1][im]
    //     v2[1] = v2[0] = tmp[1][col=1][re]
    //     v2[2] = v2[2] = tmp[1][col=2][im]
    //     v2[3] = v2[3] = tmp[1][col=2][re]
    // 
    // We can do this with v2 = shuf( v2, v2 ). Shuf code is 10 11 00 01 = xB1
    v2.vector = _mm_shuffle_ps(v2.vector, v2.vector, 0xB1);

    // Mask in -1, +1, -1, +1 from v7
    v2.vector = _mm_mul_ps(v2.vector, v7.vector);

    // so I can  them out already - bypass cache.
    _mm_store_ps(dst_shadow+12, v0.vector);
    _mm_store_ps(dst_shadow+16, v1.vector);
    _mm_store_ps(dst_shadow+20, v2.vector);

    dst_shadow+=24;

    v0.vector = _mm_load_ps(src_shadow);
    v1.vector = _mm_load_ps(src_shadow+4);
    v2.vector = _mm_load_ps(src_shadow+8);

    
  }
  
  // Now we have
  // v0.floats[0] = tmp[0][col=0][re]
  // v0.floats[1] = tmp[0][col=0][im]
  // v0.floats[2] = tmp[0][col=1][re]
  // v0.floats[3] = tmp[0][col=1][im]
  //
  // v1.floats[0] = tmp[0][col=2][re]
  // v1.floats[1] = tmp[0][col=2][im]
  // v1.floats[2] = tmp[1][col=0][re]
  // v1.floats[3] = tmp[1][col=0][im]
  //
  // v2.floats[0] = tmp[1][col=1][re]
  // v2.floats[1] = tmp[1][col=1][im]
  // v2.floats[2] = tmp[1][col=2][re]
  // v2.floats[3] = tmp[1][col=2][im]
  
  // These are the top 2 components of the result
  // so I can store them out already.
  _mm_store_ps(dst_shadow, v0.vector);
  _mm_store_ps(dst_shadow+4, v1.vector);
  _mm_store_ps(dst_shadow+8, v2.vector);
  
  // I want to set up 
  //
  //     v0[0] = v0[1] = tmp[0][col=0][im]
  //     v0[1] = v0[0] = tmp[0][col=0][re]
  //     v0[2] = v0[3] = tmp[0][col=1][im]
  //     v0[3] = v0[2] = tmp[0][col=1][re]
  // 
  // We can do this with v0 = shuf( v0, v0 ). Shuf code is 10 11 00 01 = xB1
  v0.vector = _mm_shuffle_ps(v0.vector, v0.vector, 0xB1);
  
  // Mask (+1,-1,+1,-1)
  v0.vector = _mm_mul_ps(v0.vector, v5.vector);
  
  // I want to set up 
  //
  //     v1[0] = v1[1] = tmp[0][col=2][im]
  //     v1[1] = v1[0] = tmp[0][col=2][re]
  //     v1[2] = v1[3] = tmp[1][col=0][im]
  //     v1[3] = v1[2] = tmp[1][col=0][re]
  // 
  // We can do this with v1 = shuf( v1, v1 ). Shuf code is 10 11 00 01 = xB1
  v1.vector = _mm_shuffle_ps(v1.vector, v1.vector, 0xB1);
  
  // Need to multiply in mask (+1,-1,-1,+1) 
  v1.vector = _mm_mul_ps(v1.vector, v6.vector);
  
  // I want to set up 
  //
  //     v2[0] = v2[1] = tmp[1][col=1][im]
  //     v2[1] = v2[0] = tmp[1][col=1][re]
  //     v2[2] = v2[2] = tmp[1][col=2][im]
  //     v2[3] = v2[3] = tmp[1][col=2][re]
  // 
  // We can do this with v2 = shuf( v2, v2 ). Shuf code is 10 11 00 01 = xB1
  v2.vector = _mm_shuffle_ps(v2.vector, v2.vector, 0xB1);
  
  // Mask in -1, +1, -1, +1 from v7
  v2.vector = _mm_mul_ps(v2.vector, v7.vector);
  
  // so I can  them out already - bypass cache.
  _mm_store_ps(dst_shadow+12, v0.vector);
  _mm_store_ps(dst_shadow+16, v1.vector);
  _mm_store_ps(dst_shadow+20, v2.vector);


}

/** \brief Spin recon (1/2)(1-\gamma_2)
 *
 * \param src (pointer to 2 vector source)
 * \param dst (pointer to 4 vector dest)
 * \param n_vec (number of vectors to reconstruct)

 * It is assumeed that src points to an array of  2 spin
 * 3 colour vectors with complex components. 
 * It is assumed that dst points to an array of 4 spin
 * 3 colour vectors with complex components. In both cases
 * ordering is that spin is slowest, and complex components
 * are fastest. */
inline  
void inlineSpinReconDir2Minus(const REAL32* src, REAL32 *dst, unsigned int n_vec) 
{

#ifdef DEBUG_GENERIC_SPIN_RECON_INLINES
  QDPIO::cout << "inlineSpinReconDir2Minus" << endl;
#endif

  SSEVec v0, v1, v2, v5, v6, v7;
  const REAL32* src_shadow = src;
  REAL32* dst_shadow = dst;

  v5.floats[0] = -1;
  v5.floats[1] = +1;
  v5.floats[2] = -1;
  v5.floats[3] = +1;

  v6.floats[0] = -1;
  v6.floats[1] = +1;
  v6.floats[2] = +1;
  v6.floats[3] = -1;

  v7.floats[0] = +1;
  v7.floats[1] = -1;
  v7.floats[2] = +1;
  v7.floats[3] = -1;

  v0.vector = _mm_load_ps(src_shadow);
  v1.vector = _mm_load_ps(src_shadow+4);
  v2.vector = _mm_load_ps(src_shadow+8);

  for(unsigned int site=0; site < n_vec-1; site++) { 

    src_shadow += 12;
    _mm_prefetch((const char *)src_shadow, _MM_HINT_T0);

    // Now we have
    // v0.floats[0] = tmp[0][col=0][re]
    // v0.floats[1] = tmp[0][col=0][im]
    // v0.floats[2] = tmp[0][col=1][re]
    // v0.floats[3] = tmp[0][col=1][im]
    //
    // v1.floats[0] = tmp[0][col=2][re]
    // v1.floats[1] = tmp[0][col=2][im]
    // v1.floats[2] = tmp[1][col=0][re]
    // v1.floats[3] = tmp[1][col=0][im]
    //
    // v2.floats[0] = tmp[1][col=1][re]
    // v2.floats[1] = tmp[1][col=1][im]
    // v2.floats[2] = tmp[1][col=2][re]
    // v2.floats[3] = tmp[1][col=2][im]
    
    // These are the top 2 components of the result
    // so I can store them out already.
    _mm_store_ps(dst_shadow, v0.vector);
    _mm_store_ps(dst_shadow+4, v1.vector);
    _mm_store_ps(dst_shadow+8, v2.vector);

    // I want to set up 
    //
    //     v0[0] = v0[1] = tmp[0][col=0][im]
    //     v0[1] = v0[0] = tmp[0][col=0][re]
    //     v0[2] = v0[3] = tmp[0][col=1][im]
    //     v0[3] = v0[2] = tmp[0][col=1][re]
    // 
    // We can do this with v0 = shuf( v0, v0 ). Shuf code is 10 11 00 01 = xB1
    v0.vector = _mm_shuffle_ps(v0.vector, v0.vector, 0xB1);

    // Mask (+1,-1,+1,-1)
    v0.vector = _mm_mul_ps(v0.vector, v5.vector);

    // I want to set up 
    //
    //     v1[0] = v1[1] = tmp[0][col=2][im]
    //     v1[1] = v1[0] = tmp[0][col=2][re]
    //     v1[2] = v1[3] = tmp[1][col=0][im]
    //     v1[3] = v1[2] = tmp[1][col=0][re]
    // 
    // We can do this with v1 = shuf( v1, v1 ). Shuf code is 10 11 00 01 = xB1
    v1.vector = _mm_shuffle_ps(v1.vector, v1.vector, 0xB1);

    // Need to multiply in mask (+1,-1,-1,+1) 
    v1.vector = _mm_mul_ps(v1.vector, v6.vector);

    // I want to set up 
    //
    //     v2[0] = v2[1] = tmp[1][col=1][im]
    //     v2[1] = v2[0] = tmp[1][col=1][re]
    //     v2[2] = v2[2] = tmp[1][col=2][im]
    //     v2[3] = v2[3] = tmp[1][col=2][re]
    // 
    // We can do this with v2 = shuf( v2, v2 ). Shuf code is 10 11 00 01 = xB1
    v2.vector = _mm_shuffle_ps(v2.vector, v2.vector, 0xB1);

    // Mask in -1, +1, -1, +1 from v7
    v2.vector = _mm_mul_ps(v2.vector, v7.vector);

    // so I can  them out already - bypass cache.
    _mm_store_ps(dst_shadow+12, v0.vector);
    _mm_store_ps(dst_shadow+16, v1.vector);
    _mm_store_ps(dst_shadow+20, v2.vector);

    dst_shadow+=24;

    v0.vector = _mm_load_ps(src_shadow);
    v1.vector = _mm_load_ps(src_shadow+4);
    v2.vector = _mm_load_ps(src_shadow+8);

    
  }
  
  // Now we have
  // v0.floats[0] = tmp[0][col=0][re]
  // v0.floats[1] = tmp[0][col=0][im]
  // v0.floats[2] = tmp[0][col=1][re]
  // v0.floats[3] = tmp[0][col=1][im]
  //
  // v1.floats[0] = tmp[0][col=2][re]
  // v1.floats[1] = tmp[0][col=2][im]
  // v1.floats[2] = tmp[1][col=0][re]
  // v1.floats[3] = tmp[1][col=0][im]
  //
  // v2.floats[0] = tmp[1][col=1][re]
  // v2.floats[1] = tmp[1][col=1][im]
  // v2.floats[2] = tmp[1][col=2][re]
  // v2.floats[3] = tmp[1][col=2][im]
  
  // These are the top 2 components of the result
  // so I can store them out already.
  _mm_store_ps(dst_shadow, v0.vector);
  _mm_store_ps(dst_shadow+4, v1.vector);
  _mm_store_ps(dst_shadow+8, v2.vector);
  
  // I want to set up 
  //
  //     v0[0] = v0[1] = tmp[0][col=0][im]
  //     v0[1] = v0[0] = tmp[0][col=0][re]
  //     v0[2] = v0[3] = tmp[0][col=1][im]
  //     v0[3] = v0[2] = tmp[0][col=1][re]
  // 
  // We can do this with v0 = shuf( v0, v0 ). Shuf code is 10 11 00 01 = xB1
  v0.vector = _mm_shuffle_ps(v0.vector, v0.vector, 0xB1);
  
  // Mask (+1,-1,+1,-1)
  v0.vector = _mm_mul_ps(v0.vector, v5.vector);
  
  // I want to set up 
  //
  //     v1[0] = v1[1] = tmp[0][col=2][im]
  //     v1[1] = v1[0] = tmp[0][col=2][re]
  //     v1[2] = v1[3] = tmp[1][col=0][im]
  //     v1[3] = v1[2] = tmp[1][col=0][re]
  // 
  // We can do this with v1 = shuf( v1, v1 ). Shuf code is 10 11 00 01 = xB1
  v1.vector = _mm_shuffle_ps(v1.vector, v1.vector, 0xB1);
  
  // Need to multiply in mask (+1,-1,-1,+1) 
  v1.vector = _mm_mul_ps(v1.vector, v6.vector);
  
  // I want to set up 
  //
  //     v2[0] = v2[1] = tmp[1][col=1][im]
  //     v2[1] = v2[0] = tmp[1][col=1][re]
  //     v2[2] = v2[2] = tmp[1][col=2][im]
  //     v2[3] = v2[3] = tmp[1][col=2][re]
  // 
  // We can do this with v2 = shuf( v2, v2 ). Shuf code is 10 11 00 01 = xB1
  v2.vector = _mm_shuffle_ps(v2.vector, v2.vector, 0xB1);
  
  // Mask in -1, +1, -1, +1 from v7
  v2.vector = _mm_mul_ps(v2.vector, v7.vector);
  
  // so I can  them out already - bypass cache.
  _mm_store_ps(dst_shadow+12, v0.vector);
  _mm_store_ps(dst_shadow+16, v1.vector);
  _mm_store_ps(dst_shadow+20, v2.vector);

}


/** \brief Spin recon (1/2)(1+\gamma3)
 *
 * \param src (pointer to 2 vector source)
 * \param dst (pointer to 4 vector dest)
 * \param n_vec (number of vectors to reconstruct)

 * It is assumeed that src points to an array of  2 spin
 * 3 colour vectors with complex components. 
 * It is assumed that dst points to an array of 4 spin
 * 3 colour vectors with complex components. In both cases
 * ordering is that spin is slowest, and complex components
 * are fastest. */
inline  
void inlineSpinReconDir3Plus(const REAL32* src, REAL32 *dst, unsigned int n_vec) 
{

#ifdef DEBUG_GENERIC_SPIN_RECON_INLINES
  QDPIO::cout << "inlineSpinReconDir0Plus" << endl;
#endif

  SSEVec v0, v1, v2;
  const REAL32* src_shadow = src;
  REAL32* dst_shadow = dst;

  /*                              ( 1  0  1  0)  ( a0 )    ( a0 + a2 )
   *  B  :=  ( 1 + Gamma  ) A  =  ( 0  1  0  1)  ( a1 )  = ( a1 + a3 )
   *                    3         ( 1  0  1  0)  ( a2 )    ( a2 + a0 )
   *                              ( 0  1  0  1)  ( a3 )    ( a3 + a1 )
   
   * The bottom components of be may be reconstructed using the formula
   
   *   ( b2r + i b2i )  =  ( {a2r + a0r} + i{a2i + a0i} )  =  ( b0r + i b0i )
   *   ( b3r + i b3i )     ( {a3r + a1r} + i{a3i + a1i} )     ( b1r + i b1i ) 
   */

  v0.vector = _mm_load_ps(src_shadow);
  v1.vector = _mm_load_ps(src_shadow+4);
  v2.vector = _mm_load_ps(src_shadow+8);

  for(unsigned int site=0; site < n_vec-1; site++) { 

    src_shadow += 12;
    _mm_prefetch((const char *)src_shadow, _MM_HINT_T0);

    // Now we have
    // v0.floats[0] = tmp[0][col=0][re]
    // v0.floats[1] = tmp[0][col=0][im]
    // v0.floats[2] = tmp[0][col=1][re]
    // v0.floats[3] = tmp[0][col=1][im]
    //
    // v1.floats[0] = tmp[0][col=2][re]
    // v1.floats[1] = tmp[0][col=2][im]
    // v1.floats[2] = tmp[1][col=0][re]
    // v1.floats[3] = tmp[1][col=0][im]
    //
    // v2.floats[0] = tmp[1][col=1][re]
    // v2.floats[1] = tmp[1][col=1][im]
    // v2.floats[2] = tmp[1][col=2][re]
    // v2.floats[3] = tmp[1][col=2][im]
    
    // This one is easy no shufs needed
    _mm_store_ps(dst_shadow, v0.vector);
    _mm_store_ps(dst_shadow+4, v1.vector);
    _mm_store_ps(dst_shadow+8, v2.vector);

    _mm_store_ps(dst_shadow+12, v0.vector);
    _mm_store_ps(dst_shadow+16, v1.vector);
    _mm_store_ps(dst_shadow+20, v2.vector);

    dst_shadow+=24;

    v0.vector = _mm_load_ps(src_shadow);
    v1.vector = _mm_load_ps(src_shadow+4);
    v2.vector = _mm_load_ps(src_shadow+8);

    
  }
  
  // Now we have
  // v0.floats[0] = tmp[0][col=0][re]
  // v0.floats[1] = tmp[0][col=0][im]
  // v0.floats[2] = tmp[0][col=1][re]
  // v0.floats[3] = tmp[0][col=1][im]
  //
  // v1.floats[0] = tmp[0][col=2][re]
  // v1.floats[1] = tmp[0][col=2][im]
  // v1.floats[2] = tmp[1][col=0][re]
  // v1.floats[3] = tmp[1][col=0][im]
  //
  // v2.floats[0] = tmp[1][col=1][re]
  // v2.floats[1] = tmp[1][col=1][im]
  // v2.floats[2] = tmp[1][col=2][re]
  // v2.floats[3] = tmp[1][col=2][im]
  
  // These are the top 2 components of the result
  // so I can store them out already.
  // Again no shufs needed
  _mm_store_ps(dst_shadow, v0.vector);
  _mm_store_ps(dst_shadow+4, v1.vector);
  _mm_store_ps(dst_shadow+8, v2.vector);
    
  // so I can  them out already - bypass cache.
  _mm_store_ps(dst_shadow+12, v0.vector);
  _mm_store_ps(dst_shadow+16, v1.vector);
  _mm_store_ps(dst_shadow+20, v2.vector);
}

/** \brief Spin recon (1/2)(1-\gamma3)
 *
 * \param src (pointer to 2 vector source)
 * \param dst (pointer to 4 vector dest)
 * \param n_vec (number of vectors to reconstruct)

 * It is assumeed that src points to an array of  2 spin
 * 3 colour vectors with complex components. 
 * It is assumed that dst points to an array of 4 spin
 * 3 colour vectors with complex components. In both cases
 * ordering is that spin is slowest, and complex components
 * are fastest. */
inline  
void inlineSpinReconDir3Minus(const REAL32* src, REAL32 *dst, unsigned int n_vec) 
{

#ifdef DEBUG_GENERIC_SPIN_RECON_INLINES
  QDPIO::cout << "inlineSpinReconDir0Minus" << endl;
#endif

  SSEVec v0, v1, v2, v7;
  const REAL32* src_shadow = src;
  REAL32* dst_shadow = dst;


  v7.floats[0] = -1;
  v7.floats[1] = -1;
  v7.floats[2] = -1;
  v7.floats[3] = -1;
  /*                              ( 1  0 -1  0)  ( a0 )    ( a0 - a2 )
   *  B  :=  ( 1 - Gamma  ) A  =  ( 0  1  0 -1)  ( a1 )  = ( a1 - a3 )
   *                    3         (-1  0  1  0)  ( a2 )    ( a2 - a0 )
   *                              ( 0 -1  0  1)  ( a3 )    ( a3 - a1 )
   
   * The bottom components of be may be reconstructed using the formula
   *  ( b2r + i b2i )  =  ( {a2r - a0r} + i{a2i - a0i} )  =  ( - b0r - i b0i )
   *  ( b3r + i b3i )     ( {a3r - a1r} + i{a3i - a1i} )     ( - b1r - i b1i ) 
   */    
  
  v0.vector = _mm_load_ps(src_shadow);
  v1.vector = _mm_load_ps(src_shadow+4);
  v2.vector = _mm_load_ps(src_shadow+8);

  for(unsigned int site=0; site < n_vec-1; site++) { 

    src_shadow += 12;
    _mm_prefetch((const char *)src_shadow, _MM_HINT_T0);

    // Now we have
    // v0.floats[0] = tmp[0][col=0][re]
    // v0.floats[1] = tmp[0][col=0][im]
    // v0.floats[2] = tmp[0][col=1][re]
    // v0.floats[3] = tmp[0][col=1][im]
    //
    // v1.floats[0] = tmp[0][col=2][re]
    // v1.floats[1] = tmp[0][col=2][im]
    // v1.floats[2] = tmp[1][col=0][re]
    // v1.floats[3] = tmp[1][col=0][im]
    //
    // v2.floats[0] = tmp[1][col=1][re]
    // v2.floats[1] = tmp[1][col=1][im]
    // v2.floats[2] = tmp[1][col=2][re]
    // v2.floats[3] = tmp[1][col=2][im]
    _mm_store_ps(dst_shadow, v0.vector);
    _mm_store_ps(dst_shadow+4, v1.vector);
    _mm_store_ps(dst_shadow+8, v2.vector);
    
    // This one is easy no shufs needed
    // multiply by -1,-1,-1,-1
    v0.vector = _mm_mul_ps(v0.vector, v7.vector);
    v1.vector = _mm_mul_ps(v1.vector, v7.vector);
    v2.vector = _mm_mul_ps(v2.vector, v7.vector);

    _mm_store_ps(dst_shadow+12, v0.vector);
    _mm_store_ps(dst_shadow+16, v1.vector);
    _mm_store_ps(dst_shadow+20, v2.vector);

    dst_shadow+=24;

    v0.vector = _mm_load_ps(src_shadow);
    v1.vector = _mm_load_ps(src_shadow+4);
    v2.vector = _mm_load_ps(src_shadow+8);

    
  }
  // Now we have
  // v0.floats[0] = tmp[0][col=0][re]
  // v0.floats[1] = tmp[0][col=0][im]
  // v0.floats[2] = tmp[0][col=1][re]
  // v0.floats[3] = tmp[0][col=1][im]
  //
  // v1.floats[0] = tmp[0][col=2][re]
  // v1.floats[1] = tmp[0][col=2][im]
  // v1.floats[2] = tmp[1][col=0][re]
  // v1.floats[3] = tmp[1][col=0][im]
  //
  // v2.floats[0] = tmp[1][col=1][re]
  // v2.floats[1] = tmp[1][col=1][im]
  // v2.floats[2] = tmp[1][col=2][re]
  // v2.floats[3] = tmp[1][col=2][im]
  _mm_store_ps(dst_shadow, v0.vector);
  _mm_store_ps(dst_shadow+4, v1.vector);
  _mm_store_ps(dst_shadow+8, v2.vector);
  
  // This one is easy no shufs needed
  // multiply by -1,-1,-1,-1
  v0.vector = _mm_mul_ps(v0.vector, v7.vector);
  v1.vector = _mm_mul_ps(v1.vector, v7.vector);
  v2.vector = _mm_mul_ps(v2.vector, v7.vector);
  
  _mm_store_ps(dst_shadow+12, v0.vector);
  _mm_store_ps(dst_shadow+16, v1.vector);
  _mm_store_ps(dst_shadow+20, v2.vector);
  
}



/** \brief Spin recon (1/2)(1+\gamma_0)
 *
 * \param src (pointer to 2 vector source)
 * \param dst (pointer to 4 vector dest)
 * \param n_vec (number of vectors to reconstruct)

 * It is assumeed that src points to an array of  2 spin
 * 3 colour vectors with complex components. 
 * It is assumed that dst points to an array of 4 spin
 * 3 colour vectors with complex components. In both cases
 * ordering is that spin is slowest, and complex components
 * are fastest. */
inline  
void inlineAddSpinReconDir0Plus(const REAL32* src, REAL32 *dst, unsigned int n_vec) 
{

#ifdef DEBUG_GENERIC_SPIN_RECON_INLINES
  QDPIO::cout << "inlinaAddSpinReconDir0Plus" << endl;
#endif

  /* 1 + \gamma_0 =  1  0  0  i 
                     0  1  i  0
                     0 -i  1  0
                    -i  0  0  1 
 
    *  ( b2r + i b2i )  =  ( {a2r + a1i} + i{a2i - a1r} )  =  ( b1i - i b1r )
    *  ( b3r + i b3i )     ( {a3r + a0i} + i{a3i - a0r} )     ( b0i - i b0r ) 
   */

  // 12 components in source, 24 in dest
  // 3 vectors for source

  const REAL32* src_shadow = src;
  REAL32* dst_shadow = dst;

  
  SSEVec v0, v1, v2, v3, v4, v5, v6;
  v6.floats[0] = +1;
  v6.floats[1] = -1;
  v6.floats[2] = +1;
  v6.floats[3] = -1;


  // Load source 
  v0.vector = _mm_load_ps(src_shadow);
  v1.vector = _mm_load_ps(src_shadow+4);
  v2.vector = _mm_load_ps(src_shadow+8);
  
  // Load top half of result
  v3.vector = _mm_load_ps(dst_shadow);
  v4.vector = _mm_load_ps(dst_shadow+4);
  v5.vector = _mm_load_ps(dst_shadow+8);

  // Add source to result 
  v3.vector = _mm_add_ps(v3.vector, v0.vector);
  v4.vector = _mm_add_ps(v4.vector, v1.vector);
  v5.vector = _mm_add_ps(v5.vector, v2.vector);

  // Store out result
  _mm_store_ps(dst_shadow, v3.vector);
  _mm_store_ps(dst_shadow+4, v4.vector);
  _mm_store_ps(dst_shadow+8, v5.vector);

  for(unsigned int site=0; site < n_vec-1; site++) { 


    src_shadow += 12;
    _mm_prefetch((const char *)src_shadow, _MM_HINT_T0);

    // Now the bottom half.  -- preload v4 and v5 
    v4.vector = _mm_load_ps(dst_shadow+12);
    v5.vector = _mm_load_ps(dst_shadow+16);
  
    // I want to set up 
    //
    //     v3[0] = v1[3] = tmp[1][col=0][im]
    //     v3[1] = v1[2] = tmp[1][col=0][re]
    //     v3[2] = v2[1] = tmp[1][col=1][im]
    //     v3[3] = v2[0] = tmp[1][col=1][re]
    // 
    // We can do this with v3 = v1,  v3 = shuf( v3, v1 ). Shuf code is 00 01 10 11 = x1B
    v3.vector = v1.vector;
    v3.vector = _mm_shuffle_ps(v3.vector, v2.vector, 0x1B);
    
    // Now I need to do multiply the mask (+1, -1, +1, -1)
    v3.vector = _mm_mul_ps(v3.vector, v6.vector);

    // Add to the destination
    v4.vector = _mm_add_ps(v3.vector, v4.vector);

    // v3 now free reuse it with the last element of the destination
    v3.vector = _mm_load_ps(dst_shadow+20);
    _mm_prefetch((const char *) dst_shadow+20, _MM_HINT_T0);
    // Store result
    _mm_store_ps(dst_shadow+12, v4.vector);
 

    // V4 is now free?
    // I want to set up 
    //
    //     v4[0] = v2[3] = tmp[1][col=2][im]
    //     v4[1] = v2[2] = tmp[1][col=2][re]
    //     v4[2] = v0[1] = tmp[0][col=0][im]
    //     v4[3] = v0[0] = tmp[0][col=0][re]
    // 
    // We can do this with v4 = v2,  v4 = shuf( v4, v2 ). Shuf code is 00 01 10 11 = x1B
    v4.vector = v2.vector;
    v4.vector = _mm_shuffle_ps(v4.vector, v0.vector, 0x1B);
    
    // Now I need to do multiply the mask (+1, -1, +1, -1)
    v4.vector = _mm_mul_ps(v4.vector, v6.vector);
    v5.vector = _mm_add_ps(v5.vector, v4.vector);
    _mm_store_ps(dst_shadow+16, v5.vector);

    // I want to set up 
    //
    //     v5[0] = v0[3] = tmp[0][col=1][im]
    //     v5[1] = v0[2] = tmp[0][col=1][re]
    //     v5[2] = v1[1] = tmp[0][col=2][im]
    //     v5[3] = v1[0] = tmp[0][col=2][re]
    // 
    // We can do this with v5 = v0,  v5 = shuf( v5, v1 ). Shuf code is 00 01 10 11 = x1B
    v5.vector = v0.vector;
    v5.vector = _mm_shuffle_ps(v5.vector, v1.vector, 0x1B);
    
    // Now I need to do multiply the mask (+1, -1, +1, -1)
    v5.vector = _mm_mul_ps(v5.vector, v6.vector);
    v3.vector = _mm_add_ps(v3.vector, v5.vector);
    _mm_store_ps(dst_shadow+20, v3.vector);

    dst_shadow+=24;

    // Load source 
    v0.vector = _mm_load_ps(src_shadow);
    v1.vector = _mm_load_ps(src_shadow+4);
    v2.vector = _mm_load_ps(src_shadow+8);
    
    // Load top half of result
    v3.vector = _mm_load_ps(dst_shadow);
    v4.vector = _mm_load_ps(dst_shadow+4);
    v5.vector = _mm_load_ps(dst_shadow+8);
    
    // Add source to result 
    v3.vector = _mm_add_ps(v3.vector, v0.vector);
    v4.vector = _mm_add_ps(v4.vector, v1.vector);
    v5.vector = _mm_add_ps(v5.vector, v2.vector);

    // Store out result
    _mm_store_ps(dst_shadow, v3.vector);
    _mm_store_ps(dst_shadow+4, v4.vector);
    _mm_store_ps(dst_shadow+8, v5.vector);
    
  }

  // Now the bottom half.  -- preload v4 and v5 
  v4.vector = _mm_load_ps(dst_shadow+12);
  v5.vector = _mm_load_ps(dst_shadow+16);
  
  // I want to set up 
  //
  //     v3[0] = v1[3] = tmp[1][col=0][im]
  //     v3[1] = v1[2] = tmp[1][col=0][re]
  //     v3[2] = v2[1] = tmp[1][col=1][im]
  //     v3[3] = v2[0] = tmp[1][col=1][re]
  // 
  // We can do this with v3 = v1,  v3 = shuf( v3, v1 ). Shuf code is 00 01 10 11 = x1B
  v3.vector = v1.vector;
  v3.vector = _mm_shuffle_ps(v3.vector, v2.vector, 0x1B);
  
  // Now I need to do multiply the mask (+1, -1, +1, -1)
  v3.vector = _mm_mul_ps(v3.vector, v6.vector);
  
  // Add to the destination
  v4.vector = _mm_add_ps(v3.vector, v4.vector);
  
  // v3 now free reuse it with the last element of the destination
  v3.vector = _mm_load_ps(dst_shadow+20);
  _mm_prefetch((const char *) dst_shadow+20, _MM_HINT_T0);
  
  // Store result
  _mm_store_ps(dst_shadow+12, v4.vector);
  
  
  // V4 is now free?
  // I want to set up 
  //
  //     v4[0] = v2[3] = tmp[1][col=2][im]
  //     v4[1] = v2[2] = tmp[1][col=2][re]
  //     v4[2] = v0[1] = tmp[0][col=0][im]
  //     v4[3] = v0[0] = tmp[0][col=0][re]
  // 
  // We can do this with v4 = v2,  v4 = shuf( v4, v2 ). Shuf code is 00 01 10 11 = x1B
  v4.vector = v2.vector;
  v4.vector = _mm_shuffle_ps(v4.vector, v0.vector, 0x1B);
  
  // Now I need to do multiply the mask (+1, -1, +1, -1)
  v4.vector = _mm_mul_ps(v4.vector, v6.vector);
  v5.vector = _mm_add_ps(v5.vector, v4.vector);
  _mm_store_ps(dst_shadow+16, v5.vector);
  
  // I want to set up 
  //
  //     v5[0] = v0[3] = tmp[0][col=1][im]
  //     v5[1] = v0[2] = tmp[0][col=1][re]
  //     v5[2] = v1[1] = tmp[0][col=2][im]
  //     v5[3] = v1[0] = tmp[0][col=2][re]
  // 
  // We can do this with v5 = v0,  v5 = shuf( v5, v1 ). Shuf code is 00 01 10 11 = x1B
  v5.vector = v0.vector;
  v5.vector = _mm_shuffle_ps(v5.vector, v1.vector, 0x1B);
  
  // Now I need to do multiply the mask (+1, -1, +1, -1)
  v5.vector = _mm_mul_ps(v5.vector, v6.vector);
  v3.vector = _mm_add_ps(v3.vector, v5.vector);
  _mm_store_ps(dst_shadow+20, v3.vector);

}

/** \brief Spin recon (1/2)(1-\gamma_0)
 *
 * \param src (pointer to 2 vector source)
 * \param dst (pointer to 4 vector dest)
 * \param n_vec (number of vectors to reconstruct)

 * It is assumeed that src points to an array of  2 spin
 * 3 colour vectors with complex components. 
 * It is assumed that dst points to an array of 4 spin
 * 3 colour vectors with complex components. In both cases
 * ordering is that spin is slowest, and complex components
 * are fastest. */
inline  
void inlineAddSpinReconDir0Minus(const REAL32* src, REAL32 *dst, unsigned int n_vec) 
{

#ifdef DEBUG_GENERIC_SPIN_RECON_INLINES
  QDPIO::cout << "inlinaAddSpinReconDir0Minus" << endl;
#endif


   /*                              ( 1  0  0 -i)  ( a0 )    ( a0 - i a3 )
    *  B  :=  ( 1 - Gamma  ) A  =  ( 0  1 -i  0)  ( a1 )  = ( a1 - i a2 )
    *                    0         ( 0  i  1  0)  ( a2 )    ( a2 + i a1 )
    *                              ( i  0  0  1)  ( a3 )    ( a3 + i a0 )
    
    * The bottom components of be may be reconstructed using the formula
    *   ( b2r + i b2i )  =  ( {a2r - a1i} + i{a2i + a1r} )  =  ( - b1i + i b1r )
    *   ( b3r + i b3i )     ( {a3r - a0i} + i{a3i + a0r} )     ( - b0i + i b0r ) 
    */

  const REAL32* src_shadow = src;
  REAL32* dst_shadow = dst;

  
  SSEVec v0, v1, v2, v3, v4, v5, v6;
  v6.floats[0] = -1;
  v6.floats[1] = +1;
  v6.floats[2] = -1;
  v6.floats[3] = +1;


  // Load source 
  v0.vector = _mm_load_ps(src_shadow);
  v1.vector = _mm_load_ps(src_shadow+4);
  v2.vector = _mm_load_ps(src_shadow+8);
  
  // Load top half of result
  v3.vector = _mm_load_ps(dst_shadow);
  v4.vector = _mm_load_ps(dst_shadow+4);
  v5.vector = _mm_load_ps(dst_shadow+8);

  // Add source to result 
  v3.vector = _mm_add_ps(v3.vector, v0.vector);
  v4.vector = _mm_add_ps(v4.vector, v1.vector);
  v5.vector = _mm_add_ps(v5.vector, v2.vector);

  // Store out result
  _mm_store_ps(dst_shadow, v3.vector);
  _mm_store_ps(dst_shadow+4, v4.vector);
  _mm_store_ps(dst_shadow+8, v5.vector);

  for(unsigned int site=0; site < n_vec-1; site++) { 


    src_shadow += 12;
    _mm_prefetch((const char *)src_shadow, _MM_HINT_T0);

    // Now the bottom half.  -- preload v4 and v5 
    v4.vector = _mm_load_ps(dst_shadow+12);
    v5.vector = _mm_load_ps(dst_shadow+16);
  
    // I want to set up 
    //
    //     v3[0] = v1[3] = tmp[1][col=0][im]
    //     v3[1] = v1[2] = tmp[1][col=0][re]
    //     v3[2] = v2[1] = tmp[1][col=1][im]
    //     v3[3] = v2[0] = tmp[1][col=1][re]
    // 
    // We can do this with v3 = v1,  v3 = shuf( v3, v1 ). Shuf code is 00 01 10 11 = x1B
    v3.vector = v1.vector;
    v3.vector = _mm_shuffle_ps(v3.vector, v2.vector, 0x1B);
    
    // Now I need to do multiply the mask (+1, -1, +1, -1)
    v3.vector = _mm_mul_ps(v3.vector, v6.vector);

    // Add to the destination
    v4.vector = _mm_add_ps(v3.vector, v4.vector);

    // v3 now free reuse it with the last element of the destination
    v3.vector = _mm_load_ps(dst_shadow+20);
    _mm_prefetch((const char *) dst_shadow+20, _MM_HINT_T0);

    // Store result
    _mm_store_ps(dst_shadow+12, v4.vector);
 

    // V4 is now free?
    // I want to set up 
    //
    //     v4[0] = v2[3] = tmp[1][col=2][im]
    //     v4[1] = v2[2] = tmp[1][col=2][re]
    //     v4[2] = v0[1] = tmp[0][col=0][im]
    //     v4[3] = v0[0] = tmp[0][col=0][re]
    // 
    // We can do this with v4 = v2,  v4 = shuf( v4, v2 ). Shuf code is 00 01 10 11 = x1B
    v4.vector = v2.vector;
    v4.vector = _mm_shuffle_ps(v4.vector, v0.vector, 0x1B);
    
    // Now I need to do multiply the mask (+1, -1, +1, -1)
    v4.vector = _mm_mul_ps(v4.vector, v6.vector);
    v5.vector = _mm_add_ps(v5.vector, v4.vector);
    _mm_store_ps(dst_shadow+16, v5.vector);

    // I want to set up 
    //
    //     v5[0] = v0[3] = tmp[0][col=1][im]
    //     v5[1] = v0[2] = tmp[0][col=1][re]
    //     v5[2] = v1[1] = tmp[0][col=2][im]
    //     v5[3] = v1[0] = tmp[0][col=2][re]
    // 
    // We can do this with v5 = v0,  v5 = shuf( v5, v1 ). Shuf code is 00 01 10 11 = x1B
    v5.vector = v0.vector;
    v5.vector = _mm_shuffle_ps(v5.vector, v1.vector, 0x1B);
    
    // Now I need to do multiply the mask (+1, -1, +1, -1)
    v5.vector = _mm_mul_ps(v5.vector, v6.vector);
    v3.vector = _mm_add_ps(v3.vector, v5.vector);
    _mm_store_ps(dst_shadow+20, v3.vector);

    dst_shadow+=24;

    // Load source 
    v0.vector = _mm_load_ps(src_shadow);
    v1.vector = _mm_load_ps(src_shadow+4);
    v2.vector = _mm_load_ps(src_shadow+8);
    
    // Load top half of result
    v3.vector = _mm_load_ps(dst_shadow);
    v4.vector = _mm_load_ps(dst_shadow+4);
    v5.vector = _mm_load_ps(dst_shadow+8);
    
    // Add source to result 
    v3.vector = _mm_add_ps(v3.vector, v0.vector);
    v4.vector = _mm_add_ps(v4.vector, v1.vector);
    v5.vector = _mm_add_ps(v5.vector, v2.vector);

    // Store out result
    _mm_store_ps(dst_shadow, v3.vector);
    _mm_store_ps(dst_shadow+4, v4.vector);
    _mm_store_ps(dst_shadow+8, v5.vector);
    
  }

  // Now the bottom half.  -- preload v4 and v5 
  v4.vector = _mm_load_ps(dst_shadow+12);
  v5.vector = _mm_load_ps(dst_shadow+16);
  
  // I want to set up 
  //
  //     v3[0] = v1[3] = tmp[1][col=0][im]
  //     v3[1] = v1[2] = tmp[1][col=0][re]
  //     v3[2] = v2[1] = tmp[1][col=1][im]
  //     v3[3] = v2[0] = tmp[1][col=1][re]
  // 
  // We can do this with v3 = v1,  v3 = shuf( v3, v1 ). Shuf code is 00 01 10 11 = x1B
  v3.vector = v1.vector;
  v3.vector = _mm_shuffle_ps(v3.vector, v2.vector, 0x1B);
  
  // Now I need to do multiply the mask (+1, -1, +1, -1)
  v3.vector = _mm_mul_ps(v3.vector, v6.vector);
  
  // Add to the destination
  v4.vector = _mm_add_ps(v3.vector, v4.vector);
  
  // v3 now free reuse it with the last element of the destination
  v3.vector = _mm_load_ps(dst_shadow+20);
  _mm_prefetch((const char *) dst_shadow+20, _MM_HINT_T0);
  
  // Store result
  _mm_store_ps(dst_shadow+12, v4.vector);
  
  
  // V4 is now free?
  // I want to set up 
  //
  //     v4[0] = v2[3] = tmp[1][col=2][im]
  //     v4[1] = v2[2] = tmp[1][col=2][re]
  //     v4[2] = v0[1] = tmp[0][col=0][im]
  //     v4[3] = v0[0] = tmp[0][col=0][re]
  // 
  // We can do this with v4 = v2,  v4 = shuf( v4, v2 ). Shuf code is 00 01 10 11 = x1B
  v4.vector = v2.vector;
  v4.vector = _mm_shuffle_ps(v4.vector, v0.vector, 0x1B);
  
  // Now I need to do multiply the mask (+1, -1, +1, -1)
  v4.vector = _mm_mul_ps(v4.vector, v6.vector);
  v5.vector = _mm_add_ps(v5.vector, v4.vector);
  _mm_store_ps(dst_shadow+16, v5.vector);
  
  // I want to set up 
  //
  //     v5[0] = v0[3] = tmp[0][col=1][im]
  //     v5[1] = v0[2] = tmp[0][col=1][re]
  //     v5[2] = v1[1] = tmp[0][col=2][im]
  //     v5[3] = v1[0] = tmp[0][col=2][re]
  // 
  // We can do this with v5 = v0,  v5 = shuf( v5, v1 ). Shuf code is 00 01 10 11 = x1B
  v5.vector = v0.vector;
  v5.vector = _mm_shuffle_ps(v5.vector, v1.vector, 0x1B);
  
  // Now I need to do multiply the mask (+1, -1, +1, -1)
  v5.vector = _mm_mul_ps(v5.vector, v6.vector);
  v3.vector = _mm_add_ps(v3.vector, v5.vector);
  _mm_store_ps(dst_shadow+20, v3.vector);

  
  
}

/** \brief Spin recon (1/2)(1+\gamma_1)
 *
 * \param src (pointer to 2 vector source)
 * \param dst (pointer to 4 vector dest)
 * \param n_vec (number of vectors to reconstruct)

 * It is assumeed that src points to an array of  2 spin
 * 3 colour vectors with complex components. 
 * It is assumed that dst points to an array of 4 spin
 * 3 colour vectors with complex components. In both cases
 * ordering is that spin is slowest, and complex components
 * are fastest. */
inline  
void inlineAddSpinReconDir1Plus(const REAL32* src, REAL32 *dst, unsigned int n_vec) 
{

#ifdef DEBUG_GENERIC_SPIN_RECON_INLINES
  QDPIO::cout << "inlinaAddSpinReconDir0Plus" << endl;
#endif

  const REAL32* src_shadow = src;
  REAL32* dst_shadow = dst;

    /* 1 + \gamma_1 =  1  0  0 -1 
     *                 0  1  1  0
     *                 0  1  1  0
     *                -1  0  0  1 
 
     *   ( b2r + i b2i )  =  ( {a2r + a1r} + i{a2i + a1i} )  =  (   b1r + i b1i )
     *   ( b3r + i b3i )     ( {a3r - a0r} + i{a3i - a0i} )     ( - b0r - i b0i ) 
  
    */
  
  
  SSEVec v0, v1, v2, v3, v4, v5, v6;
  v6.floats[0] = +1;
  v6.floats[1] = +1;
  v6.floats[2] = -1;
  v6.floats[3] = -1;


  // Load source 
  v0.vector = _mm_load_ps(src_shadow);
  v1.vector = _mm_load_ps(src_shadow+4);
  v2.vector = _mm_load_ps(src_shadow+8);
  
  // Load top half of result
  v3.vector = _mm_load_ps(dst_shadow);
  v4.vector = _mm_load_ps(dst_shadow+4);
  v5.vector = _mm_load_ps(dst_shadow+8);

  // Add source to result 
  v3.vector = _mm_add_ps(v3.vector, v0.vector);
  v4.vector = _mm_add_ps(v4.vector, v1.vector);
  v5.vector = _mm_add_ps(v5.vector, v2.vector);

  // Store out result
  _mm_store_ps(dst_shadow, v3.vector);
  _mm_store_ps(dst_shadow+4, v4.vector);
  _mm_store_ps(dst_shadow+8, v5.vector);

  for(unsigned int site=0; site < n_vec-1; site++) { 


    src_shadow += 12;
    _mm_prefetch((const char *)src_shadow, _MM_HINT_T0);

    // Now the bottom half.  -- preload v4 and v5 
    v4.vector = _mm_load_ps(dst_shadow+12);
    v5.vector = _mm_load_ps(dst_shadow+16);
 
    // I want to set up 
    //
    //     v3[0] = v1[2] = tmp[1][col=0][re]
    //     v3[1] = v1[3] = tmp[1][col=0][im]
    //     v3[2] = v2[0] = tmp[1][col=1][re]
    //     v3[3] = v2[1] = tmp[1][col=1][im]
    // 
    // We can do this with v3 = v1,  v3 = shuf( v3, v2 ). Shuf code is 01 00 11 10 = x4E
    v3.vector = v1.vector;
    v3.vector = _mm_shuffle_ps(v3.vector, v2.vector, 0x4E);
    
    // Add to result
    v4.vector = _mm_add_ps(v4.vector, v3.vector);
    _mm_store_ps( dst_shadow+12, v4.vector);

    v3.vector = _mm_load_ps(dst_shadow+20);
    _mm_prefetch((const char *) dst_shadow+20, _MM_HINT_T0);

    // I want to set up 
    //
    //     v4[0] = v2[2] = tmp[1][col=2][re]
    //     v4[1] = v2[3] = tmp[1][col=2][im]
    //     v4[2] = v0[0] = tmp[0][col=0][re]
    //     v4[3] = v0[1] = tmp[0][col=0][im]
    // 
    // We can do this with v4 = v2,  v4 = shuf( v4, v0 ). Shuf code is 01 00 11 10 = x4E
    v4.vector = v2.vector;
    v4.vector = _mm_shuffle_ps(v4.vector, v0.vector, 0x4E);

    // Need to multiply in mask (+1,+1,-1,-1) 
    v4.vector = _mm_mul_ps(v4.vector, v6.vector);
   
    // Add to result 
    v5.vector = _mm_add_ps(v5.vector, v4.vector);
    _mm_store_ps( dst_shadow+16, v5.vector);

    // I want to set up 
    //
    //     v5[0] = v0[2] = tmp[0][col=1][re]
    //     v5[1] = v0[3] = tmp[0][col=1][im]
    //     v5[2] = v1[0] = tmp[0][col=2][re]
    //     v5[3] = v1[1] = tmp[0][col=2][im]
    // 
    // We can do this with v5 = v0,  v5 = shuf( v5, v1 ). Shuf code is 01 00 11 10 = x4E
    v5.vector = v0.vector;
    v5.vector = _mm_shuffle_ps(v5.vector, v1.vector, 0x4E);
    
    // Now I need to do multiply the mask (-1,-1,-1,-11)
    v3.vector = _mm_sub_ps(v3.vector, v5.vector);

    // so I can  them out already - bypass cache.
    _mm_store_ps(dst_shadow+20, v3.vector);
 
    dst_shadow+=24;

    // Load source 
    v0.vector = _mm_load_ps(src_shadow);
    v1.vector = _mm_load_ps(src_shadow+4);
    v2.vector = _mm_load_ps(src_shadow+8);
    
    // Load top half of result
    v3.vector = _mm_load_ps(dst_shadow);
    v4.vector = _mm_load_ps(dst_shadow+4);
    v5.vector = _mm_load_ps(dst_shadow+8);
    
    // Add source to result 
    v3.vector = _mm_add_ps(v3.vector, v0.vector);
    v4.vector = _mm_add_ps(v4.vector, v1.vector);
    v5.vector = _mm_add_ps(v5.vector, v2.vector);

    // Store out result
    _mm_store_ps(dst_shadow, v3.vector);
    _mm_store_ps(dst_shadow+4, v4.vector);
    _mm_store_ps(dst_shadow+8, v5.vector);
    
  }

    // Now the bottom half.  -- preload v4 and v5 
    v4.vector = _mm_load_ps(dst_shadow+12);
    v5.vector = _mm_load_ps(dst_shadow+16);
 
    // I want to set up 
    //
    //     v3[0] = v1[2] = tmp[1][col=0][re]
    //     v3[1] = v1[3] = tmp[1][col=0][im]
    //     v3[2] = v2[0] = tmp[1][col=1][re]
    //     v3[3] = v2[1] = tmp[1][col=1][im]
    // 
    // We can do this with v3 = v1,  v3 = shuf( v3, v2 ). Shuf code is 01 00 11 10 = x4E
    v3.vector = v1.vector;
    v3.vector = _mm_shuffle_ps(v3.vector, v2.vector, 0x4E);
    
    // Add to result
    v4.vector = _mm_add_ps(v4.vector, v3.vector);
    _mm_store_ps( dst_shadow+12, v4.vector);

    v3.vector = _mm_load_ps(dst_shadow+20);
    _mm_prefetch((const char *) dst_shadow+20, _MM_HINT_T0);
    // I want to set up 
    //
    //     v4[0] = v2[2] = tmp[1][col=2][re]
    //     v4[1] = v2[3] = tmp[1][col=2][im]
    //     v4[2] = v0[0] = tmp[0][col=0][re]
    //     v4[3] = v0[1] = tmp[0][col=0][im]
    // 
    // We can do this with v4 = v2,  v4 = shuf( v4, v0 ). Shuf code is 01 00 11 10 = x4E
    v4.vector = v2.vector;
    v4.vector = _mm_shuffle_ps(v4.vector, v0.vector, 0x4E);

    // Need to multiply in mask (+1,+1,-1,-1) 
    v4.vector = _mm_mul_ps(v4.vector, v6.vector);
   
    // Add to result 
    v5.vector = _mm_add_ps(v5.vector, v4.vector);
    _mm_store_ps( dst_shadow+16, v5.vector);

    // I want to set up 
    //
    //     v5[0] = v0[2] = tmp[0][col=1][re]
    //     v5[1] = v0[3] = tmp[0][col=1][im]
    //     v5[2] = v1[0] = tmp[0][col=2][re]
    //     v5[3] = v1[1] = tmp[0][col=2][im]
    // 
    // We can do this with v5 = v0,  v5 = shuf( v5, v1 ). Shuf code is 01 00 11 10 = x4E
    v5.vector = v0.vector;
    v5.vector = _mm_shuffle_ps(v5.vector, v1.vector, 0x4E);
    
    // Now I need to do multiply the mask (-1,-1,-1,-11)
    v3.vector = _mm_sub_ps(v3.vector, v5.vector);

    // so I can  them out already - bypass cache.
    _mm_store_ps(dst_shadow+20, v3.vector);

#if 0

  REAL32 tmp_hspinor[4][3][2];

  const REAL32* src_shadow = src;
  REAL32* dst_shadow = dst;

  const int re = 0;
  const int im = 1;
  const int Ncmpx = 2;
  const int Nsby2 = 2;

  for(int site=0; site < n_vec; site++) {
    REAL32* tmp_shadow = &(tmp_hspinor[0][0][0]);
    
    // Store in the half spinor - write out the first two components
    for(int store=0; store < Nsby2*Nc*Ncmpx; store++) {
      REAL32 tmp = *(src_shadow++);
      *(tmp_shadow++) = tmp;
      *(dst_shadow++) += tmp;
    }
    
    // Reconstruct Spin 2
    for(int col=0; col < Nc; col++) { 
      *(dst_shadow++) += tmp_hspinor[1][col][re]; 
      *(dst_shadow++) += tmp_hspinor[1][col][im];
    }

    // Reconstruct Spin 3;
    for(int col=0; col < Nc; col++) { 
      *(dst_shadow++) -= tmp_hspinor[0][col][re];
      *(dst_shadow++) -= tmp_hspinor[0][col][im];
    }


  }
#endif
}

/** \brief Spin recon (1/2)(1-\gamma_1)
 *
 * \param src (pointer to 2 vector source)
 * \param dst (pointer to 4 vector dest)
 * \param n_vec (number of vectors to reconstruct)

 * It is assumeed that src points to an array of  2 spin
 * 3 colour vectors with complex components. 
 * It is assumed that dst points to an array of 4 spin
 * 3 colour vectors with complex components. In both cases
 * ordering is that spin is slowest, and complex components
 * are fastest. */
inline  
void inlineAddSpinReconDir1Minus(const REAL32* src, REAL32 *dst, unsigned int n_vec) 
{

#ifdef DEBUG_GENERIC_SPIN_RECON_INLINES
  QDPIO::cout << "inlinaAddSpinReconDir0Minus" << endl;
#endif

  const REAL32* src_shadow = src;
  REAL32* dst_shadow = dst;

    /* 1 + \gamma_1 =  1  0  0 -1 
     *                 0  1  1  0
     *                 0  1  1  0
     *                -1  0  0  1 
 
     *   ( b2r + i b2i )  =  ( {a2r + a1r} + i{a2i + a1i} )  =  (   b1r + i b1i )
     *   ( b3r + i b3i )     ( {a3r - a0r} + i{a3i - a0i} )     ( - b0r - i b0i ) 
  
    */
  
  
  SSEVec v0, v1, v2, v3, v4, v5, v6;
  v6.floats[0] = -1;
  v6.floats[1] = -1;
  v6.floats[2] = +1;
  v6.floats[3] = +1;


  // Load source 
  v0.vector = _mm_load_ps(src_shadow);
  v1.vector = _mm_load_ps(src_shadow+4);
  v2.vector = _mm_load_ps(src_shadow+8);
  
  // Load top half of result
  v3.vector = _mm_load_ps(dst_shadow);
  v4.vector = _mm_load_ps(dst_shadow+4);
  v5.vector = _mm_load_ps(dst_shadow+8);

  // Add source to result 
  v3.vector = _mm_add_ps(v3.vector, v0.vector);
  v4.vector = _mm_add_ps(v4.vector, v1.vector);
  v5.vector = _mm_add_ps(v5.vector, v2.vector);

  // Store out result
  _mm_store_ps(dst_shadow, v3.vector);
  _mm_store_ps(dst_shadow+4, v4.vector);
  _mm_store_ps(dst_shadow+8, v5.vector);

  for(unsigned int site=0; site < n_vec-1; site++) { 

    src_shadow += 12;
    _mm_prefetch((const char *)src_shadow, _MM_HINT_T0);

    // Now the bottom half.  -- preload v4 and v5 
    v4.vector = _mm_load_ps(dst_shadow+12);
    v5.vector = _mm_load_ps(dst_shadow+16);
 
    // I want to set up 
    //
    //     v3[0] = v1[2] = tmp[1][col=0][re]
    //     v3[1] = v1[3] = tmp[1][col=0][im]
    //     v3[2] = v2[0] = tmp[1][col=1][re]
    //     v3[3] = v2[1] = tmp[1][col=1][im]
    // 
    // We can do this with v3 = v1,  v3 = shuf( v3, v2 ). Shuf code is 01 00 11 10 = x4E
    v3.vector = v1.vector;
    v3.vector = _mm_shuffle_ps(v3.vector, v2.vector, 0x4E);
    
    // Add to result
    v4.vector = _mm_sub_ps(v4.vector, v3.vector);
    _mm_store_ps( dst_shadow+12, v4.vector);

    v3.vector = _mm_load_ps(dst_shadow+20);
    _mm_prefetch((const char *) dst_shadow+20, _MM_HINT_T0);
    // I want to set up 
    //
    //     v4[0] = v2[2] = tmp[1][col=2][re]
    //     v4[1] = v2[3] = tmp[1][col=2][im]
    //     v4[2] = v0[0] = tmp[0][col=0][re]
    //     v4[3] = v0[1] = tmp[0][col=0][im]
    // 
    // We can do this with v4 = v2,  v4 = shuf( v4, v0 ). Shuf code is 01 00 11 10 = x4E
    v4.vector = v2.vector;
    v4.vector = _mm_shuffle_ps(v4.vector, v0.vector, 0x4E);

    // Need to multiply in mask (+1,+1,-1,-1) 
    v4.vector = _mm_mul_ps(v4.vector, v6.vector);
   
    // Add to result 
    v5.vector = _mm_add_ps(v5.vector, v4.vector);
    _mm_store_ps( dst_shadow+16, v5.vector);

    // I want to set up 
    //
    //     v5[0] = v0[2] = tmp[0][col=1][re]
    //     v5[1] = v0[3] = tmp[0][col=1][im]
    //     v5[2] = v1[0] = tmp[0][col=2][re]
    //     v5[3] = v1[1] = tmp[0][col=2][im]
    // 
    // We can do this with v5 = v0,  v5 = shuf( v5, v1 ). Shuf code is 01 00 11 10 = x4E
    v5.vector = v0.vector;
    v5.vector = _mm_shuffle_ps(v5.vector, v1.vector, 0x4E);
    
    // Now I need to do multiply the mask (-1,-1,-1,-11)
    v3.vector = _mm_add_ps(v3.vector, v5.vector);

    // so I can  them out already - bypass cache.
    _mm_store_ps(dst_shadow+20, v3.vector);
 
    dst_shadow+=24;

    // Load source 
    v0.vector = _mm_load_ps(src_shadow);
    v1.vector = _mm_load_ps(src_shadow+4);
    v2.vector = _mm_load_ps(src_shadow+8);
    
    // Load top half of result
    v3.vector = _mm_load_ps(dst_shadow);
    v4.vector = _mm_load_ps(dst_shadow+4);
    v5.vector = _mm_load_ps(dst_shadow+8);
    
    // Add source to result 
    v3.vector = _mm_add_ps(v3.vector, v0.vector);
    v4.vector = _mm_add_ps(v4.vector, v1.vector);
    v5.vector = _mm_add_ps(v5.vector, v2.vector);

    // Store out result
    _mm_store_ps(dst_shadow, v3.vector);
    _mm_store_ps(dst_shadow+4, v4.vector);
    _mm_store_ps(dst_shadow+8, v5.vector);
    
  }

  // Now the bottom half.  -- preload v4 and v5 
  v4.vector = _mm_load_ps(dst_shadow+12);
  v5.vector = _mm_load_ps(dst_shadow+16);
  
  // I want to set up 
  //
  //     v3[0] = v1[2] = tmp[1][col=0][re]
  //     v3[1] = v1[3] = tmp[1][col=0][im]
  //     v3[2] = v2[0] = tmp[1][col=1][re]
  //     v3[3] = v2[1] = tmp[1][col=1][im]
  // 
  // We can do this with v3 = v1,  v3 = shuf( v3, v2 ). Shuf code is 01 00 11 10 = x4E
  v3.vector = v1.vector;
  v3.vector = _mm_shuffle_ps(v3.vector, v2.vector, 0x4E);
  
  // Add to result
  v4.vector = _mm_sub_ps(v4.vector, v3.vector);
  _mm_store_ps( dst_shadow+12, v4.vector);
  
  v3.vector = _mm_load_ps(dst_shadow+20);
  _mm_prefetch((const char *) dst_shadow+20, _MM_HINT_T0);
  // I want to set up 
  //
  //     v4[0] = v2[2] = tmp[1][col=2][re]
  //     v4[1] = v2[3] = tmp[1][col=2][im]
  //     v4[2] = v0[0] = tmp[0][col=0][re]
  //     v4[3] = v0[1] = tmp[0][col=0][im]
  // 
  // We can do this with v4 = v2,  v4 = shuf( v4, v0 ). Shuf code is 01 00 11 10 = x4E
  v4.vector = v2.vector;
  v4.vector = _mm_shuffle_ps(v4.vector, v0.vector, 0x4E);
  
  // Need to multiply in mask (+1,+1,-1,-1) 
  v4.vector = _mm_mul_ps(v4.vector, v6.vector);
  
  // Add to result 
  v5.vector = _mm_add_ps(v5.vector, v4.vector);
  _mm_store_ps( dst_shadow+16, v5.vector);
  
  // I want to set up 
  //
  //     v5[0] = v0[2] = tmp[0][col=1][re]
  //     v5[1] = v0[3] = tmp[0][col=1][im]
  //     v5[2] = v1[0] = tmp[0][col=2][re]
  //     v5[3] = v1[1] = tmp[0][col=2][im]
  // 
  // We can do this with v5 = v0,  v5 = shuf( v5, v1 ). Shuf code is 01 00 11 10 = x4E
  v5.vector = v0.vector;
  v5.vector = _mm_shuffle_ps(v5.vector, v1.vector, 0x4E);
  
  // Now I need to do multiply the mask (-1,-1,-1,-11)
  v3.vector = _mm_add_ps(v3.vector, v5.vector);
  
  // so I can  them out already - bypass cache.
  _mm_store_ps(dst_shadow+20, v3.vector);

}


/** \brief Spin recon (1/2)(1+\gamma_2)
 *
 * \param src (pointer to 2 vector source)
 * \param dst (pointer to 4 vector dest)
 * \param n_vec (number of vectors to reconstruct)

 * It is assumeed that src points to an array of  2 spin
 * 3 colour vectors with complex components. 
 * It is assumed that dst points to an array of 4 spin
 * 3 colour vectors with complex components. In both cases
 * ordering is that spin is slowest, and complex components
 * are fastest. */
inline  
void inlineAddSpinReconDir2Plus(const REAL32* src, REAL32 *dst, unsigned int n_vec) 
{

#ifdef DEBUG_GENERIC_SPIN_RECON_INLINES
  QDPIO::cout << "inlinaAddSpinReconDir0Plus" << endl;
#endif


  const REAL32* src_shadow = src;
  REAL32* dst_shadow = dst;


  /* 1 + \gamma_2 =  1  0  i  0 
   *                 0  1  0 -i
   *                -i  0  1  0
   *                 0  i  0  1 
   *	     
   *  ( b2r + i b2i )  =  ( {a2r + a0i} + i{a2i - a0r} )  =  (   b0i - i b0r )
   *  ( b3r + i b3i )     ( {a3r - a1i} + i{a3i + a1r} )     ( - b1i + i b1r ) 
  */    
  SSEVec v0, v1, v2, v3, v4, v5, v6, v7;
  v6.floats[0] = +1;
  v6.floats[1] = -1;
  v6.floats[2] = +1;
  v6.floats[3] = -1;

  v7.floats[0] = +1;
  v7.floats[1] = -1;
  v7.floats[2] = -1;
  v7.floats[3] = +1;

  // Load source 
  v0.vector = _mm_load_ps(src_shadow);
  v1.vector = _mm_load_ps(src_shadow+4);
  v2.vector = _mm_load_ps(src_shadow+8);
  
  // Load top half of result
  v3.vector = _mm_load_ps(dst_shadow);
  v4.vector = _mm_load_ps(dst_shadow+4);
  v5.vector = _mm_load_ps(dst_shadow+8);

  // Add source to result 
  v3.vector = _mm_add_ps(v3.vector, v0.vector);
  v4.vector = _mm_add_ps(v4.vector, v1.vector);
  v5.vector = _mm_add_ps(v5.vector, v2.vector);

  // Store out result
  _mm_store_ps(dst_shadow, v3.vector);
  _mm_store_ps(dst_shadow+4, v4.vector);
  _mm_store_ps(dst_shadow+8, v5.vector);

  for(unsigned int site=0; site < n_vec-1; site++) { 
    
    src_shadow += 12;
    _mm_prefetch((const char *)src_shadow, _MM_HINT_T0);

    // Now the bottom half.  -- preload v4 and v5 
    v3.vector = _mm_load_ps(dst_shadow+12);
    v4.vector = _mm_load_ps(dst_shadow+16);
    v5.vector = _mm_load_ps(dst_shadow+20);
    _mm_prefetch((const char *) dst_shadow+20, _MM_HINT_T0);

    // I want to set up 
    //
    //     v0[0] = v0[1] = tmp[0][col=0][im]
    //     v0[1] = v0[0] = tmp[0][col=0][re]
    //     v0[2] = v0[3] = tmp[0][col=1][im]
    //     v0[3] = v0[2] = tmp[0][col=1][re]
    // 
    // We can do this with v0 = shuf( v0, v0 ). Shuf code is 10 11 00 01 = xB1
    v0.vector = _mm_shuffle_ps(v0.vector, v0.vector, 0xB1);
    
    // Mask (+1,-1,+1,-1)
    v0.vector = _mm_mul_ps(v0.vector, v6.vector);
    v3.vector = _mm_add_ps(v3.vector, v0.vector);
    _mm_store_ps(dst_shadow+12, v3.vector);
    

    // I want to set up 
    //
    //     v1[0] = v1[1] = tmp[0][col=2][im]
    //     v1[1] = v1[0] = tmp[0][col=2][re]
    //     v1[2] = v1[3] = tmp[1][col=0][im]
    //     v1[3] = v1[2] = tmp[1][col=0][re]
    // 
    // We can do this with v1 = shuf( v1, v1 ). Shuf code is 10 11 00 01 = xB1
    v1.vector = _mm_shuffle_ps(v1.vector, v1.vector, 0xB1);

    // Need to multiply in mask (+1,-1,-1,+1) 
    v1.vector = _mm_mul_ps(v1.vector, v7.vector);
    v4.vector = _mm_add_ps(v4.vector, v1.vector);
    _mm_store_ps(dst_shadow+16, v4.vector);

    // I want to set up 
    //
    //     v2[0] = v2[1] = tmp[1][col=1][im]
    //     v2[1] = v2[0] = tmp[1][col=1][re]
    //     v2[2] = v2[2] = tmp[1][col=2][im]
    //     v2[3] = v2[3] = tmp[1][col=2][re]
    // 
    // We can do this with v2 = shuf( v2, v2 ). Shuf code is 10 11 00 01 = xB1
    v2.vector = _mm_shuffle_ps(v2.vector, v2.vector, 0xB1);

    // Mask in -1, +1, -1, +1 from v7
    v2.vector = _mm_mul_ps(v2.vector, v6.vector);
    v5.vector = _mm_sub_ps(v5.vector, v2.vector);
    _mm_store_ps(dst_shadow+20, v5.vector);

    dst_shadow+=24;

    // Load source 
    v0.vector = _mm_load_ps(src_shadow);
    v1.vector = _mm_load_ps(src_shadow+4);
    v2.vector = _mm_load_ps(src_shadow+8);
    
    // Load top half of result
    v3.vector = _mm_load_ps(dst_shadow);
    v4.vector = _mm_load_ps(dst_shadow+4);
    v5.vector = _mm_load_ps(dst_shadow+8);
    
    // Add source to result 
    v3.vector = _mm_add_ps(v3.vector, v0.vector);
    v4.vector = _mm_add_ps(v4.vector, v1.vector);
    v5.vector = _mm_add_ps(v5.vector, v2.vector);

    // Store out result
    _mm_store_ps(dst_shadow, v3.vector);
    _mm_store_ps(dst_shadow+4, v4.vector);
    _mm_store_ps(dst_shadow+8, v5.vector);
    
  }

  // Now the bottom half.  -- preload v4 and v5 
  // Now the bottom half.  -- preload v4 and v5 
  v3.vector = _mm_load_ps(dst_shadow+12);
  v4.vector = _mm_load_ps(dst_shadow+16);
  v5.vector = _mm_load_ps(dst_shadow+20);
  _mm_prefetch((const char *) dst_shadow+20, _MM_HINT_T0);
  
  // I want to set up 
  //
  //     v0[0] = v0[1] = tmp[0][col=0][im]
  //     v0[1] = v0[0] = tmp[0][col=0][re]
  //     v0[2] = v0[3] = tmp[0][col=1][im]
  //     v0[3] = v0[2] = tmp[0][col=1][re]
  // 
  // We can do this with v0 = shuf( v0, v0 ). Shuf code is 10 11 00 01 = xB1
  v0.vector = _mm_shuffle_ps(v0.vector, v0.vector, 0xB1);
  
  // Mask (+1,-1,+1,-1)
  v0.vector = _mm_mul_ps(v0.vector, v6.vector);
  v3.vector = _mm_add_ps(v3.vector, v0.vector);
  _mm_store_ps(dst_shadow+12,v3.vector);
  
  
  // I want to set up 
  //
  //     v1[0] = v1[1] = tmp[0][col=2][im]
  //     v1[1] = v1[0] = tmp[0][col=2][re]
  //     v1[2] = v1[3] = tmp[1][col=0][im]
  //     v1[3] = v1[2] = tmp[1][col=0][re]
  // 
  // We can do this with v1 = shuf( v1, v1 ). Shuf code is 10 11 00 01 = xB1
  v1.vector = _mm_shuffle_ps(v1.vector, v1.vector, 0xB1);
  
  // Need to multiply in mask (+1,-1,-1,+1) 
  v1.vector = _mm_mul_ps(v1.vector, v7.vector);
  v4.vector = _mm_add_ps(v4.vector, v1.vector);
  _mm_store_ps(dst_shadow+16, v4.vector);

  // I want to set up 
  //
  //     v2[0] = v2[1] = tmp[1][col=1][im]
  //     v2[1] = v2[0] = tmp[1][col=1][re]
  //     v2[2] = v2[2] = tmp[1][col=2][im]
  //     v2[3] = v2[3] = tmp[1][col=2][re]
  // 
  // We can do this with v2 = shuf( v2, v2 ). Shuf code is 10 11 00 01 = xB1
  v2.vector = _mm_shuffle_ps(v2.vector, v2.vector, 0xB1);
  
  // Mask in -1, +1, -1, +1 from v7
  v2.vector = _mm_mul_ps(v2.vector, v6.vector);
  v5.vector = _mm_sub_ps(v5.vector, v2.vector);
  _mm_store_ps(dst_shadow+20, v5.vector);

}

/** \brief Spin recon (1/2)(1-\gamma_2)
 *
 * \param src (pointer to 2 vector source)
 * \param dst (pointer to 4 vector dest)
 * \param n_vec (number of vectors to reconstruct)

 * It is assumeed that src points to an array of  2 spin
 * 3 colour vectors with complex components. 
 * It is assumed that dst points to an array of 4 spin
 * 3 colour vectors with complex components. In both cases
 * ordering is that spin is slowest, and complex components
 * are fastest. */
inline  
void inlineAddSpinReconDir2Minus(const REAL32* src, REAL32 *dst, unsigned int n_vec) 
{

#ifdef DEBUG_GENERIC_SPIN_RECON_INLINES
  QDPIO::cout << "inlinaAddSpinReconDir0Minus" << endl;
#endif

  /*                              ( 1  0 -i  0)  ( a0 )    ( a0 - i a2 )
   *  B  :=  ( 1 - Gamma  ) A  =  ( 0  1  0  i)  ( a1 )  = ( a1 + i a3 )
   *                    2         ( i  0  1  0)  ( a2 )    ( a2 + i a0 )
   *                              ( 0 -i  0  1)  ( a3 )    ( a3 - i a1 )

   * The bottom components of be may be reconstructed using the formula
   *  ( b2r + i b2i )  =  ( {a2r - a0i} + i{a2i + a0r} )  =  ( - b0i + i b0r )
   *  ( b3r + i b3i )     ( {a3r + a1i} + i{a3i - a1r} )     (   b1i - i b1r )
   */

  const REAL32* src_shadow = src;
  REAL32* dst_shadow = dst;

  SSEVec v0, v1, v2, v3, v4, v5, v6, v7;
  v6.floats[0] = -1;
  v6.floats[1] = +1;
  v6.floats[2] = -1;
  v6.floats[3] = +1;

  v7.floats[0] = -1;
  v7.floats[1] = +1;
  v7.floats[2] = +1;
  v7.floats[3] = -1;

  // Load source 
  v0.vector = _mm_load_ps(src_shadow);
  v1.vector = _mm_load_ps(src_shadow+4);
  v2.vector = _mm_load_ps(src_shadow+8);
  
  // Load top half of result
  v3.vector = _mm_load_ps(dst_shadow);
  v4.vector = _mm_load_ps(dst_shadow+4);
  v5.vector = _mm_load_ps(dst_shadow+8);

  // Add source to result 
  v3.vector = _mm_add_ps(v3.vector, v0.vector);
  v4.vector = _mm_add_ps(v4.vector, v1.vector);
  v5.vector = _mm_add_ps(v5.vector, v2.vector);

  // Store out result
  _mm_store_ps(dst_shadow, v3.vector);
  _mm_store_ps(dst_shadow+4, v4.vector);
  _mm_store_ps(dst_shadow+8, v5.vector);

  for(unsigned int site=0; site < n_vec-1; site++) { 
    
    src_shadow += 12;
    _mm_prefetch((const char *)src_shadow, _MM_HINT_T0);

    // Now the bottom half.  -- preload v4 and v5 
    v3.vector = _mm_load_ps(dst_shadow+12);
    v4.vector = _mm_load_ps(dst_shadow+16);
    v5.vector = _mm_load_ps(dst_shadow+20);
    _mm_prefetch((const char *) dst_shadow+20, _MM_HINT_T0);

    // I want to set up 
    //
    //     v0[0] = v0[1] = tmp[0][col=0][im]
    //     v0[1] = v0[0] = tmp[0][col=0][re]
    //     v0[2] = v0[3] = tmp[0][col=1][im]
    //     v0[3] = v0[2] = tmp[0][col=1][re]
    // 
    // We can do this with v0 = shuf( v0, v0 ). Shuf code is 10 11 00 01 = xB1
    v0.vector = _mm_shuffle_ps(v0.vector, v0.vector, 0xB1);
    
    // Mask (+1,-1,+1,-1)
    v0.vector = _mm_mul_ps(v0.vector, v6.vector);
    v3.vector = _mm_add_ps(v3.vector, v0.vector);
    _mm_store_ps(dst_shadow+12, v3.vector);
    

    // I want to set up 
    //
    //     v1[0] = v1[1] = tmp[0][col=2][im]
    //     v1[1] = v1[0] = tmp[0][col=2][re]
    //     v1[2] = v1[3] = tmp[1][col=0][im]
    //     v1[3] = v1[2] = tmp[1][col=0][re]
    // 
    // We can do this with v1 = shuf( v1, v1 ). Shuf code is 10 11 00 01 = xB1
    v1.vector = _mm_shuffle_ps(v1.vector, v1.vector, 0xB1);

    // Need to multiply in mask (+1,-1,-1,+1) 
    v1.vector = _mm_mul_ps(v1.vector, v7.vector);
    v4.vector = _mm_add_ps(v4.vector, v1.vector);
    _mm_store_ps(dst_shadow+16, v4.vector);

    // I want to set up 
    //
    //     v2[0] = v2[1] = tmp[1][col=1][im]
    //     v2[1] = v2[0] = tmp[1][col=1][re]
    //     v2[2] = v2[2] = tmp[1][col=2][im]
    //     v2[3] = v2[3] = tmp[1][col=2][re]
    // 
    // We can do this with v2 = shuf( v2, v2 ). Shuf code is 10 11 00 01 = xB1
    v2.vector = _mm_shuffle_ps(v2.vector, v2.vector, 0xB1);

    // Mask in -1, +1, -1, +1 from v7
    v2.vector = _mm_mul_ps(v2.vector, v6.vector);
    v5.vector = _mm_sub_ps(v5.vector, v2.vector);
    _mm_store_ps(dst_shadow+20, v5.vector);

    dst_shadow+=24;

    // Load source 
    v0.vector = _mm_load_ps(src_shadow);
    v1.vector = _mm_load_ps(src_shadow+4);
    v2.vector = _mm_load_ps(src_shadow+8);
    
    // Load top half of result
    v3.vector = _mm_load_ps(dst_shadow);
    v4.vector = _mm_load_ps(dst_shadow+4);
    v5.vector = _mm_load_ps(dst_shadow+8);
    
    // Add source to result 
    v3.vector = _mm_add_ps(v3.vector, v0.vector);
    v4.vector = _mm_add_ps(v4.vector, v1.vector);
    v5.vector = _mm_add_ps(v5.vector, v2.vector);

    // Store out result
    _mm_store_ps(dst_shadow, v3.vector);
    _mm_store_ps(dst_shadow+4, v4.vector);
    _mm_store_ps(dst_shadow+8, v5.vector);
    
  }

  // Now the bottom half.  -- preload v4 and v5 
  // Now the bottom half.  -- preload v4 and v5 
  v3.vector = _mm_load_ps(dst_shadow+12);
  v4.vector = _mm_load_ps(dst_shadow+16);
  v5.vector = _mm_load_ps(dst_shadow+20);
  _mm_prefetch((const char *) dst_shadow+20, _MM_HINT_T0);

  // I want to set up 
  //
  //     v0[0] = v0[1] = tmp[0][col=0][im]
  //     v0[1] = v0[0] = tmp[0][col=0][re]
  //     v0[2] = v0[3] = tmp[0][col=1][im]
  //     v0[3] = v0[2] = tmp[0][col=1][re]
  // 
  // We can do this with v0 = shuf( v0, v0 ). Shuf code is 10 11 00 01 = xB1
  v0.vector = _mm_shuffle_ps(v0.vector, v0.vector, 0xB1);
  
  // Mask (+1,-1,+1,-1)
  v0.vector = _mm_mul_ps(v0.vector, v6.vector);
  v3.vector = _mm_add_ps(v3.vector, v0.vector);
  _mm_store_ps(dst_shadow+12, v3.vector);
  
  
  // I want to set up 
  //
  //     v1[0] = v1[1] = tmp[0][col=2][im]
  //     v1[1] = v1[0] = tmp[0][col=2][re]
  //     v1[2] = v1[3] = tmp[1][col=0][im]
  //     v1[3] = v1[2] = tmp[1][col=0][re]
  // 
  // We can do this with v1 = shuf( v1, v1 ). Shuf code is 10 11 00 01 = xB1
  v1.vector = _mm_shuffle_ps(v1.vector, v1.vector, 0xB1);
  
  // Need to multiply in mask (+1,-1,-1,+1) 
  v1.vector = _mm_mul_ps(v1.vector, v7.vector);
  v4.vector = _mm_add_ps(v4.vector, v1.vector);
  _mm_store_ps(dst_shadow+16, v4.vector);

  // I want to set up 
  //
  //     v2[0] = v2[1] = tmp[1][col=1][im]
  //     v2[1] = v2[0] = tmp[1][col=1][re]
  //     v2[2] = v2[2] = tmp[1][col=2][im]
  //     v2[3] = v2[3] = tmp[1][col=2][re]
  // 
  // We can do this with v2 = shuf( v2, v2 ). Shuf code is 10 11 00 01 = xB1
  v2.vector = _mm_shuffle_ps(v2.vector, v2.vector, 0xB1);
  
  // Mask in -1, +1, -1, +1 from v7
  v2.vector = _mm_mul_ps(v2.vector, v6.vector);
  v5.vector = _mm_sub_ps(v5.vector, v2.vector);
  _mm_store_ps(dst_shadow+20, v5.vector);

    

}


/** \brief Spin recon (1/2)(1+\gamma3)
 *
 * \param src (pointer to 2 vector source)
 * \param dst (pointer to 4 vector dest)
 * \param n_vec (number of vectors to reconstruct)

 * It is assumeed that src points to an array of  2 spin
 * 3 colour vectors with complex components. 
 * It is assumed that dst points to an array of 4 spin
 * 3 colour vectors with complex components. In both cases
 * ordering is that spin is slowest, and complex components
 * are fastest. */
inline  
void inlineAddSpinReconDir3Plus(const REAL32* src, REAL32 *dst, unsigned int n_vec) 
{

#ifdef DEBUG_GENERIC_SPIN_RECON_INLINES
  QDPIO::cout << "inlinaAddSpinReconDir0Plus" << endl;
#endif

  /*                              ( 1  0  1  0)  ( a0 )    ( a0 + a2 )
   *  B  :=  ( 1 + Gamma  ) A  =  ( 0  1  0  1)  ( a1 )  = ( a1 + a3 )
   *                    3         ( 1  0  1  0)  ( a2 )    ( a2 + a0 )
   *                              ( 0  1  0  1)  ( a3 )    ( a3 + a1 )
   
   * The bottom components of be may be reconstructed using the formula
   
   *   ( b2r + i b2i )  =  ( {a2r + a0r} + i{a2i + a0i} )  =  ( b0r + i b0i )
   *   ( b3r + i b3i )     ( {a3r + a1r} + i{a3i + a1i} )     ( b1r + i b1i ) 
   */


  const REAL32* src_shadow = src;
  REAL32* dst_shadow = dst;

  SSEVec v0, v1, v2, v3, v4, v5;

  // Load source 
  v0.vector = _mm_load_ps(src_shadow);
  v1.vector = _mm_load_ps(src_shadow+4);
  v2.vector = _mm_load_ps(src_shadow+8);
  
  // Load top half of result
  v3.vector = _mm_load_ps(dst_shadow);
  v4.vector = _mm_load_ps(dst_shadow+4);
  v5.vector = _mm_load_ps(dst_shadow+8);

  // Add source to result 
  v3.vector = _mm_add_ps(v3.vector, v0.vector);
  v4.vector = _mm_add_ps(v4.vector, v1.vector);
  v5.vector = _mm_add_ps(v5.vector, v2.vector);

  // Store out result
  _mm_store_ps(dst_shadow, v3.vector);
  _mm_store_ps(dst_shadow+4, v4.vector);
  _mm_store_ps(dst_shadow+8, v5.vector);

  for(unsigned int site=0; site < n_vec-1; site++) { 
    
    src_shadow += 12;
    _mm_prefetch((const char *)src_shadow, _MM_HINT_T0);


    // Now the bottom half.  -- preload v4 and v5 
    v3.vector = _mm_load_ps(dst_shadow+12);
    v4.vector = _mm_load_ps(dst_shadow+16);
    v5.vector = _mm_load_ps(dst_shadow+20);

    v3.vector = _mm_add_ps(v3.vector, v0.vector);
    v4.vector = _mm_add_ps(v4.vector, v1.vector);
    v5.vector = _mm_add_ps(v5.vector, v2.vector);

    _mm_store_ps(dst_shadow+12, v3.vector);
    _mm_store_ps(dst_shadow+16, v4.vector);
    _mm_store_ps(dst_shadow+20, v5.vector);

    dst_shadow+=24;

    // Load source 
    v0.vector = _mm_load_ps(src_shadow);
    v1.vector = _mm_load_ps(src_shadow+4);
    v2.vector = _mm_load_ps(src_shadow+8);
    
    // Load top half of result
    v3.vector = _mm_load_ps(dst_shadow);
    v4.vector = _mm_load_ps(dst_shadow+4);
    v5.vector = _mm_load_ps(dst_shadow+8);
    
    // Add source to result (no shufs needed, no mask needed)
    v3.vector = _mm_add_ps(v3.vector, v0.vector);
    v4.vector = _mm_add_ps(v4.vector, v1.vector);
    v5.vector = _mm_add_ps(v5.vector, v2.vector);

    // Store out result
    _mm_store_ps(dst_shadow, v3.vector);
    _mm_store_ps(dst_shadow+4, v4.vector);
    _mm_store_ps(dst_shadow+8, v5.vector);
    
  }

  // Now the bottom half.  -- preload v4 and v5 
  // Now the bottom half.  -- preload v4 and v5 
  v3.vector = _mm_load_ps(dst_shadow+12);
  v4.vector = _mm_load_ps(dst_shadow+16);
  v5.vector = _mm_load_ps(dst_shadow+20);
  
  v3.vector = _mm_add_ps(v3.vector, v0.vector);
  v4.vector = _mm_add_ps(v4.vector, v1.vector);
  v5.vector = _mm_add_ps(v5.vector, v2.vector);
  
  _mm_store_ps(dst_shadow+12, v3.vector);
  _mm_store_ps(dst_shadow+16, v4.vector);
  _mm_store_ps(dst_shadow+20, v5.vector);

}

/** \brief Spin recon (1/2)(1-\gamma3)
 *
 * \param src (pointer to 2 vector source)
 * \param dst (pointer to 4 vector dest)
 * \param n_vec (number of vectors to reconstruct)

 * It is assumeed that src points to an array of  2 spin
 * 3 colour vectors with complex components. 
 * It is assumed that dst points to an array of 4 spin
 * 3 colour vectors with complex components. In both cases
 * ordering is that spin is slowest, and complex components
 * are fastest. */
inline  
void inlineAddSpinReconDir3Minus(const REAL32* src, REAL32 *dst, unsigned int n_vec) 
{

#ifdef DEBUG_GENERIC_SPIN_RECON_INLINES
  QDPIO::cout << "inlinaAddSpinReconDir0Minus" << endl;
#endif

  /*                              ( 1  0 -1  0)  ( a0 )    ( a0 - a2 )
   *  B  :=  ( 1 - Gamma  ) A  =  ( 0  1  0 -1)  ( a1 )  = ( a1 - a3 )
   *                    3         (-1  0  1  0)  ( a2 )    ( a2 - a0 )
   *                              ( 0 -1  0  1)  ( a3 )    ( a3 - a1 )
   
   * The bottom components of be may be reconstructed using the formula
   *  ( b2r + i b2i )  =  ( {a2r - a0r} + i{a2i - a0i} )  =  ( - b0r - i b0i )
   *  ( b3r + i b3i )     ( {a3r - a1r} + i{a3i - a1i} )     ( - b1r - i b1i ) 
   */    
  
  const REAL32* src_shadow = src;
  REAL32* dst_shadow = dst;
  
  SSEVec v0, v1, v2, v3, v4, v5;
  
  // Load source 
  v0.vector = _mm_load_ps(src_shadow);
  v1.vector = _mm_load_ps(src_shadow+4);
  v2.vector = _mm_load_ps(src_shadow+8);
  
  // Load top half of result
  v3.vector = _mm_load_ps(dst_shadow);
  v4.vector = _mm_load_ps(dst_shadow+4);
  v5.vector = _mm_load_ps(dst_shadow+8);

  // Add source to result 
  v3.vector = _mm_add_ps(v3.vector, v0.vector);
  v4.vector = _mm_add_ps(v4.vector, v1.vector);
  v5.vector = _mm_add_ps(v5.vector, v2.vector);

  // Store out result
  _mm_store_ps(dst_shadow, v3.vector);
  _mm_store_ps(dst_shadow+4, v4.vector);
  _mm_store_ps(dst_shadow+8, v5.vector);

  for(unsigned int site=0; site < n_vec-1; site++) { 
    
    src_shadow += 12;
    _mm_prefetch((const char *)src_shadow, _MM_HINT_T0);

    // Now the bottom half.  -- preload v4 and v5 
    v3.vector = _mm_load_ps(dst_shadow+12);
    v4.vector = _mm_load_ps(dst_shadow+16);
    v5.vector = _mm_load_ps(dst_shadow+20);

    v3.vector = _mm_sub_ps(v3.vector, v0.vector);
    v4.vector = _mm_sub_ps(v4.vector, v1.vector);
    v5.vector = _mm_sub_ps(v5.vector, v2.vector);

    _mm_store_ps(dst_shadow+12, v3.vector);
    _mm_store_ps(dst_shadow+16, v4.vector);
    _mm_store_ps(dst_shadow+20, v5.vector);

    dst_shadow+=24;

    // Load source 
    v0.vector = _mm_load_ps(src_shadow);
    v1.vector = _mm_load_ps(src_shadow+4);
    v2.vector = _mm_load_ps(src_shadow+8);
    
    // Load top half of result
    v3.vector = _mm_load_ps(dst_shadow);
    v4.vector = _mm_load_ps(dst_shadow+4);
    v5.vector = _mm_load_ps(dst_shadow+8);
    
    // Add source to result (no shufs needed, no mask needed)
    v3.vector = _mm_add_ps(v3.vector, v0.vector);
    v4.vector = _mm_add_ps(v4.vector, v1.vector);
    v5.vector = _mm_add_ps(v5.vector, v2.vector);

    // Store out result
    _mm_store_ps(dst_shadow, v3.vector);
    _mm_store_ps(dst_shadow+4, v4.vector);
    _mm_store_ps(dst_shadow+8, v5.vector);
    
  }

  // Now the bottom half.  -- preload v4 and v5 
  // Now the bottom half.  -- preload v4 and v5 
  v3.vector = _mm_load_ps(dst_shadow+12);
  v4.vector = _mm_load_ps(dst_shadow+16);
  v5.vector = _mm_load_ps(dst_shadow+20);
  
  v3.vector = _mm_sub_ps(v3.vector, v0.vector);
  v4.vector = _mm_sub_ps(v4.vector, v1.vector);
  v5.vector = _mm_sub_ps(v5.vector, v2.vector);
  
  _mm_store_ps(dst_shadow+12, v3.vector);
  _mm_store_ps(dst_shadow+16, v4.vector);
  _mm_store_ps(dst_shadow+20, v5.vector);
  
}


} // namespace QDP

#endif 
