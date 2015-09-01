#ifndef GENERIC_SPIN_RECON_INLINES_H
#define GENERIC_SPIN_RECON_INLINES_H

/* File: generic_spin_recon_inlines.h
   Purpose: Supply inline functions to do spin reconstruction
   Author: $Id: generic_spin_recon_inlines.h,v 1.3 2007-07-17 16:56:09 bjoo Exp $
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
void inlineSpinReconDir0Plus(const REAL* src, REAL *dst, unsigned int n_vec) 
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

  REAL tmp_hspinor[4][3][2];

  const REAL* src_shadow = src;
  REAL* dst_shadow = dst;

  const int re = 0;
  const int im = 1;
  const int Ncmpx = 2;
  const int Nsby2 = 2;

  for(unsigned int site=0; site < n_vec; site++) {
    REAL* tmp_shadow = &(tmp_hspinor[0][0][0]);
    
    // Stream in the half spinor - write out the first two components
    for(int stream=0; stream < Nsby2*Nc*Ncmpx; stream++) {
      REAL tmp = *(src_shadow++);
      *(tmp_shadow++) = tmp;
      *(dst_shadow++) = tmp;
    }
    
    // Reconstruct Spin 2
    for(int col=0; col < Nc; col++) { 
      *(dst_shadow++) =  tmp_hspinor[1][col][im];
      *(dst_shadow++) = -tmp_hspinor[1][col][re];
    }

    // Reconstruct Spin 3;
    for(int col=0; col < Nc; col++) { 
      *(dst_shadow++) =  tmp_hspinor[0][col][im];
      *(dst_shadow++) = -tmp_hspinor[0][col][re];
    }


  }
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
void inlineSpinReconDir0Minus(const REAL* src, REAL *dst, unsigned int n_vec) 
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
  
  REAL tmp_hspinor[4][3][2];

  const REAL* src_shadow = src;
  REAL* dst_shadow = dst;

  const int re = 0;
  const int im = 1;
  const int Ncmpx = 2;
  const int Nsby2 = 2;

  for(unsigned int site=0; site < n_vec; site++) {
    REAL* tmp_shadow = &(tmp_hspinor[0][0][0]);
    
    // Stream in the half spinor - write out the first two components
    for(int stream=0; stream < Nsby2*Nc*Ncmpx; stream++) {
      REAL tmp = *(src_shadow++);
      *(tmp_shadow++) = tmp;
      *(dst_shadow++) = tmp;
    }
    
    // Reconstruct Spin 2
    for(int col=0; col < Nc; col++) { 
      *(dst_shadow++) = -tmp_hspinor[1][col][im];
      *(dst_shadow++) = tmp_hspinor[1][col][re];
    }

    // Reconstruct Spin 3;
    for(int col=0; col < Nc; col++) { 
      *(dst_shadow++) = -tmp_hspinor[0][col][im];
      *(dst_shadow++) = tmp_hspinor[0][col][re];
    }


  }
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
void inlineSpinReconDir1Plus(const REAL* src, REAL *dst, unsigned int n_vec) 
{

#ifdef DEBUG_GENERIC_SPIN_RECON_INLINES
  QDPIO::cout << "inlineSpinReconDir0Plus" << endl;
#endif



  REAL tmp_hspinor[4][3][2];

  const REAL* src_shadow = src;
  REAL* dst_shadow = dst;

  const int re = 0;
  const int im = 1;
  const int Ncmpx = 2;
  const int Nsby2 = 2;

  for(unsigned int site=0; site < n_vec; site++) {
    REAL* tmp_shadow = &(tmp_hspinor[0][0][0]);
    
    // Stream in the half spinor - write out the first two components
    for(int stream=0; stream < Nsby2*Nc*Ncmpx; stream++) {
      REAL tmp = *(src_shadow++);
      *(tmp_shadow++) = tmp;
      *(dst_shadow++) = tmp;
    }
    
    /* 1 + \gamma_1 =  1  0  0 -1 
     *                 0  1  1  0
     *                 0  1  1  0
     *                -1  0  0  1 
 
     *   ( b2r + i b2i )  =  ( {a2r + a1r} + i{a2i + a1i} )  =  (   b1r + i b1i )
     *   ( b3r + i b3i )     ( {a3r - a0r} + i{a3i - a0i} )     ( - b0r - i b0i ) 
  
    */
    // Reconstruct Spin 2
    for(int col=0; col < Nc; col++) { 
      *(dst_shadow++) = tmp_hspinor[1][col][re]; 
      *(dst_shadow++) = tmp_hspinor[1][col][im];
    }

    // Reconstruct Spin 3;
    for(int col=0; col < Nc; col++) { 
      *(dst_shadow++) = -tmp_hspinor[0][col][re];
      *(dst_shadow++) = -tmp_hspinor[0][col][im];
    }


  }
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
void inlineSpinReconDir1Minus(const REAL* src, REAL *dst, unsigned int n_vec) 
{

#ifdef DEBUG_GENERIC_SPIN_RECON_INLINES
  QDPIO::cout << "inlineSpinReconDir0Minus" << endl;
#endif

  
  REAL tmp_hspinor[4][3][2];

  const REAL* src_shadow = src;
  REAL* dst_shadow = dst;

  const int re = 0;
  const int im = 1;
  const int Ncmpx = 2;
  const int Nsby2 = 2;

  for(unsigned int site=0; site < n_vec; site++) {
    REAL* tmp_shadow = &(tmp_hspinor[0][0][0]);
    
    // Stream in the half spinor - write out the first two components
    for(int stream=0; stream < Nsby2*Nc*Ncmpx; stream++) {
      REAL tmp = *(src_shadow++);
      *(tmp_shadow++) = tmp;
      *(dst_shadow++) = tmp;
    }


    /*                              ( 1  0  0  1)  ( a0 )    ( a0 + a3 )
     *  B  :=  ( 1 - Gamma  ) A  =  ( 0  1 -1  0)  ( a1 )  = ( a1 - a2 )
     *                    1         ( 0 -1  1  0)  ( a2 )    ( a2 - a1 )
     *                              ( 1  0  0  1)  ( a3 )    ( a3 + a0 )
     
     * The bottom components of be may be reconstructed using the formula

     *  ( b2r + i b2i )  =  ( {a2r - a1r} + i{a2i - a1i} )  =  ( - b1r - i b1i )
     *  ( b3r + i b3i )     ( {a3r + a0r} + i{a3i + a0i} )     (   b0r + i b0i ) 
     */

    
    // Reconstruct Spin 2
    for(int col=0; col < Nc; col++) { 
      *(dst_shadow++) = -tmp_hspinor[1][col][re];
      *(dst_shadow++) = -tmp_hspinor[1][col][im];
    }

    // Reconstruct Spin 3;
    for(int col=0; col < Nc; col++) { 
      *(dst_shadow++) = tmp_hspinor[0][col][re];
      *(dst_shadow++) = tmp_hspinor[0][col][im];
    }


  }
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
void inlineSpinReconDir2Plus(const REAL* src, REAL *dst, unsigned int n_vec) 
{

#ifdef DEBUG_GENERIC_SPIN_RECON_INLINES
  QDPIO::cout << "inlineSpinReconDir0Plus" << endl;
#endif



  REAL tmp_hspinor[4][3][2];

  const REAL* src_shadow = src;
  REAL* dst_shadow = dst;

  const int re = 0;
  const int im = 1;
  const int Ncmpx = 2;
  const int Nsby2 = 2;

  for(unsigned int site=0; site < n_vec; site++) {
    REAL* tmp_shadow = &(tmp_hspinor[0][0][0]);
    
    // Stream in the half spinor - write out the first two components
    for(int stream=0; stream < Nsby2*Nc*Ncmpx; stream++) {
      REAL tmp = *(src_shadow++);
      *(tmp_shadow++) = tmp;
      *(dst_shadow++) = tmp;
    }

  /* 1 + \gamma_2 =  1  0  i  0 
   *                 0  1  0 -i
   *                -i  0  1  0
   *                 0  i  0  1 
   *	     
   *  ( b2r + i b2i )  =  ( {a2r + a0i} + i{a2i - a0r} )  =  (   b0i - i b0r )
   *  ( b3r + i b3i )     ( {a3r - a1i} + i{a3i + a1r} )     ( - b1i + i b1r ) 
  */    

    // Reconstruct Spin 2
    for(int col=0; col < Nc; col++) { 
      *(dst_shadow++) =  tmp_hspinor[0][col][im]; 
      *(dst_shadow++) = -tmp_hspinor[0][col][re];
    }

    // Reconstruct Spin 3;
    for(int col=0; col < Nc; col++) { 
      *(dst_shadow++) = -tmp_hspinor[1][col][im];
      *(dst_shadow++) =  tmp_hspinor[1][col][re];
    }


  }
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
void inlineSpinReconDir2Minus(const REAL* src, REAL *dst, unsigned int n_vec) 
{

#ifdef DEBUG_GENERIC_SPIN_RECON_INLINES
  QDPIO::cout << "inlineSpinReconDir0Minus" << endl;
#endif

  
  REAL tmp_hspinor[4][3][2];

  const REAL* src_shadow = src;
  REAL* dst_shadow = dst;

  const int re = 0;
  const int im = 1;
  const int Ncmpx = 2;
  const int Nsby2 = 2;

  for(unsigned int site=0; site < n_vec; site++) {
    REAL* tmp_shadow = &(tmp_hspinor[0][0][0]);
    
    // Stream in the half spinor - write out the first two components
    for(int stream=0; stream < Nsby2*Nc*Ncmpx; stream++) {
      REAL tmp = *(src_shadow++);
      *(tmp_shadow++) = tmp;
      *(dst_shadow++) = tmp;
    }

  /*                              ( 1  0 -i  0)  ( a0 )    ( a0 - i a2 )
   *  B  :=  ( 1 - Gamma  ) A  =  ( 0  1  0  i)  ( a1 )  = ( a1 + i a3 )
   *                    2         ( i  0  1  0)  ( a2 )    ( a2 + i a0 )
   *                              ( 0 -i  0  1)  ( a3 )    ( a3 - i a1 )

   * The bottom components of be may be reconstructed using the formula
   *  ( b2r + i b2i )  =  ( {a2r - a0i} + i{a2i + a0r} )  =  ( - b0i + i b0r )
   *  ( b3r + i b3i )     ( {a3r + a1i} + i{a3i - a1r} )     (   b1i - i b1r )
   */
    
    // Reconstruct Spin 2
    for(int col=0; col < Nc; col++) { 
      *(dst_shadow++) = -tmp_hspinor[0][col][im];
      *(dst_shadow++) =  tmp_hspinor[0][col][re];
    }

    // Reconstruct Spin 3;
    for(int col=0; col < Nc; col++) { 
      *(dst_shadow++) =  tmp_hspinor[1][col][im];
      *(dst_shadow++) = -tmp_hspinor[1][col][re];
    }


  }
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
void inlineSpinReconDir3Plus(const REAL* src, REAL *dst, unsigned int n_vec) 
{

#ifdef DEBUG_GENERIC_SPIN_RECON_INLINES
  QDPIO::cout << "inlineSpinReconDir0Plus" << endl;
#endif



  REAL tmp_hspinor[4][3][2];

  const REAL* src_shadow = src;
  REAL* dst_shadow = dst;

  const int re = 0;
  const int im = 1;
  const int Ncmpx = 2;
  const int Nsby2 = 2;

  for(unsigned int site=0; site < n_vec; site++) {
    REAL* tmp_shadow = &(tmp_hspinor[0][0][0]);
    
    // Stream in the half spinor - write out the first two components
    for(int stream=0; stream < Nsby2*Nc*Ncmpx; stream++) {
      REAL tmp = *(src_shadow++);
      *(tmp_shadow++) = tmp;
      *(dst_shadow++) = tmp;
    }

  /*                              ( 1  0  1  0)  ( a0 )    ( a0 + a2 )
   *  B  :=  ( 1 + Gamma  ) A  =  ( 0  1  0  1)  ( a1 )  = ( a1 + a3 )
   *                    3         ( 1  0  1  0)  ( a2 )    ( a2 + a0 )
   *                              ( 0  1  0  1)  ( a3 )    ( a3 + a1 )
   
   * The bottom components of be may be reconstructed using the formula
   
   *   ( b2r + i b2i )  =  ( {a2r + a0r} + i{a2i + a0i} )  =  ( b0r + i b0i )
   *   ( b3r + i b3i )     ( {a3r + a1r} + i{a3i + a1i} )     ( b1r + i b1i ) 
   */
    // Reconstruct Spin 2
    for(int col=0; col < Nc; col++) { 
      *(dst_shadow++) =  tmp_hspinor[0][col][re]; 
      *(dst_shadow++) =  tmp_hspinor[0][col][im];
    }

    // Reconstruct Spin 3;
    for(int col=0; col < Nc; col++) { 
      *(dst_shadow++) =  tmp_hspinor[1][col][re];
      *(dst_shadow++) =  tmp_hspinor[1][col][im];
    }


  }
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
void inlineSpinReconDir3Minus(const REAL* src, REAL *dst, unsigned int n_vec) 
{

#ifdef DEBUG_GENERIC_SPIN_RECON_INLINES
  QDPIO::cout << "inlineSpinReconDir0Minus" << endl;
#endif

  
  REAL tmp_hspinor[4][3][2];

  const REAL* src_shadow = src;
  REAL* dst_shadow = dst;

  const int re = 0;
  const int im = 1;
  const int Ncmpx = 2;
  const int Nsby2 = 2;

  for(unsigned int site=0; site < n_vec; site++) {
    REAL* tmp_shadow = &(tmp_hspinor[0][0][0]);
    
    // Stream in the half spinor - write out the first two components
    for(int stream=0; stream < Nsby2*Nc*Ncmpx; stream++) {
      REAL tmp = *(src_shadow++);
      *(tmp_shadow++) = tmp;
      *(dst_shadow++) = tmp;
    }

    /*                              ( 1  0 -1  0)  ( a0 )    ( a0 - a2 )
     *  B  :=  ( 1 - Gamma  ) A  =  ( 0  1  0 -1)  ( a1 )  = ( a1 - a3 )
     *                    3         (-1  0  1  0)  ( a2 )    ( a2 - a0 )
     *                              ( 0 -1  0  1)  ( a3 )    ( a3 - a1 )
     
     * The bottom components of be may be reconstructed using the formula
     *  ( b2r + i b2i )  =  ( {a2r - a0r} + i{a2i - a0i} )  =  ( - b0r - i b0i )
     *  ( b3r + i b3i )     ( {a3r - a1r} + i{a3i - a1i} )     ( - b1r - i b1i ) 
     */    
    // Reconstruct Spin 2
    for(int col=0; col < Nc; col++) { 
      *(dst_shadow++) = -tmp_hspinor[0][col][re];
      *(dst_shadow++) = -tmp_hspinor[0][col][im];
    }

    // Reconstruct Spin 3;
    for(int col=0; col < Nc; col++) { 
      *(dst_shadow++) = -tmp_hspinor[1][col][re];
      *(dst_shadow++) = -tmp_hspinor[1][col][im];
    }


  }
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
void inlineAddSpinReconDir0Plus(const REAL* src, REAL *dst, unsigned int n_vec) 
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

  REAL tmp_hspinor[4][3][2];

  const REAL* src_shadow = src;
  REAL* dst_shadow = dst;

  const int re = 0;
  const int im = 1;
  const int Ncmpx = 2;
  const int Nsby2 = 2;

  for(unsigned int site=0; site < n_vec; site++) {
    REAL* tmp_shadow = &(tmp_hspinor[0][0][0]);
    
    // Stream in the half spinor - write out the first two components
    for(int stream=0; stream < Nsby2*Nc*Ncmpx; stream++) {
      REAL tmp = *(src_shadow++);
      *(tmp_shadow++) = tmp;
      *(dst_shadow++) += tmp;
    }
    
    // Reconstruct Spin 2
    for(int col=0; col < Nc; col++) { 
      *(dst_shadow++) +=  tmp_hspinor[1][col][im];
      *(dst_shadow++) -=  tmp_hspinor[1][col][re];
    }

    // Reconstruct Spin 3;
    for(int col=0; col < Nc; col++) { 
      *(dst_shadow++) +=  tmp_hspinor[0][col][im];
      *(dst_shadow++) -=  tmp_hspinor[0][col][re];
    }


  }
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
void inlineAddSpinReconDir0Minus(const REAL* src, REAL *dst, unsigned int n_vec) 
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
  
  REAL tmp_hspinor[4][3][2];

  const REAL* src_shadow = src;
  REAL* dst_shadow = dst;

  const int re = 0;
  const int im = 1;
  const int Ncmpx = 2;
  const int Nsby2 = 2;

  for(unsigned int site=0; site < n_vec; site++) {
    REAL* tmp_shadow = &(tmp_hspinor[0][0][0]);
    
    // Stream in the half spinor - write out the first two components
    for(int stream=0; stream < Nsby2*Nc*Ncmpx; stream++) {
      REAL tmp = *(src_shadow++);
      *(tmp_shadow++) = tmp;
      *(dst_shadow++) += tmp;
    }
    
    // Reconstruct Spin 2
    for(int col=0; col < Nc; col++) { 
      *(dst_shadow++) -= tmp_hspinor[1][col][im];
      *(dst_shadow++) += tmp_hspinor[1][col][re];
    }

    // Reconstruct Spin 3;
    for(int col=0; col < Nc; col++) { 
      *(dst_shadow++) -= tmp_hspinor[0][col][im];
      *(dst_shadow++) += tmp_hspinor[0][col][re];
    }


  }
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
void inlineAddSpinReconDir1Plus(const REAL* src, REAL *dst, unsigned int n_vec) 
{

#ifdef DEBUG_GENERIC_SPIN_RECON_INLINES
  QDPIO::cout << "inlinaAddSpinReconDir0Plus" << endl;
#endif



  REAL tmp_hspinor[4][3][2];

  const REAL* src_shadow = src;
  REAL* dst_shadow = dst;

  const int re = 0;
  const int im = 1;
  const int Ncmpx = 2;
  const int Nsby2 = 2;

  for(unsigned int site=0; site < n_vec; site++) {
    REAL* tmp_shadow = &(tmp_hspinor[0][0][0]);
    
    // Stream in the half spinor - write out the first two components
    for(int stream=0; stream < Nsby2*Nc*Ncmpx; stream++) {
      REAL tmp = *(src_shadow++);
      *(tmp_shadow++) = tmp;
      *(dst_shadow++) += tmp;
    }
    
    /* 1 + \gamma_1 =  1  0  0 -1 
     *                 0  1  1  0
     *                 0  1  1  0
     *                -1  0  0  1 
 
     *   ( b2r + i b2i )  =  ( {a2r + a1r} + i{a2i + a1i} )  =  (   b1r + i b1i )
     *   ( b3r + i b3i )     ( {a3r - a0r} + i{a3i - a0i} )     ( - b0r - i b0i ) 
  
    */
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
void inlineAddSpinReconDir1Minus(const REAL* src, REAL *dst, unsigned int n_vec) 
{

#ifdef DEBUG_GENERIC_SPIN_RECON_INLINES
  QDPIO::cout << "inlinaAddSpinReconDir0Minus" << endl;
#endif

  
  REAL tmp_hspinor[4][3][2];

  const REAL* src_shadow = src;
  REAL* dst_shadow = dst;

  const int re = 0;
  const int im = 1;
  const int Ncmpx = 2;
  const int Nsby2 = 2;

  for(unsigned int site=0; site < n_vec; site++) {
    REAL* tmp_shadow = &(tmp_hspinor[0][0][0]);
    
    // Stream in the half spinor - write out the first two components
    for(int stream=0; stream < Nsby2*Nc*Ncmpx; stream++) {
      REAL tmp = *(src_shadow++);
      *(tmp_shadow++) = tmp;
      *(dst_shadow++) += tmp;
    }


    /*                              ( 1  0  0  1)  ( a0 )    ( a0 + a3 )
     *  B  :=  ( 1 - Gamma  ) A  =  ( 0  1 -1  0)  ( a1 )  = ( a1 - a2 )
     *                    1         ( 0 -1  1  0)  ( a2 )    ( a2 - a1 )
     *                              ( 1  0  0  1)  ( a3 )    ( a3 + a0 )
     
     * The bottom components of be may be reconstructed using the formula

     *  ( b2r + i b2i )  =  ( {a2r - a1r} + i{a2i - a1i} )  =  ( - b1r - i b1i )
     *  ( b3r + i b3i )     ( {a3r + a0r} + i{a3i + a0i} )     (   b0r + i b0i ) 
     */

    
    // Reconstruct Spin 2
    for(int col=0; col < Nc; col++) { 
      *(dst_shadow++) -= tmp_hspinor[1][col][re];
      *(dst_shadow++) -= tmp_hspinor[1][col][im];
    }

    // Reconstruct Spin 3;
    for(int col=0; col < Nc; col++) { 
      *(dst_shadow++) += tmp_hspinor[0][col][re];
      *(dst_shadow++) += tmp_hspinor[0][col][im];
    }


  }
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
void inlineAddSpinReconDir2Plus(const REAL* src, REAL *dst, unsigned int n_vec) 
{

#ifdef DEBUG_GENERIC_SPIN_RECON_INLINES
  QDPIO::cout << "inlinaAddSpinReconDir0Plus" << endl;
#endif



  REAL tmp_hspinor[4][3][2];

  const REAL* src_shadow = src;
  REAL* dst_shadow = dst;

  const int re = 0;
  const int im = 1;
  const int Ncmpx = 2;
  const int Nsby2 = 2;

  for(unsigned int site=0; site < n_vec; site++) {
    REAL* tmp_shadow = &(tmp_hspinor[0][0][0]);
    
    // Stream in the half spinor - write out the first two components
    for(int stream=0; stream < Nsby2*Nc*Ncmpx; stream++) {
      REAL tmp = *(src_shadow++);
      *(tmp_shadow++) = tmp;
      *(dst_shadow++) += tmp;
    }

  /* 1 + \gamma_2 =  1  0  i  0 
   *                 0  1  0 -i
   *                -i  0  1  0
   *                 0  i  0  1 
   *	     
   *  ( b2r + i b2i )  =  ( {a2r + a0i} + i{a2i - a0r} )  =  (   b0i - i b0r )
   *  ( b3r + i b3i )     ( {a3r - a1i} + i{a3i + a1r} )     ( - b1i + i b1r ) 
  */    

    // Reconstruct Spin 2
    for(int col=0; col < Nc; col++) { 
      *(dst_shadow++) +=  tmp_hspinor[0][col][im]; 
      *(dst_shadow++) -=  tmp_hspinor[0][col][re];
    }

    // Reconstruct Spin 3;
    for(int col=0; col < Nc; col++) { 
      *(dst_shadow++) -= tmp_hspinor[1][col][im];
      *(dst_shadow++) += tmp_hspinor[1][col][re];
    }


  }
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
void inlineAddSpinReconDir2Minus(const REAL* src, REAL *dst, unsigned int n_vec) 
{

#ifdef DEBUG_GENERIC_SPIN_RECON_INLINES
  QDPIO::cout << "inlinaAddSpinReconDir0Minus" << endl;
#endif

  
  REAL tmp_hspinor[4][3][2];

  const REAL* src_shadow = src;
  REAL* dst_shadow = dst;

  const int re = 0;
  const int im = 1;
  const int Ncmpx = 2;
  const int Nsby2 = 2;

  for(unsigned int site=0; site < n_vec; site++) {
    REAL* tmp_shadow = &(tmp_hspinor[0][0][0]);
    
    // Stream in the half spinor - write out the first two components
    for(int stream=0; stream < Nsby2*Nc*Ncmpx; stream++) {
      REAL tmp = *(src_shadow++);
      *(tmp_shadow++) = tmp;
      *(dst_shadow++) += tmp;
    }

  /*                              ( 1  0 -i  0)  ( a0 )    ( a0 - i a2 )
   *  B  :=  ( 1 - Gamma  ) A  =  ( 0  1  0  i)  ( a1 )  = ( a1 + i a3 )
   *                    2         ( i  0  1  0)  ( a2 )    ( a2 + i a0 )
   *                              ( 0 -i  0  1)  ( a3 )    ( a3 - i a1 )

   * The bottom components of be may be reconstructed using the formula
   *  ( b2r + i b2i )  =  ( {a2r - a0i} + i{a2i + a0r} )  =  ( - b0i + i b0r )
   *  ( b3r + i b3i )     ( {a3r + a1i} + i{a3i - a1r} )     (   b1i - i b1r )
   */
    
    // Reconstruct Spin 2
    for(int col=0; col < Nc; col++) { 
      *(dst_shadow++) -= tmp_hspinor[0][col][im];
      *(dst_shadow++) +=  tmp_hspinor[0][col][re];
    }

    // Reconstruct Spin 3;
    for(int col=0; col < Nc; col++) { 
      *(dst_shadow++) +=  tmp_hspinor[1][col][im];
      *(dst_shadow++) -=  tmp_hspinor[1][col][re];
    }


  }
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
void inlineAddSpinReconDir3Plus(const REAL* src, REAL *dst, unsigned int n_vec) 
{

#ifdef DEBUG_GENERIC_SPIN_RECON_INLINES
  QDPIO::cout << "inlinaAddSpinReconDir0Plus" << endl;
#endif



  REAL tmp_hspinor[4][3][2];

  const REAL* src_shadow = src;
  REAL* dst_shadow = dst;

  const int re = 0;
  const int im = 1;
  const int Ncmpx = 2;
  const int Nsby2 = 2;

  for(unsigned int site=0; site < n_vec; site++) {
    REAL* tmp_shadow = &(tmp_hspinor[0][0][0]);
    
    // Stream in the half spinor - write out the first two components
    for(int stream=0; stream < Nsby2*Nc*Ncmpx; stream++) {
      REAL tmp = *(src_shadow++);
      *(tmp_shadow++) = tmp;
      *(dst_shadow++) += tmp;
    }

  /*                              ( 1  0  1  0)  ( a0 )    ( a0 + a2 )
   *  B  :=  ( 1 + Gamma  ) A  =  ( 0  1  0  1)  ( a1 )  = ( a1 + a3 )
   *                    3         ( 1  0  1  0)  ( a2 )    ( a2 + a0 )
   *                              ( 0  1  0  1)  ( a3 )    ( a3 + a1 )
   
   * The bottom components of be may be reconstructed using the formula
   
   *   ( b2r + i b2i )  =  ( {a2r + a0r} + i{a2i + a0i} )  =  ( b0r + i b0i )
   *   ( b3r + i b3i )     ( {a3r + a1r} + i{a3i + a1i} )     ( b1r + i b1i ) 
   */
    // Reconstruct Spin 2
    for(int col=0; col < Nc; col++) { 
      *(dst_shadow++) +=  tmp_hspinor[0][col][re]; 
      *(dst_shadow++) +=  tmp_hspinor[0][col][im];
    }

    // Reconstruct Spin 3;
    for(int col=0; col < Nc; col++) { 
      *(dst_shadow++) +=  tmp_hspinor[1][col][re];
      *(dst_shadow++) +=  tmp_hspinor[1][col][im];
    }


  }
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
void inlineAddSpinReconDir3Minus(const REAL* src, REAL *dst, unsigned int n_vec) 
{

#ifdef DEBUG_GENERIC_SPIN_RECON_INLINES
  QDPIO::cout << "inlinaAddSpinReconDir0Minus" << endl;
#endif

  
  REAL tmp_hspinor[4][3][2];

  const REAL* src_shadow = src;
  REAL* dst_shadow = dst;

  const int re = 0;
  const int im = 1;
  const int Ncmpx = 2;
  const int Nsby2 = 2;

  for(unsigned int site=0; site < n_vec; site++) {
    REAL* tmp_shadow = &(tmp_hspinor[0][0][0]);
    
    // Stream in the half spinor - write out the first two components
    for(int stream=0; stream < Nsby2*Nc*Ncmpx; stream++) {
      REAL tmp = *(src_shadow++);
      *(tmp_shadow++) = tmp;
      *(dst_shadow++) += tmp;
    }

    /*                              ( 1  0 -1  0)  ( a0 )    ( a0 - a2 )
     *  B  :=  ( 1 - Gamma  ) A  =  ( 0  1  0 -1)  ( a1 )  = ( a1 - a3 )
     *                    3         (-1  0  1  0)  ( a2 )    ( a2 - a0 )
     *                              ( 0 -1  0  1)  ( a3 )    ( a3 - a1 )
     
     * The bottom components of be may be reconstructed using the formula
     *  ( b2r + i b2i )  =  ( {a2r - a0r} + i{a2i - a0i} )  =  ( - b0r - i b0i )
     *  ( b3r + i b3i )     ( {a3r - a1r} + i{a3i - a1i} )     ( - b1r - i b1i ) 
     */    
    // Reconstruct Spin 2
    for(int col=0; col < Nc; col++) { 
      *(dst_shadow++) -= tmp_hspinor[0][col][re];
      *(dst_shadow++) -= tmp_hspinor[0][col][im];
    }

    // Reconstruct Spin 3;
    for(int col=0; col < Nc; col++) { 
      *(dst_shadow++) -= tmp_hspinor[1][col][re];
      *(dst_shadow++) -= tmp_hspinor[1][col][im];
    }


  }
}


} // namespace QDP

#endif 
