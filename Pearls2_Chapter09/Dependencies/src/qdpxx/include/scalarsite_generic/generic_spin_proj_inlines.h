#ifndef GENERIC_SPIN_PROJ_INLINES_H
#define GENERIC_SPIN_PROJ_INLINES_H

/* File: generic_spin_proj_inlines.h
   Purpose: Supply inline functions to do spin projection
   Author: $Id: generic_spin_proj_inlines.h,v 1.4 2008-12-22 17:42:57 bjoo Exp $
*/
namespace QDP {


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
void inlineSpinProjDir0Plus(const REAL* src, REAL *dst, unsigned int n_vec) 
{

#ifdef DEBUG_GENERIC_SPIN_PROJ_INLINES
  QDPIO::cout << "inlineSpinProjDir0Plus" << endl;
#endif

  /* 1 + \gamma_0 =  1  0  0  i 
                     0  1  i  0
                     0 -i  1  0
                    -i  0  0  1 
 
   *      ( d0r + i d0i )  =  ( {x0r - x3i} + i{x0i + x3r} )
   *      ( d1r + i d1i )     ( {x1r - x2i} + i{x1i + x2r} )
   */
  REAL tmp_spinor[4][3][2];

  const REAL* src_shadow = src;
  REAL* dst_shadow = dst;

  const int re = 0;
  const int im = 1;
  const int Ncmpx = 2;
 
  for(unsigned int site=0; site < n_vec; site++) {
    REAL* tmp_shadow = &(tmp_spinor[0][0][0]);    

    // Stream in the spinor
    //ut << "site = "<< site << ", Nc = " << Nc <<", Ns * Nc * Ncmpx = " << Ns*Nc*Ncmpx << endl;

    for(int stream=0; stream < Ns*Nc*Ncmpx; stream++) {
      *(tmp_shadow++) = *(src_shadow++);
    }

     // Project and store
    // Spin 0
    for(int col=0; col < Nc; col++) { 
      *(dst_shadow++) = tmp_spinor[0][col][re] - tmp_spinor[3][col][im];
      *(dst_shadow++) = tmp_spinor[0][col][im] + tmp_spinor[3][col][re];
    }
    
    // Spin 1
    for(int col=0; col < Nc; col++) { 
      *(dst_shadow++) = tmp_spinor[1][col][re] - tmp_spinor[2][col][im];
      *(dst_shadow++) = tmp_spinor[1][col][im] + tmp_spinor[2][col][re];
    }
  }
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
void inlineSpinProjDir0Minus(const REAL* src, REAL *dst, unsigned int n_vec) 
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
  REAL tmp_spinor[4][3][2];

  const REAL* src_shadow = src;
  REAL* dst_shadow = dst;

  const int re = 0;
  const int im = 1;
  const int Ncmpx = 2;
  
  for(unsigned int site=0; site < n_vec; site++) {
    REAL* tmpptr = &(tmp_spinor[0][0][0]);        
    // Stream in the spinor
    for(int stream=0; stream < Ns*Nc*Ncmpx; stream++) {
      *(tmpptr++) = *(src_shadow++);
    }
    
    // Project and store
    // Spin 0
    for(int col=0; col < Nc; col++) { 
      *(dst_shadow++) = tmp_spinor[0][col][re] + tmp_spinor[3][col][im];
      *(dst_shadow++) = tmp_spinor[0][col][im] - tmp_spinor[3][col][re];
    }
    
    // Spin 1
    for(int col=0; col < Nc; col++) { 
      *(dst_shadow++) = tmp_spinor[1][col][re] + tmp_spinor[2][col][im];
      *(dst_shadow++) = tmp_spinor[1][col][im] - tmp_spinor[2][col][re];
    }
  }
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
void inlineSpinProjDir1Plus(const REAL* src, REAL *dst, unsigned int n_vec) 
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
  REAL tmp_spinor[4][3][2];

  const REAL* src_shadow = src;
  REAL* dst_shadow = dst;

  const int re = 0;
  const int im = 1;
  const int Ncmpx = 2;

  for(unsigned int site=0; site < n_vec; site++) {
    REAL* tmp_shadow = &(tmp_spinor[0][0][0]);
    
    // Stream in the spinor
    for(int stream=0; stream < Ns*Nc*Ncmpx; stream++) {
      *(tmp_shadow++) = *(src_shadow++);
    }
    
    // Project and store
    // Spin 0
    for(int col=0; col < Nc; col++) { 
      *(dst_shadow++) = tmp_spinor[0][col][re] - tmp_spinor[3][col][re];
      *(dst_shadow++) = tmp_spinor[0][col][im] - tmp_spinor[3][col][im];
    }
    
    // Spin 1
    for(int col=0; col < Nc; col++) { 
      *(dst_shadow++) = tmp_spinor[1][col][re] + tmp_spinor[2][col][re];
      *(dst_shadow++) = tmp_spinor[1][col][im] + tmp_spinor[2][col][im];
    }
  }
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
void inlineSpinProjDir1Minus(const REAL* src, REAL *dst, unsigned int n_vec) 
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

  REAL tmp_spinor[4][3][2];

  const REAL* src_shadow = src;
  REAL* dst_shadow = dst;

  const int re = 0;
  const int im = 1;
  const int Ncmpx = 2;
  
  for(unsigned int site=0; site < n_vec; site++) {
    REAL* tmpptr = &(tmp_spinor[0][0][0]);
    
    // Stream in the spinor
    for(int stream=0; stream < Ns*Nc*Ncmpx; stream++) {
      *(tmpptr++) = *(src_shadow++);
    }
    
    // Project and store
    // Spin 0
    for(int col=0; col < Nc; col++) { 
      *(dst_shadow++) = tmp_spinor[0][col][re] + tmp_spinor[3][col][re];
      *(dst_shadow++) = tmp_spinor[0][col][im] + tmp_spinor[3][col][im];
    }
    
    // Spin 1
    for(int col=0; col < Nc; col++) { 
      *(dst_shadow++) = tmp_spinor[1][col][re] - tmp_spinor[2][col][re];
      *(dst_shadow++) = tmp_spinor[1][col][im] - tmp_spinor[2][col][im];
    }
  }
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
void inlineSpinProjDir2Plus(const REAL* src, REAL *dst, unsigned int n_vec) 
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

  REAL tmp_spinor[4][3][2];

  const REAL* src_shadow = src;
  REAL* dst_shadow = dst;

  const int re = 0;
  const int im = 1;
  const int Ncmpx = 2;

  for(unsigned int site=0; site < n_vec; site++) {
    REAL* tmp_shadow = &(tmp_spinor[0][0][0]);
    
    // Stream in the spinor
    for(int stream=0; stream < Ns*Nc*Ncmpx; stream++) {
      *(tmp_shadow++) = *(src_shadow++);
    }
    
    // Project and store
    // Spin 0
    for(int col=0; col < Nc; col++) { 
      *(dst_shadow++) = tmp_spinor[0][col][re] - tmp_spinor[2][col][im];
      *(dst_shadow++) = tmp_spinor[0][col][im] + tmp_spinor[2][col][re];
    }
    
    // Spin 1
    for(int col=0; col < Nc; col++) { 
      *(dst_shadow++) = tmp_spinor[1][col][re] + tmp_spinor[3][col][im];
      *(dst_shadow++) = tmp_spinor[1][col][im] - tmp_spinor[3][col][re];
    }
  }
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
void inlineSpinProjDir2Minus(const REAL* src, REAL *dst, unsigned int n_vec) 
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

  REAL tmp_spinor[4][3][2];

  const REAL* src_shadow = src;
  REAL* dst_shadow = dst;

  const int re = 0;
  const int im = 1;
  const int Ncmpx = 2;
  
  for(unsigned int site=0; site < n_vec; site++) {
    REAL* tmpptr = &(tmp_spinor[0][0][0]);
    
    // Stream in the spinor
    for(int stream=0; stream < Ns*Nc*Ncmpx; stream++) {
      *(tmpptr++) = *(src_shadow++);
    }
    
    // Project and store
    // Spin 0
    for(int col=0; col < Nc; col++) { 
      *(dst_shadow++) = tmp_spinor[0][col][re] + tmp_spinor[2][col][im];
      *(dst_shadow++) = tmp_spinor[0][col][im] - tmp_spinor[2][col][re];
    }
    
    // Spin 1
    for(int col=0; col < Nc; col++) { 
      *(dst_shadow++) = tmp_spinor[1][col][re] - tmp_spinor[3][col][im];
      *(dst_shadow++) = tmp_spinor[1][col][im] + tmp_spinor[3][col][re];
    }
  }
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
void inlineSpinProjDir3Plus(const REAL* src, REAL *dst, unsigned int n_vec) 
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

  REAL tmp_spinor[4][3][2];

  const REAL* src_shadow = src;
  REAL* dst_shadow = dst;

  const int re = 0;
  const int im = 1;
  const int Ncmpx = 2;

  for(unsigned int site=0; site < n_vec; site++) {
    REAL* tmp_shadow = &(tmp_spinor[0][0][0]);
    
    // Stream in the spinor
    for(int stream=0; stream < Ns*Nc*Ncmpx; stream++) {
      *(tmp_shadow++) = *(src_shadow++);
    }
    
    // Project and store
    // Spin 0
    for(int col=0; col < Nc; col++) { 
      *(dst_shadow++) = tmp_spinor[0][col][re] + tmp_spinor[2][col][re];
      *(dst_shadow++) = tmp_spinor[0][col][im] + tmp_spinor[2][col][im];
    }
    
    // Spin 1
    for(int col=0; col < Nc; col++) { 
      *(dst_shadow++) = tmp_spinor[1][col][re] + tmp_spinor[3][col][re];
      *(dst_shadow++) = tmp_spinor[1][col][im] + tmp_spinor[3][col][im];
    }
  }
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
void inlineSpinProjDir3Minus(const REAL* src, REAL *dst, unsigned int n_vec) 
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


  REAL tmp_spinor[4][3][2];

  const REAL* src_shadow = src;
  REAL* dst_shadow = dst;

  const int re = 0;
  const int im = 1;
  const int Ncmpx = 2;
  
  for(unsigned int site=0; site < n_vec; site++) {
    REAL* tmpptr = &(tmp_spinor[0][0][0]);
    
    // Stream in the spinor
    for(int stream=0; stream < Ns*Nc*Ncmpx; stream++) {
      *(tmpptr++) = *(src_shadow++);
    }
    
    // Project and store
    // Spin 0
    for(int col=0; col < Nc; col++) { 
      *(dst_shadow++) = tmp_spinor[0][col][re] - tmp_spinor[2][col][re];
      *(dst_shadow++) = tmp_spinor[0][col][im] - tmp_spinor[2][col][im];
    }
    
    // Spin 1
    for(int col=0; col < Nc; col++) { 
      *(dst_shadow++) = tmp_spinor[1][col][re] - tmp_spinor[3][col][re];
      *(dst_shadow++) = tmp_spinor[1][col][im] - tmp_spinor[3][col][im];
    }
  }
}

} // namespace QDP;

#endif 
