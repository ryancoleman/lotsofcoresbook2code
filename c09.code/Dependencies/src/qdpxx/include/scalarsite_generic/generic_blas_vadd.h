// $Id: generic_blas_vadd.h,v 1.2 2007-06-10 14:32:10 edwards Exp $

/*! @file
 *  @brief Generic Scalar VADD routine
 *
 */

#ifndef QDP_GENERIC_BLAS_VADD
#define QDP_GENERIC_BLAS_VADD

namespace QDP {

// (Vector) Out = (Vector) In1 + (Vector) In2
inline
void vadd(REAL *Out, REAL *In1, REAL *In2, int n_3vec)
{
  
  int len = n_3vec*24;
  for(int i=0; i < len; i++) { 
    Out[i] = In1[i]+In2[i];
  }
#if 0
   double in10r;
   double in10i;
   double in11r;
   double in11i;
   double in12r;
   double in12i;

   double in20r;
   double in20i;
   double in21r;
   double in21i;
   double in22r;
   double in22i;

   double out0r;
   double out0i;
   double out1r;
   double out1i;
   double out2r;
   double out2i;

   int counter =0;
   int in1ptr =0;
   int in2ptr =0;
   int outptr =0;

  if( n_3vec > 0 ) {
    in10r = (double)In1[in1ptr++];
    in20r = (double)In2[in2ptr++];
    in10i = (double)In1[in1ptr++];
    in20i = (double)In2[in2ptr++];
    for(counter = 0; counter < n_3vec-1; counter++) { 
      out0r = in10r + in20r;
      Out[outptr++] = (REAL)out0r;

      in11r = (double)In1[in1ptr++];
      in21r = (double)In2[in2ptr++];
      out0i = in10i + in20i;
      Out[outptr++] = (REAL)out0i;

      in11i = (double)In1[in1ptr++];
      in21i = (double)In2[in2ptr++];
      out1r = in11r + in21r;
      Out[outptr++] = (REAL)out1r;

      in12r = (double)In1[in1ptr++];
      in22r = (double)In2[in2ptr++];
      out1i = in11i + in21i;
      Out[outptr++] = (REAL)out1i;

      in12i = (double)In1[in1ptr++];
      in22i = (double)In2[in2ptr++];
      out2r = in12r + in22r;
      Out[outptr++] = (REAL)out2r;

      in10r = (double)In1[in1ptr++];
      in20r = (double)In2[in2ptr++];     
      out2i = in12i + in22i;
      Out[outptr++] = (REAL)out2i;

      in10i = (double)In1[in1ptr++];
      in20i = (double)In2[in2ptr++];
    }
    out0r = in10r + in20r;
    Out[outptr++] = (REAL)out0r;

    in11r = (double)In1[in1ptr++];
    in21r = (double)In2[in2ptr++];
    out0i = in10i + in20i;
    Out[outptr++] = (REAL)out0i;

    in11i = (double)In1[in1ptr++];
    in21i = (double)In2[in2ptr++];
    out1r = in11r + in21r;
    Out[outptr++] = (REAL)out1r;

    in12r = (double)In1[in1ptr++];
    in22r = (double)In2[in2ptr++];
    out1i = in11i + in21i;
    Out[outptr++] = (REAL)out1i;

    in12i = (double)In1[in1ptr++];
    in22i = (double)In2[in2ptr++];
    out2r = in12r + in22r;
    Out[outptr++] = (REAL)out2r;
    out2i = in12i + in22i;
    Out[outptr++] = (REAL)out2i;
  }
#endif

}


} // namespace QDP;

#endif // guard
