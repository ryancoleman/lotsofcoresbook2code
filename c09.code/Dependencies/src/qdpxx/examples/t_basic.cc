// $Id: t_basic.cc,v 1.4 2004-08-11 18:53:10 edwards Exp $
/*! \file
 *  \brief Test some simple basic routines
 */

#include "qdp.h"

using namespace QDP;


int main(int argc, char *argv[])
{
  // Put the machine into a known state
  QDP_initialize(&argc, &argv);

  QDP_PUSH_PROFILE(QDP::getProfileLevel());

  // Setup the layout
  const int foo[] = {2,2,2,1};
  multi1d<int> nrow(Nd);
  nrow = foo;  // Use only Nd elements
  Layout::setLattSize(nrow);
  Layout::create();

  XMLFileWriter xml_out("t_basic.xml");
  push(xml_out, "t_basic");

  push(xml_out,"lattice");
  write(xml_out,"nrow", nrow);
  write(xml_out, "Nd", Nd);
  write(xml_out, "Nc", Nc);
  write(xml_out,"logicalSize",Layout::logicalSize());
  pop(xml_out);

  multi1d<LatticeColorMatrix> uu(Nd);
  LatticeColorMatrix   g;

  LatticeColorMatrix u;
  LatticeColorMatrix lctmp1;
  LatticeColorMatrix lctmp2;
  LatticeColorMatrix lctmp3;
  LatticePropagator q;
  LatticePropagator lqtmp1;
  LatticePropagator lqtmp2;
  LatticePropagator qtmp3;
  LatticeFermion    lftmp1;
  LatticeFermion    lftmp2;
  LatticeFermion    lftmp3;
  LatticeReal R1;
  LatticeReal Rtmp;
  LatticeComplex c;
  LatticeComplex ctmp;
  SpinMatrix s1;
  SpinMatrix s2;
  SpinMatrix s3;
  SpinMatrix s4;
  Real r1 = 1.7;
  Real r2 = 2.1;
  int i;
  int mu;
  int nu;
  Double dsum1;
  Double dsum2;
  Double dsum3;

  // RNG
  random(u);
  gaussian(lctmp1);
  push(xml_out,"LATTICE_COLORMATRIX_Site_Variables");
  write(xml_out, "u", u);
  write(xml_out, "lctmp1", lctmp1);
  pop(xml_out);

  // Colormat ops
  random(lctmp1); random(lctmp2); random(R1); random(c);
  push(xml_out,"Lattice_ColorMat_ops");
  random(lctmp3); lctmp3 = lctmp1 * lctmp2;
  write(xml_out, "C_X_C",lctmp3);
  random(lctmp3); lctmp3 += lctmp1 * lctmp2;
  write(xml_out, "C_peq_C_X_C",lctmp3);
  random(lctmp3); lctmp3 -= lctmp1 * lctmp2;
  write(xml_out, "C_meq_C_X_C",lctmp3);

  write(xml_out, "C_eq_C_x_C_pl_C",lctmp1*lctmp2+lctmp1);
  write(xml_out, "C_eq_aC_x_C",adj(lctmp1)*lctmp2);
  write(xml_out, "C_eq_C_x_aC",lctmp1*adj(lctmp2));
  write(xml_out, "C_eq_aC_x_aC",adj(lctmp1)*adj(lctmp2));

  write(xml_out, "C_eq_r_x_C",r1*lctmp1);

  random(lctmp3); lctmp3 += r1;
  write(xml_out, "C_peq_r",lctmp3);
  random(lctmp3); lctmp3 -= r1;
  write(xml_out, "C_meq_r",lctmp3);
  random(lctmp3); lctmp3 *=r1;
  write(xml_out, "C_teq_r",lctmp3);
  random(lctmp3); lctmp3 /= r1;
  write(xml_out, "C_deq_r",lctmp3);

  random(lctmp3); lctmp3 += R1;
  write(xml_out, "C_peq_R",lctmp3);
  random(lctmp3); lctmp3 -= R1;
  write(xml_out, "C_meq_R",lctmp3);
  random(lctmp3); lctmp3 *= R1;
  write(xml_out, "C_teq_R",lctmp3);
  random(lctmp3); lctmp3 /= R1;
  write(xml_out, "C_deq_R",lctmp3);

  random(lctmp3); lctmp3 += c;
  write(xml_out, "C_peq_c",lctmp3);
  random(lctmp3); lctmp3 -= c;
  write(xml_out, "C_meq_c",lctmp3);
  random(lctmp3); lctmp3 *= c;
  write(xml_out, "C_teq_c",lctmp3);
  random(lctmp3); lctmp3 /= c;
  write(xml_out, "C_deq_c",lctmp3);

  random(lctmp3); lctmp3 += lctmp1;
  write(xml_out, "C_peq_C",lctmp3);
  random(lctmp3); lctmp3 -= lctmp1;
  write(xml_out, "C_meq_C",lctmp3);

  write(xml_out, "C_eq_R_pl_aC_x_C",R1+adj(lctmp1)*lctmp2);
  write(xml_out, "C_eq_R_pl_C_x_aC",R1+lctmp1*adj(lctmp2));
  write(xml_out, "C_eq_R_pl_aC_x_aC",R1+adj(lctmp1)*adj(lctmp2));

  pop(xml_out);

  // Ferm mult
  random(lctmp1); random(lftmp1);
  push(xml_out,"Lattice_Ferm_ops");
  random(lftmp2); lftmp2 = lctmp1 * lftmp1;
  write(xml_out, "D_eq_C_x_D",lftmp2);
  random(lftmp2); lftmp2 += lctmp1 * lftmp1;
  write(xml_out, "D_peq_C_x_D",lftmp2);
  random(lftmp2); lftmp2 -= lctmp1 * lftmp1;
  write(xml_out, "D_meq_C_x_D",lftmp2);
  write(xml_out, "D_eq_aC_x_D",adj(lctmp1)*lftmp1);
  pop(xml_out);

  // Various site ops
  random(lftmp1); random(lftmp2);
  push(xml_out,"Site_functions");
  write(xml_out, "outerProduct",outerProduct(lftmp1,lftmp2));
  lctmp1 = traceSpin(outerProduct(lftmp1,lftmp2));
  write(xml_out, "C_eq_traceSpin_outerProduct",lctmp1);
  pop(xml_out);


#if 0

  mu = 0;
  nu = 1;

  // test 1
  lctmp2 = shift(u, FORWARD, mu) * lctmp1;
  push(xml_out,"MULTIPLY_MATRIX_Forward1_Fetched");
  write(xml_out, "mu", mu);
  write(xml_out, "lctmp2", lctmp2);
  pop(xml_out);

  /* test 2 */
  lctmp2 = shift(adj(u), FORWARD, mu) * lctmp1;
  push(xml_out,"MULTIPLY_MATRIX_Conj1_Forward1_Fetched");
  write(xml_out, "mu", mu);
  write(xml_out, "lctmp2", lctmp2);
  pop(xml_out);

  /* test 3 */
  lctmp2 = shift(adj(u), FORWARD, mu) * adj(lctmp1);
  push(xml_out,"MULTIPLY_MATRIX_Conj12_Forward1_Fetched");
  write(xml_out, "mu", mu);
  write(xml_out, "lctmp2", lctmp2);
  pop(xml_out);

  /* test 4 */
  lctmp2 = shift(adj(u), BACKWARD, nu) * shift(u, FORWARD, mu);
  push(xml_out,"MULTIPLY_MATRIX_Conj1_Comm12_Fetched");
  write(xml_out, "mu", mu);
  write(xml_out, "lctmp2", lctmp2);
  pop(xml_out);

  /* test 5 */
  rtmp = real(trace(lctmp1 * u));
  push(xml_out,"TRACE_MULTIPLY_MATRIX_realpart");
  write(xml_out, "rtmp", rtmp);
  pop(xml_out);

  /* test 6 */
  ctmp = trace(lctmp1 * u);
  push(xml_out,"TRACE_MULTIPLY_MATRIX_complexpart");
  write(xml_out, "ctmp", ctmp);
  pop(xml_out);

  /* test 7 */
  rtmp = real(trace(u));
  push(xml_out,"TRACE_MATRIX_realpart");
  write(xml_out, "rtmp", rtmp);
  pop(xml_out);

  /* test 8 */
  ctmp = trace(u);
  push(xml_out,"TRACE_MATRIX_complexpart");
  write(xml_out, "ctmp", ctmp);
  pop(xml_out);


  /* Now do tests on propagators */
  gaussian(q);
  gaussian(lqtmp1);
  gaussian(lqtmp2);
  push(xml_out,"LATTICE_PROPAGATOR_Site_Variables");
  write(xml_out, "q", q);
  write(xml_out, "lqtmp1", lqtmp1);
  pop(xml_out);

  /* test 9 */
  lqtmp2 = q * lqtmp1;
  push(xml_out,"MULTIPLY_PROP_replace");
  write(xml_out, "lqtmp2", lqtmp2);
  pop(xml_out);

  /* test 10 */
  gaussian(lqtmp2);
  lqtmp2 += q * lqtmp1;
  push(xml_out,"MULTIPLY_PROP_add");
  write(xml_out, "lqtmp2", lqtmp2);
  pop(xml_out);

  /* test 11 */
  lqtmp2 = adj(q) * lqtmp1;
  push(xml_out,"MULTIPLY_PROP_Conj1_replace");
  write(xml_out, "lqtmp2", lqtmp2);
  pop(xml_out);

  /* test 12 */
  gaussian(lqtmp2);
  lqtmp2 += adj(q) * lqtmp1;
  push(xml_out,"MULTIPLY_PROP_Conj1_add");
  write(xml_out, "lqtmp2", lqtmp2);
  pop(xml_out);

  /* test 13 */
  lqtmp2 = shift(q, FORWARD, mu) * lqtmp1;
  push(xml_out,"MULTIPLY_PROP_Forward1_Fetched");
  write(xml_out, "mu", mu);
  write(xml_out, "lqtmp2", lqtmp2);
  pop(xml_out);

  /* test 14 */
  lqtmp2 = shift(adj(q), FORWARD, mu) * adj(lqtmp1);
  push(xml_out,"MULTIPLY_PROP_Conj12_Forward1_Fetched");
  write(xml_out, "mu", mu);
  write(xml_out, "lqtmp2", lqtmp2);
  pop(xml_out);

  /* test 15 */
  lqtmp2 = shift(adj(q), BACKWARD, nu) * shift(lqtmp1, FORWARD, mu);
  push(xml_out,"MULTIPLY_PROP_Conj1_Comm12_Fetched");
  write(xml_out, "mu", mu);
  write(xml_out, "lqtmp2", lqtmp2);
  pop(xml_out);

  /* test 16 */
  lqtmp2 = u * q;
  push(xml_out,"MULTIPLY_PROP_U_front");
  write(xml_out, "lqtmp2", lqtmp2);
  pop(xml_out);

  /* test 17 */
  lqtmp2 = q * u;
  push(xml_out,"MULTIPLY_PROP_U_back");
  write(xml_out, "lqtmp2", lqtmp2);
  pop(xml_out);

  /* test 18 */
  dsum1 = norm2(q);
  push(xml_out,"SUMSQ_PROP");
  write(xml_out, "dsum1", dsum1);
  pop(xml_out);

  /* test 19 */
  int m = 1;
  lqtmp2 = Gamma(m) * q;
  dsum2 = norm2(lqtmp2);
  push(xml_out,"SPIN_PRODUCT_PROP");
  write(xml_out, "m", m);
  write(xml_out, "dsum2", dsum2);
  write(xml_out, "lqtmp2", lqtmp2);
  pop(xml_out);

  /* test 20 */
  lqtmp1 = Gamma(m) * lqtmp2;
  dsum1 = norm2(lqtmp1);
  push(xml_out,"SPIN_PRODUCT_PROP_again");
  write(xml_out, "m", m);
  write(xml_out, "dsum1", dsum1);
  write(xml_out, "lqtmp1", lqtmp1);
  pop(xml_out);

  /* test 21 */
  lqtmp1 += q;
  dsum1 = norm2(lqtmp1);
  push(xml_out,"SPIN_PROP_COPY");
  write(xml_out, "m", m);
  write(xml_out, "dsum1", dsum1);
  write(xml_out, "lqtmp1", lqtmp1);
  pop(xml_out);

  /* test 22 */
  s1 = 1;
  s2 = Gamma(m) * s1;
  s3 = Gamma(m) * s1;
  s4 = s1 * Gamma(m);
  push(xml_out,"SPIN_PROP_site_test");
  write(xml_out, "m", m);
  write(xml_out, "s2", s2);
  write(xml_out, "s3", s3);
  write(xml_out, "s4", s4);
  pop(xml_out);

  /* test 23 */
  lqtmp1 = q * s2;
  lqtmp2 = lqtmp1 * s2;
  lqtmp2 += q;
  dsum2 = norm2(lqtmp2);
  push(xml_out,"SPIN_PROP_site_mult");
  write(xml_out, "m", m);
  write(xml_out, "dsum2", dsum2);
  write(xml_out, "lqtmp2", lqtmp2);
  pop(xml_out);

  /* test 24 */
  lqtmp1 = q * Gamma(m);
  lqtmp2 = lqtmp1 * Gamma(m);
  lqtmp2 += q;
  dsum2 = norm2(lqtmp2);
  push(xml_out,"SPIN_PROP_site_mult_use_mult");
  write(xml_out, "m", m);
  write(xml_out, "dsum2", dsum2);
  write(xml_out, "lqtmp2", lqtmp2);
  pop(xml_out);

  /* test 25 */
  for(int m=0; m < Ns*Ns; ++m)
  {
    c = trace(q);

    s2 = Gamma(m) * s1;
    lqtmp1 = Gamma(m) * q;
    lqtmp2 = lqtmp1 * s2;
    lqtmp1 = lqtmp2;
    ctmp = trace(lqtmp2);
    ctmp -= c;
    dsum2 = norm2(ctmp);

    ctmp = trace(lqtmp1);
    ctmp += c;
    dsum3 = norm2(ctmp);
    push(xml_out,"SPIN_PROP_site_test_loop");
    write(xml_out, "m", m);
    write(xml_out, "dsum2", dsum2);
    write(xml_out, "dsum3", dsum3);
    pop(xml_out);
  }

  /* test 26 */
  for(int m=0; m < Ns*Ns; ++m)
    for(int n=0; n < Ns*Ns; ++n)
    {
      s2 = Gamma(m) * s1;
      s3 = Gamma(n) * s1;

      lqtmp1 = s2 * q;
      lqtmp2 = lqtmp1 * s3;
/*    lqtmp1 = q * s2; */

      lqtmp1 = Gamma(m) * q;
      qtmp3 = lqtmp1 * Gamma(n);
/*    qtmp3 = q * Gamma(n); */

      qtmp3 -= lqtmp2;
      dsum2 = norm2(qtmp3);
      push(xml_out,"SPIN_PROP_site_test_loop_mult");
      write(xml_out, "m", m);
      write(xml_out, "n", n);
      write(xml_out, "dsum2", dsum2);
      pop(xml_out);
    }


  /* test 27 */
  gaussian(lqtmp1);
  r = real(trace(adj(q) * lqtmp1));
  push(xml_out,"TRACE_MULT_PROP");
  write(xml_out, "r", r);
  pop(xml_out);

  /* test 28 */
  gaussian(lqtmp1);
  gaussian(lqtmp2);
  push(xml_out,"CONTRACT_PROP");
  write(xml_out, "lqtmp13", quarkContract13(lqtmp1, lqtmp2));
  write(xml_out, "lqtmp14", quarkContract14(lqtmp1, lqtmp2));
  write(xml_out, "lqtmp23", quarkContract23(lqtmp1, lqtmp2));
  write(xml_out, "lqtmp24", quarkContract24(lqtmp1, lqtmp2));
  pop(xml_out);

#endif


  pop(xml_out);
  xml_out.close();

  QDP_POP_PROFILE();

  // Time to bolt
  QDP_finalize();

  exit(0);
}


