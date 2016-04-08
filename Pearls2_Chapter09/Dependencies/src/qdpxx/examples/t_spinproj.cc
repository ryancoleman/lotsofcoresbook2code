// $Id: t_spinproj.cc,v 1.6 2007-02-22 03:30:27 bjoo Exp $

#include <iostream>
#include <iomanip>
#include <cstdio>

#include <time.h>

#include "qdp.h"

using namespace QDP;

void checkSpinProjDir0Plus(const REAL* src, REAL *dst, unsigned int n_vec);
void checkSpinProjDir0Minus(const REAL* src, REAL *dst, unsigned int n_vec);
void checkSpinProjDir1Plus(const REAL* src, REAL *dst, unsigned int n_vec);
void checkSpinProjDir1Minus(const REAL* src, REAL *dst, unsigned int n_vec);
void checkSpinProjDir2Plus(const REAL* src, REAL *dst, unsigned int n_vec);
void checkSpinProjDir2Minus(const REAL* src, REAL *dst, unsigned int n_vec);
void checkSpinProjDir3Plus(const REAL* src, REAL *dst, unsigned int n_vec);
void checkSpinProjDir3Minus(const REAL* src, REAL *dst, unsigned int n_vec);
void checkSpinReconDir0Plus(const REAL* src, REAL *dst, unsigned int n_vec);
void checkSpinReconDir0Minus(const REAL* src, REAL *dst, unsigned int n_vec); 
void checkSpinReconDir1Plus(const REAL* src, REAL *dst, unsigned int n_vec);
void checkSpinReconDir1Minus(const REAL* src, REAL *dst, unsigned int n_vec);
void checkSpinReconDir2Plus(const REAL* src, REAL *dst, unsigned int n_vec);
void checkSpinReconDir2Minus(const REAL* src, REAL *dst, unsigned int n_vec); 
void checkSpinReconDir3Plus(const REAL* src, REAL *dst, unsigned int n_vec);
void checkSpinReconDir3Minus(const REAL* src, REAL *dst, unsigned int n_vec); 
void checkAddSpinReconDir0Plus(const REAL* src, REAL *dst, unsigned int n_vec);
void checkAddSpinReconDir0Minus(const REAL* src, REAL *dst, unsigned int n_vec); 
void checkAddSpinReconDir1Plus(const REAL* src, REAL *dst, unsigned int n_vec);
void checkAddSpinReconDir1Minus(const REAL* src, REAL *dst, unsigned int n_vec);
void checkAddSpinReconDir2Plus(const REAL* src, REAL *dst, unsigned int n_vec);
void checkAddSpinReconDir2Minus(const REAL* src, REAL *dst, unsigned int n_vec); 
void checkAddSpinReconDir3Plus(const REAL* src, REAL *dst, unsigned int n_vec);
void checkAddSpinReconDir3Minus(const REAL* src, REAL *dst, unsigned int n_vec); 


int main(int argc, char *argv[])
{
  // Put the machine into a known state
  QDP_initialize(&argc, &argv);
  
  // Setup the layout
  const int foo[] = {4,4,4,4};
  multi1d<int> nrow(Nd);
  nrow = foo;  // Use only Nd elements
  Layout::setLattSize(nrow);
  Layout::create();

  LatticeHalfFermion H;
  HalfFermion h QDP_ALIGN16;
  Fermion v QDP_ALIGN16 ;

  LatticeHalfFermion H2,H3;
  LatticeFermion V;
  LatticeFermion V2;
  LatticeFermion V3;
  gaussian(V);
  LatticeHalfFermion diff;
  LatticeFermion diff_v;

  /*
   *  Proj 0 Plus
   */
  // Old Way
  checkSpinProjDir0Plus(&(V.elem(0).elem(0).elem(0).real()),
			&(H2.elem(0).elem(0).elem(0).real()),
			all.end()-all.start()+1);

  // New Way
  H = spinProjectDir0Plus(V);

  diff = H2 - H;


  QDPIO::cout << "Proj0+: || old - new || / || old || = " << sqrt(norm2(diff)/norm2(H2)) << endl;

  /*
   *  Proj 0 Minus
   */

  // Old Way
  checkSpinProjDir0Minus(&(V.elem(0).elem(0).elem(0).real()),
			&(H2.elem(0).elem(0).elem(0).real()),
			all.end()-all.start()+1);

  // New Way
  H = spinProjectDir0Minus(V);
  diff = H2 - H;
  QDPIO::cout << "Proj0-: || old - new || / || old || = " << sqrt(norm2(diff)/norm2(H2)) << endl;

  /*
   *  Proj 1 Plus
   */

  // Old Way
  checkSpinProjDir1Plus(&(V.elem(0).elem(0).elem(0).real()),
			&(H2.elem(0).elem(0).elem(0).real()),
			all.end()-all.start()+1);

  // New Way
  H = spinProjectDir1Plus(V);
  diff = H2 - H;


  QDPIO::cout << "Proj1+: || old - new || / || old || = " << sqrt(norm2(diff)/norm2(H2)) << endl;

  /*
   *  Proj 1 Minus
   */

  // Old Way
  checkSpinProjDir1Minus(&(V.elem(0).elem(0).elem(0).real()),
			&(H2.elem(0).elem(0).elem(0).real()),
			all.end()-all.start()+1);

  // New Way
  H = spinProjectDir1Minus(V);
  diff = H2 - H;
  QDPIO::cout << "Proj1-: || old - new || / || old || = " << sqrt(norm2(diff)/norm2(H2)) << endl;

  /*
   *  Proj 2 Plus
   */

  // Old Way
  checkSpinProjDir2Plus(&(V.elem(0).elem(0).elem(0).real()),
			&(H2.elem(0).elem(0).elem(0).real()),
			all.end()-all.start()+1);

  // New Way
  H = spinProjectDir2Plus(V);
  diff = H2 - H;


  QDPIO::cout << "Proj2+: || old - new || / || old || = " << sqrt(norm2(diff)/norm2(H2)) << endl;

  /*
   *  Proj 2 Minus
   */

  // Old Way
  checkSpinProjDir2Minus(&(V.elem(0).elem(0).elem(0).real()),
			&(H2.elem(0).elem(0).elem(0).real()),
			all.end()-all.start()+1);

  // New Way
  H = spinProjectDir2Minus(V);
  diff = H2 - H;
  QDPIO::cout << "Proj2-: || old - new || / || old || = " << sqrt(norm2(diff)/norm2(H2)) << endl;



  /*
   *  Proj 3 Plus
   */

  // Old Way
  checkSpinProjDir3Plus(&(V.elem(0).elem(0).elem(0).real()),
			&(H2.elem(0).elem(0).elem(0).real()),
			all.end()-all.start()+1);

  // New Way
  H = spinProjectDir3Plus(V);
  diff = H2 - H;


  QDPIO::cout << "Proj3+: || old - new || / || old || = " << sqrt(norm2(diff)/norm2(H2)) << endl;

  /*
   *  Proj 3 Minus
   */

  // Old Way
  checkSpinProjDir3Minus(&(V.elem(0).elem(0).elem(0).real()),
			&(H2.elem(0).elem(0).elem(0).real()),
			all.end()-all.start()+1);

  // New Way
  H = spinProjectDir3Minus(V);
  diff = H2 - H;
  QDPIO::cout << "Proj3-: || old - new || / || old || = " << sqrt(norm2(diff)/norm2(H2)) << endl;




  gaussian(H);

  /* 
   * Recon 0 Plus
   */
  
  // Old Way
  checkSpinReconDir0Plus(&(H.elem(0).elem(0).elem(0).real()),
			 &(V.elem(0).elem(0).elem(0).real()),
			 all.end()-all.start()+1);
  
  // New Way
  V2 = spinReconstructDir0Plus(H);
  diff_v = V - V2;
  QDPIO::cout << "Recon0+: || old - new || / || old || = " << sqrt(norm2(diff_v)/norm2(V2)) << endl;

  // Old Way
  checkSpinReconDir0Minus(&(H.elem(0).elem(0).elem(0).real()),
			 &(V.elem(0).elem(0).elem(0).real()),
			 all.end()-all.start()+1);
  
  // New Way
  V2 = spinReconstructDir0Minus(H);
  diff_v = V - V2;
  QDPIO::cout << "Recon0-: || old - new || / || old || = " << sqrt(norm2(diff_v)/norm2(V2)) << endl;

  /* 
   * Recon 1 Plus
   */
  
  // Old Way
  checkSpinReconDir1Plus(&(H.elem(0).elem(0).elem(0).real()),
			 &(V.elem(0).elem(0).elem(0).real()),
			 all.end()-all.start()+1);
  
  // New Way
  V2 = spinReconstructDir1Plus(H);
  diff_v = V - V2;
  QDPIO::cout << "Recon1+: || old - new || / || old || = " << sqrt(norm2(diff_v)/norm2(V2)) << endl;

  // Old Way
  checkSpinReconDir1Minus(&(H.elem(0).elem(0).elem(0).real()),
			 &(V.elem(0).elem(0).elem(0).real()),
			 all.end()-all.start()+1);
  
  // New Way
  V2 = spinReconstructDir1Minus(H);
  diff_v = V - V2;
  QDPIO::cout << "Recon1-: || old - new || / || old || = " << sqrt(norm2(diff_v)/norm2(V2)) << endl;



  /* 
   * Recon 2 Plus
   */
  
  // Old Way
  checkSpinReconDir2Plus(&(H.elem(0).elem(0).elem(0).real()),
			 &(V.elem(0).elem(0).elem(0).real()),
			 all.end()-all.start()+1);
  
  // New Way
  V2 = spinReconstructDir2Plus(H);
  diff_v = V - V2;
  QDPIO::cout << "Recon2+: || old - new || / || old || = " << sqrt(norm2(diff_v)/norm2(V2)) << endl;

  // Old Way
  checkSpinReconDir2Minus(&(H.elem(0).elem(0).elem(0).real()),
			 &(V.elem(0).elem(0).elem(0).real()),
			 all.end()-all.start()+1);
  
  // New Way
  V2 = spinReconstructDir2Minus(H);
  diff_v = V - V2;
  QDPIO::cout << "Recon2-: || old - new || / || old || = " << sqrt(norm2(diff_v)/norm2(V2)) << endl;

  /* 
   * Recon 3 Plus
   */
  
  // Old Way
  checkSpinReconDir3Plus(&(H.elem(0).elem(0).elem(0).real()),
			 &(V.elem(0).elem(0).elem(0).real()),
			 all.end()-all.start()+1);
  
  // New Way
  V2 = spinReconstructDir3Plus(H);
  diff_v = V - V2;
  QDPIO::cout << "Recon3+: || old - new || / || old || = " << sqrt(norm2(diff_v)/norm2(V2)) << endl;

  // Old Way
  checkSpinReconDir3Minus(&(H.elem(0).elem(0).elem(0).real()),
			 &(V.elem(0).elem(0).elem(0).real()),
			 all.end()-all.start()+1);
  
  // New Way
  V2 = spinReconstructDir3Minus(H);
  diff_v = V - V2;
  QDPIO::cout << "Recon3-: || old - new || / || old || = " << sqrt(norm2(diff_v)/norm2(V2)) << endl;






  /* 
   * Add Recon 0 Plus
   */
  
  // Old Way
  gaussian(V);
  V2=V;

  checkAddSpinReconDir0Plus(&(H.elem(0).elem(0).elem(0).real()),
			    &(V2.elem(0).elem(0).elem(0).real()),
			    all.end()-all.start()+1);

  
  // New Way
  V3 = V;
  V3 += spinReconstructDir0Plus(H);
  diff_v = V2 - V3;
  QDPIO::cout << "AddRecon0+: || old - new || / || old || = " << sqrt(norm2(diff_v)/norm2(V3)) << endl;

  /* 
   * Add Recon 0 Minus
   */
  
  // Old Way
  gaussian(V);
  V2=V;

  checkAddSpinReconDir0Minus(&(H.elem(0).elem(0).elem(0).real()),
			     &(V2.elem(0).elem(0).elem(0).real()),
			     all.end()-all.start()+1);
  
  
  // New Way
  V3 = V;
  V3 += spinReconstructDir0Minus(H);
  diff_v = V2 - V3;
  QDPIO::cout << "AddRecon0-: || old - new || / || old || = " << sqrt(norm2(diff_v)/norm2(V3)) << endl;

  /* 
   * Add Recon 1 Plus
   */
  
  // Old Way
  gaussian(V);
  V2=V;

  checkAddSpinReconDir1Plus(&(H.elem(0).elem(0).elem(0).real()),
			    &(V2.elem(0).elem(0).elem(0).real()),
			    all.end()-all.start()+1);

  
  // New Way
  V3 = V;
  V3 += spinReconstructDir1Plus(H);
  diff_v = V2 - V3;
  QDPIO::cout << "AddRecon1+: || old - new || / || old || = " << sqrt(norm2(diff_v)/norm2(V3)) << endl;

  /* 
   * Add Recon 1 Minus
   */
  
  // Old Way
  gaussian(V);
  V2=V;

  checkAddSpinReconDir1Minus(&(H.elem(0).elem(0).elem(0).real()),
			     &(V2.elem(0).elem(0).elem(0).real()),
			     all.end()-all.start()+1);
  
  
  // New Way
  V3 = V;
  V3 += spinReconstructDir1Minus(H);
  diff_v = V2 - V3;
  QDPIO::cout << "AddRecon1-: || old - new || / || old || = " << sqrt(norm2(diff_v)/norm2(V3)) << endl;


  /* 
   * Add Recon 2 Plus
   */
  
  // Old Way
  gaussian(V);
  V2=V;

  checkAddSpinReconDir2Plus(&(H.elem(0).elem(0).elem(0).real()),
			    &(V2.elem(0).elem(0).elem(0).real()),
			    all.end()-all.start()+1);

  
  // New Way
  V3 = V;
  V3 += spinReconstructDir2Plus(H);
  diff_v = V2 - V3;
  QDPIO::cout << "AddRecon2+: || old - new || / || old || = " << sqrt(norm2(diff_v)/norm2(V3)) << endl;

  /* 
   * Add Recon 2 Minus
   */
  
  // Old Way
  gaussian(V);
  V2=V;

  checkAddSpinReconDir2Minus(&(H.elem(0).elem(0).elem(0).real()),
			     &(V2.elem(0).elem(0).elem(0).real()),
			     all.end()-all.start()+1);
  
  
  // New Way
  V3 = V;
  V3 += spinReconstructDir2Minus(H);
  diff_v = V2 - V3;
  QDPIO::cout << "AddRecon2-: || old - new || / || old || = " << sqrt(norm2(diff_v)/norm2(V3)) << endl;


  /* 
   * Add Recon 3 Plus
   */
  
  // Old Way
  gaussian(V);
  V2=V;

  checkAddSpinReconDir3Plus(&(H.elem(0).elem(0).elem(0).real()),
			    &(V2.elem(0).elem(0).elem(0).real()),
			    all.end()-all.start()+1);

  
  // New Way
  V3 = V;
  V3 += spinReconstructDir3Plus(H);
  diff_v = V2 - V3;
  QDPIO::cout << "AddRecon3+: || old - new || / || old || = " << sqrt(norm2(diff_v)/norm2(V3)) << endl;

  /* 
   * Add Recon 1 Minus
   */
  
  // Old Way
  gaussian(V);
  V2=V;

  checkAddSpinReconDir3Minus(&(H.elem(0).elem(0).elem(0).real()),
			     &(V2.elem(0).elem(0).elem(0).real()),
			     all.end()-all.start()+1);
  
  
  // New Way
  V3 = V;
  V3 += spinReconstructDir3Minus(H);
  diff_v = V2 - V3;
  QDPIO::cout << "AddRecon3-: || old - new || / || old || = " << sqrt(norm2(diff_v)/norm2(V3)) << endl;




  LatticeColorMatrix u;
  gaussian(u);
  gaussian(V);

  ColorMatrix colmat;
  colmat.elem() = u.elem(0).elem();
  Fermion ferm QDP_ALIGN16;
  ferm.elem() = V.elem(0);

  Fermion res QDP_ALIGN16;
  res = adj(colmat)*ferm;
  
  Fermion res2 QDP_ALIGN16;
 #if 0
  _inline_mult_adj_su3_mat_vec( colmat.elem().elem(),
				ferm.elem().elem(0),
				res2.elem().elem(0));
  
  _inline_mult_adj_su3_mat_vec( colmat.elem().elem(),
				ferm.elem().elem(1),
				res2.elem().elem(1));

  _inline_mult_adj_su3_mat_vec( colmat.elem().elem(),
				ferm.elem().elem(2),
				res2.elem().elem(2));

  _inline_mult_adj_su3_mat_vec( colmat.elem().elem(),
				ferm.elem().elem(3),
				res2.elem().elem(3));

   Fermion diff_ferm QDP_ALIGN16;
   diff_ferm= res - res2;
   QDPIO::cout << "Diff Ferm = " << diff_ferm << endl;
 #endif
  
  H = adj(u)*spinProjectDir0Plus(V);
  H2= spinProjectDir0Plus(V);
  H3= adj(u)*H2;

  diff = H3 - H;
  QDPIO::cout << "AdjProj0+: || old - new || / || old || = " << sqrt(norm2(diff)/norm2(H)) << endl;
  
  H = adj(u)*spinProjectDir0Minus(V);
  H2= spinProjectDir0Minus(V);
  H3= adj(u)*H2;

  diff = H3 - H;
  QDPIO::cout << "AdjProj0-: || old - new || / || old || = " << sqrt(norm2(diff)/norm2(H)) << endl;

  H = adj(u)*spinProjectDir1Plus(V);
  H2= spinProjectDir1Plus(V);
  H3= adj(u)*H2;

  diff = H3 - H;
  QDPIO::cout << "AdjProj1+: || old - new || / || old || = " << sqrt(norm2(diff)/norm2(H)) << endl;
  
  H = adj(u)*spinProjectDir1Minus(V);
  H2= spinProjectDir1Minus(V);
  H3= adj(u)*H2;

  diff = H3 - H;
  QDPIO::cout << "AdjProj1-: || old - new || / || old || = " << sqrt(norm2(diff)/norm2(H)) << endl;

  H = adj(u)*spinProjectDir2Plus(V);
  H2= spinProjectDir2Plus(V);
  H3= adj(u)*H2;

  diff = H3 - H;
  QDPIO::cout << "AdjProj2+: || old - new || / || old || = " << sqrt(norm2(diff)/norm2(H)) << endl;
  
  H = adj(u)*spinProjectDir2Minus(V);
  H2= spinProjectDir2Minus(V);
  H3= adj(u)*H2;

  diff = H3 - H;
  QDPIO::cout << "AdjProj2-: || old - new || / || old || = " << sqrt(norm2(diff)/norm2(H)) << endl;

  H = adj(u)*spinProjectDir3Plus(V);
  H2= spinProjectDir3Plus(V);
  H3= adj(u)*H2;

  diff = H3 - H;
  QDPIO::cout << "AdjProj3+: || old - new || / || old || = " << sqrt(norm2(diff)/norm2(H)) << endl;
  
  H = adj(u)*spinProjectDir3Minus(V);
  H2= spinProjectDir3Minus(V);
  H3= adj(u)*H2;

  diff = H3 - H;
  QDPIO::cout << "AdjProj3-: || old - new || / || old || = " << sqrt(norm2(diff)/norm2(H)) << endl;


  gaussian(H);
  V = spinReconstructDir0Plus( u*H );
  H2 = u*H;
  V2 = spinReconstructDir0Plus( H2 );
  diff_v = V2 - V;
  QDPIO::cout << "ReconUPsiDir0+: || old - new || / || old || = " << sqrt(norm2(diff_v)/norm2(V2)) << endl;

  gaussian(H);
  V = spinReconstructDir0Minus( u*H );
  H2 = u*H;
  V2 = spinReconstructDir0Minus( H2 );
  diff_v = V2 - V;
  QDPIO::cout << "ReconUPsiDir0-: || old - new || / || old || = " << sqrt(norm2(diff_v)/norm2(V2)) << endl;

  gaussian(H);
  V = spinReconstructDir1Plus( u*H );
  H2 = u*H;
  V2 = spinReconstructDir1Plus( H2 );
  diff_v = V2 - V;
  QDPIO::cout << "ReconUPsiDir1+: || old - new || / || old || = " << sqrt(norm2(diff_v)/norm2(V2)) << endl;

  gaussian(H);
  V = spinReconstructDir1Minus( u*H );
  H2 = u*H;
  V2 = spinReconstructDir1Minus( H2 );
  diff_v = V2 - V;
  QDPIO::cout << "ReconUPsiDir1-: || old - new || / || old || = " << sqrt(norm2(diff_v)/norm2(V2)) << endl;

  gaussian(H);
  V = spinReconstructDir2Plus( u*H );
  H2 = u*H;
  V2 = spinReconstructDir2Plus( H2 );
  diff_v = V2 - V;
  QDPIO::cout << "ReconUPsiDir2+: || old - new || / || old || = " << sqrt(norm2(diff_v)/norm2(V2)) << endl;

  gaussian(H);
  V = spinReconstructDir2Minus( u*H );
  H2 = u*H;
  V2 = spinReconstructDir2Minus( H2 );
  diff_v = V2 - V;
  QDPIO::cout << "ReconUPsiDir2-: || old - new || / || old || = " << sqrt(norm2(diff_v)/norm2(V2)) << endl;

  gaussian(H);
  V = spinReconstructDir3Plus( u*H );
  H2 = u*H;
  V2 = spinReconstructDir3Plus( H2 );
  diff_v = V2 - V;
  QDPIO::cout << "ReconUPsiDir3+: || old - new || / || old || = " << sqrt(norm2(diff_v)/norm2(V2)) << endl;

  gaussian(H);
  V = spinReconstructDir3Minus( u*H );
  H2 = u*H;
  V2 = spinReconstructDir3Minus( H2 );
  diff_v = V2 - V;
  QDPIO::cout << "ReconUPsiDir3-: || old - new || / || old || = " << sqrt(norm2(diff_v)/norm2(V2)) << endl;




  gaussian(H);
  gaussian(V);
  V2 = V;

  V += spinReconstructDir0Plus( u*H );
  H2 = u*H;
  V2 += spinReconstructDir0Plus( H2 );
  diff_v = V2 - V;
  QDPIO::cout << "ReconUPsiDir0+=: || old - new || / || old || = " << sqrt(norm2(diff_v)/norm2(V2)) << endl;

  gaussian(H);
  gaussian(V);
  V2=V;

  V += spinReconstructDir0Minus( u*H );
  H2 = u*H;
  V2 += spinReconstructDir0Minus( H2 );
  diff_v = V2 - V;
  QDPIO::cout << "ReconUPsiDir0-=: || old - new || / || old || = " << sqrt(norm2(diff_v)/norm2(V2)) << endl;

  gaussian(H);
  gaussian(V);
  V2 = V;

  V += spinReconstructDir1Plus( u*H );
  H2 = u*H;
  V2 += spinReconstructDir1Plus( H2 );
  diff_v = V2 - V;
  QDPIO::cout << "ReconUPsiDir1+=: || old - new || / || old || = " << sqrt(norm2(diff_v)/norm2(V2)) << endl;

  gaussian(H);
  gaussian(V);
  V2=V;

  V += spinReconstructDir1Minus( u*H );
  H2 = u*H;
  V2 += spinReconstructDir1Minus( H2 );
  diff_v = V2 - V;
  QDPIO::cout << "ReconUPsiDir1-=: || old - new || / || old || = " << sqrt(norm2(diff_v)/norm2(V2)) << endl;

  gaussian(H);
  gaussian(V);
  V2 = V;

  V += spinReconstructDir2Plus( u*H );
  H2 = u*H;
  V2 += spinReconstructDir2Plus( H2 );
  diff_v = V2 - V;
  QDPIO::cout << "ReconUPsiDir2+=: || old - new || / || old || = " << sqrt(norm2(diff_v)/norm2(V2)) << endl;

  gaussian(H);
  gaussian(V);
  V2=V;

  V += spinReconstructDir2Minus( u*H );
  H2 = u*H;
  V2 += spinReconstructDir2Minus( H2 );
  diff_v = V2 - V;
  QDPIO::cout << "ReconUPsiDir2-=: || old - new || / || old || = " << sqrt(norm2(diff_v)/norm2(V2)) << endl;

  gaussian(H);
  gaussian(V);
  V2=V;

  V += spinReconstructDir3Plus( u*H );
  H2 = u*H;
  V2 += spinReconstructDir3Plus( H2 );
  diff_v = V2 - V;
  QDPIO::cout << "ReconUPsiDir3+=: || old - new || / || old || = " << sqrt(norm2(diff_v)/norm2(V2)) << endl;

  gaussian(H);
  gaussian(V);
  V2=V;

  V += spinReconstructDir3Minus( u*H );
  H2 = u*H;
  V2 += spinReconstructDir3Minus( H2 );
  diff_v = V2 - V;
  QDPIO::cout << "ReconUPsiDir3-=: || old - new || / || old || = " << sqrt(norm2(diff_v)/norm2(V2)) << endl;


  // Time to bolt
  QDP_finalize();

  exit(0);
}


void checkSpinProjDir0Plus(const REAL* src, REAL *dst, unsigned int n_vec) 
{

#ifdef DEBUG_GENERIC_SPIN_PROJ_INLINES
  QDPIO::cout << "checkSpinProjDir0Plus" << endl;
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
void checkSpinProjDir0Minus(const REAL* src, REAL *dst, unsigned int n_vec) 
{

#ifdef DEBUG_GENERIC_SPIN_PROJ_INLINES
  QDPIO::cout << "checkSpinProjDir0Minus" << endl;
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
  
  for(int site=0; site < n_vec; site++) {
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
void checkSpinProjDir1Plus(const REAL* src, REAL *dst, unsigned int n_vec) 
{

#ifdef DEBUG_GENERIC_SPIN_PROJ_INLINES
  QDPIO::cout << "checkSpinProjDir1Plus" << endl;
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

  for(int site=0; site < n_vec; site++) {
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
void checkSpinProjDir1Minus(const REAL* src, REAL *dst, unsigned int n_vec) 
{

#ifdef DEBUG_GENERIC_SPIN_PROJ_INLINES
  QDPIO::cout << "checkSpinProjDir1Minus" << endl;
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
  
  for(int site=0; site < n_vec; site++) {
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
void checkSpinProjDir2Plus(const REAL* src, REAL *dst, unsigned int n_vec) 
{

#ifdef DEBUG_GENERIC_SPIN_PROJ_INLINES
  QDPIO::cout << "checkSpinProjDir2Plus" << endl;
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

  for(int site=0; site < n_vec; site++) {
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
void checkSpinProjDir2Minus(const REAL* src, REAL *dst, unsigned int n_vec) 
{

#ifdef DEBUG_GENERIC_SPIN_PROJ_INLINES
  QDPIO::cout << "checkSpinProjDir2Minus" << endl;
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
  
  for(int site=0; site < n_vec; site++) {
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
void checkSpinProjDir3Plus(const REAL* src, REAL *dst, unsigned int n_vec) 
{

#ifdef DEBUG_GENERIC_SPIN_PROJ_INLINES
  QDPIO::cout << "checkSpinProjDir3Plus" << endl;
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

  for(int site=0; site < n_vec; site++) {
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
void checkSpinProjDir3Minus(const REAL* src, REAL *dst, unsigned int n_vec) 
{

#ifdef DEBUG_GENERIC_SPIN_PROJ_INLINES
  QDPIO::cout << "checkSpinProjDir3Minus" << endl;
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
  
  for(int site=0; site < n_vec; site++) {
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


  
void checkSpinReconDir0Plus(const REAL* src, REAL *dst, unsigned int n_vec) 
{

#ifdef DEBUG_GENERIC_SPIN_RECON_S
  QDPIO::cout << "checkSpinProjDir0Plus" << endl;
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

  for(int site=0; site < n_vec; site++) {
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
  
void checkSpinReconDir0Minus(const REAL* src, REAL *dst, unsigned int n_vec) 
{

#ifdef DEBUG_GENERIC_SPIN_RECON_S
  QDPIO::cout << "checkSpinReconDir0Minus" << endl;
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

  for(int site=0; site < n_vec; site++) {
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
  
void checkSpinReconDir1Plus(const REAL* src, REAL *dst, unsigned int n_vec) 
{

#ifdef DEBUG_GENERIC_SPIN_RECON_S
  QDPIO::cout << "checkSpinReconDir0Plus" << endl;
#endif



  REAL tmp_hspinor[4][3][2];

  const REAL* src_shadow = src;
  REAL* dst_shadow = dst;

  const int re = 0;
  const int im = 1;
  const int Ncmpx = 2;
  const int Nsby2 = 2;

  for(int site=0; site < n_vec; site++) {
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
  
void checkSpinReconDir1Minus(const REAL* src, REAL *dst, unsigned int n_vec) 
{

#ifdef DEBUG_GENERIC_SPIN_RECON_S
  QDPIO::cout << "checkSpinReconDir0Minus" << endl;
#endif

  
  REAL tmp_hspinor[4][3][2];

  const REAL* src_shadow = src;
  REAL* dst_shadow = dst;

  const int re = 0;
  const int im = 1;
  const int Ncmpx = 2;
  const int Nsby2 = 2;

  for(int site=0; site < n_vec; site++) {
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
  
void checkSpinReconDir2Plus(const REAL* src, REAL *dst, unsigned int n_vec) 
{

#ifdef DEBUG_GENERIC_SPIN_RECON_S
  QDPIO::cout << "checkSpinReconDir0Plus" << endl;
#endif



  REAL tmp_hspinor[4][3][2];

  const REAL* src_shadow = src;
  REAL* dst_shadow = dst;

  const int re = 0;
  const int im = 1;
  const int Ncmpx = 2;
  const int Nsby2 = 2;

  for(int site=0; site < n_vec; site++) {
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
  
void checkSpinReconDir2Minus(const REAL* src, REAL *dst, unsigned int n_vec) 
{

#ifdef DEBUG_GENERIC_SPIN_RECON_S
  QDPIO::cout << "checkSpinReconDir0Minus" << endl;
#endif

  
  REAL tmp_hspinor[4][3][2];

  const REAL* src_shadow = src;
  REAL* dst_shadow = dst;

  const int re = 0;
  const int im = 1;
  const int Ncmpx = 2;
  const int Nsby2 = 2;

  for(int site=0; site < n_vec; site++) {
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
  
void checkSpinReconDir3Plus(const REAL* src, REAL *dst, unsigned int n_vec) 
{

#ifdef DEBUG_GENERIC_SPIN_RECON_S
  QDPIO::cout << "checkSpinReconDir0Plus" << endl;
#endif



  REAL tmp_hspinor[4][3][2];

  const REAL* src_shadow = src;
  REAL* dst_shadow = dst;

  const int re = 0;
  const int im = 1;
  const int Ncmpx = 2;
  const int Nsby2 = 2;

  for(int site=0; site < n_vec; site++) {
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
  
void checkSpinReconDir3Minus(const REAL* src, REAL *dst, unsigned int n_vec) 
{

#ifdef DEBUG_GENERIC_SPIN_RECON_S
  QDPIO::cout << "checkSpinReconDir0Minus" << endl;
#endif

  
  REAL tmp_hspinor[4][3][2];

  const REAL* src_shadow = src;
  REAL* dst_shadow = dst;

  const int re = 0;
  const int im = 1;
  const int Ncmpx = 2;
  const int Nsby2 = 2;

  for(int site=0; site < n_vec; site++) {
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
  
void checkAddSpinReconDir0Plus(const REAL* src, REAL *dst, unsigned int n_vec) 
{

#ifdef DEBUG_GENERIC_SPIN_RECON_S
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

  for(int site=0; site < n_vec; site++) {
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
  
void checkAddSpinReconDir0Minus(const REAL* src, REAL *dst, unsigned int n_vec) 
{

#ifdef DEBUG_GENERIC_SPIN_RECON_S
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

  for(int site=0; site < n_vec; site++) {
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
  
void checkAddSpinReconDir1Plus(const REAL* src, REAL *dst, unsigned int n_vec) 
{

#ifdef DEBUG_GENERIC_SPIN_RECON_S
  QDPIO::cout << "inlinaAddSpinReconDir0Plus" << endl;
#endif



  REAL tmp_hspinor[4][3][2];

  const REAL* src_shadow = src;
  REAL* dst_shadow = dst;

  const int re = 0;
  const int im = 1;
  const int Ncmpx = 2;
  const int Nsby2 = 2;

  for(int site=0; site < n_vec; site++) {
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
  
void checkAddSpinReconDir1Minus(const REAL* src, REAL *dst, unsigned int n_vec) 
{

#ifdef DEBUG_GENERIC_SPIN_RECON_S
  QDPIO::cout << "inlinaAddSpinReconDir0Minus" << endl;
#endif

  
  REAL tmp_hspinor[4][3][2];

  const REAL* src_shadow = src;
  REAL* dst_shadow = dst;

  const int re = 0;
  const int im = 1;
  const int Ncmpx = 2;
  const int Nsby2 = 2;

  for(int site=0; site < n_vec; site++) {
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
  
void checkAddSpinReconDir2Plus(const REAL* src, REAL *dst, unsigned int n_vec) 
{

#ifdef DEBUG_GENERIC_SPIN_RECON_S
  QDPIO::cout << "inlinaAddSpinReconDir0Plus" << endl;
#endif



  REAL tmp_hspinor[4][3][2];

  const REAL* src_shadow = src;
  REAL* dst_shadow = dst;

  const int re = 0;
  const int im = 1;
  const int Ncmpx = 2;
  const int Nsby2 = 2;

  for(int site=0; site < n_vec; site++) {
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
  
void checkAddSpinReconDir2Minus(const REAL* src, REAL *dst, unsigned int n_vec) 
{

#ifdef DEBUG_GENERIC_SPIN_RECON_S
  QDPIO::cout << "inlinaAddSpinReconDir0Minus" << endl;
#endif

  
  REAL tmp_hspinor[4][3][2];

  const REAL* src_shadow = src;
  REAL* dst_shadow = dst;

  const int re = 0;
  const int im = 1;
  const int Ncmpx = 2;
  const int Nsby2 = 2;

  for(int site=0; site < n_vec; site++) {
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
  
void checkAddSpinReconDir3Plus(const REAL* src, REAL *dst, unsigned int n_vec) 
{

#ifdef DEBUG_GENERIC_SPIN_RECON_S
  QDPIO::cout << "inlinaAddSpinReconDir0Plus" << endl;
#endif



  REAL tmp_hspinor[4][3][2];

  const REAL* src_shadow = src;
  REAL* dst_shadow = dst;

  const int re = 0;
  const int im = 1;
  const int Ncmpx = 2;
  const int Nsby2 = 2;

  for(int site=0; site < n_vec; site++) {
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
  
void checkAddSpinReconDir3Minus(const REAL* src, REAL *dst, unsigned int n_vec) 
{

#ifdef DEBUG_GENERIC_SPIN_RECON_S
  QDPIO::cout << "inlinaAddSpinReconDir0Minus" << endl;
#endif

  
  REAL tmp_hspinor[4][3][2];

  const REAL* src_shadow = src;
  REAL* dst_shadow = dst;

  const int re = 0;
  const int im = 1;
  const int Ncmpx = 2;
  const int Nsby2 = 2;

  for(int site=0; site < n_vec; site++) {
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


  
  
   
