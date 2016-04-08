// $Id: t_spinproj2.cc,v 1.3 2007-02-09 20:35:46 bjoo Exp $

#include <iostream>
#include <iomanip>
#include <cstdio>

#include <time.h>

#include "qdp.h"


using namespace QDP;

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

  // Lattice Version so I can gaussian
  LatticeColorMatrix u_lat;
  LatticeFermion ferm_lat, ferm2_lat, diff_lat;
  LatticeHalfFermion hferm_lat, hferm2_lat;

  
  gaussian(u_lat);
  gaussian(ferm_lat);
  gaussian(hferm_lat);

  // Scalar versions 

  // ColorMatrix 
  OScalar<PScalar<PColorMatrix< RComplex<REAL>, 3 > > > u, u2, u3;

  // Fermion
  OScalar<PSpinVector<PColorVector< RComplex<REAL>, 3>, 4> >  v,v2,v3,diff_v;

  // HalfFermion
  OScalar<PSpinVector<PColorVector< RComplex<REAL>, 3>, 2> > hv,hv2,hv3,diff_hv QDP_ALIGN16;

  u.elem() = u_lat.elem(0);
  u2.elem() = u_lat.elem(1);
  u3.elem() = u_lat.elem(2);
  hv.elem() = hferm_lat.elem(0);
  hv2.elem() = hferm_lat.elem(1);
  hv3.elem() = hferm_lat.elem(2);
  v.elem() = ferm_lat.elem(0);
 
  hv3 = adj(u)*spinProjectDir0Plus(v);
  hv  = spinProjectDir0Plus(v);
  hv2 = adj(u)*hv;

  diff_hv = hv2 - hv3;
  QDPIO::cout << "adjProjDir0Plus: || old - new || / || old || = " << sqrt( norm2(diff_hv)/norm2(hv2)) << endl;

  hv3 = adj(u)*spinProjectDir0Minus(v);
  hv  = spinProjectDir0Minus(v);
  hv2 = adj(u)*hv;

  diff_hv = hv2 - hv3;
  QDPIO::cout << "adjProjDir0Minus: || old - new || / || old || = " << sqrt( norm2(diff_hv)/norm2(hv2)) << endl;

  hv3 = adj(u)*spinProjectDir1Plus(v);
  hv  = spinProjectDir1Plus(v);
  hv2 = adj(u)*hv;

  diff_hv = hv2 - hv3;
  QDPIO::cout << "adjProjDir1Plus: || old - new || / || old || = " << sqrt( norm2(diff_hv)/norm2(hv2)) << endl;

  hv3 = adj(u)*spinProjectDir1Minus(v);
  hv  = spinProjectDir1Minus(v);
  hv2 = adj(u)*hv;

  diff_hv = hv2 - hv3;
  QDPIO::cout << "adjProjDir1Minus: || old - new || / || old || = " << sqrt( norm2(diff_hv)/norm2(hv2)) << endl;

  hv3 = adj(u)*spinProjectDir2Plus(v);
  hv  = spinProjectDir2Plus(v);
  hv2 = adj(u)*hv;

  diff_hv = hv2 - hv3;
  QDPIO::cout << "adjProjDir2Plus: || old - new || / || old || = " << sqrt( norm2(diff_hv)/norm2(hv2)) << endl;

  hv3 = adj(u)*spinProjectDir2Minus(v);
  hv  = spinProjectDir2Minus(v);
  hv2 = adj(u)*hv;

  diff_hv = hv2 - hv3;
  QDPIO::cout << "adjProjDir2Minus: || old - new || / || old || = " << sqrt( norm2(diff_hv)/norm2(hv2)) << endl;

  hv3 = adj(u)*spinProjectDir3Plus(v);
  hv  = spinProjectDir3Plus(v);
  hv2 = adj(u)*hv;

  diff_hv = hv2 - hv3;
  QDPIO::cout << "adjProjDir3Plus: || old - new || / || old || = " << sqrt( norm2(diff_hv)/norm2(hv2)) << endl;

  hv3 = adj(u)*spinProjectDir3Minus(v);
  hv  = spinProjectDir3Minus(v);
  hv2 = adj(u)*hv;

  diff_hv = hv2 - hv3;
  QDPIO::cout << "adjProjDir3Minus: || old - new || / || old || = " << sqrt( norm2(diff_hv)/norm2(hv2)) << endl;
  
  hv.elem() = hferm_lat.elem(0);
  v=spinReconstructDir0Minus(u*hv);
  hv2  = u *hv;
  v2   = spinReconstructDir0Minus(hv2);
  diff_v = v2 - v;
  QDPIO::cout << "reconProdDir0Minus: || old - new || / || old || = " << sqrt( norm2(diff_v)/norm2(v2)) << endl;


  hv.elem() = hferm_lat.elem(0);
  v=spinReconstructDir0Plus(u*hv);
  hv2  = u *hv;
  v2   = spinReconstructDir0Plus(hv2);
  diff_v = v2 - v;
  QDPIO::cout << "reconProdDir0Plus: || old - new || / || old || = " << sqrt( norm2(diff_v)/norm2(v2)) << endl;

  hv.elem() = hferm_lat.elem(0);
  v=spinReconstructDir1Minus(u*hv);
  hv2  = u *hv;
  v2   = spinReconstructDir1Minus(hv2);
  diff_v = v2 - v;
  QDPIO::cout << "reconProdDir1Minus: || old - new || / || old || = " << sqrt( norm2(diff_v)/norm2(v2)) << endl;


  hv.elem() = hferm_lat.elem(0);
  v=spinReconstructDir1Plus(u*hv);
  hv2  = u *hv;
  v2   = spinReconstructDir1Plus(hv2);
  diff_v = v2 - v;
  QDPIO::cout << "reconProdDir1Plus: || old - new || / || old || = " << sqrt( norm2(diff_v)/norm2(v2)) << endl;

  hv.elem() = hferm_lat.elem(0);
  v=spinReconstructDir2Minus(u*hv);
  hv2  = u *hv;
  v2   = spinReconstructDir2Minus(hv2);
  diff_v = v2 - v;
  QDPIO::cout << "reconProdDir2Minus: || old - new || / || old || = " << sqrt( norm2(diff_v)/norm2(v2)) << endl;


  hv.elem() = hferm_lat.elem(0);
  v=spinReconstructDir2Plus(u*hv);
  hv2  = u *hv;
  v2   = spinReconstructDir2Plus(hv2);
  diff_v = v2 - v;
  QDPIO::cout << "reconProdDir2Plus: || old - new || / || old || = " << sqrt( norm2(diff_v)/norm2(v2)) << endl;


  hv.elem() = hferm_lat.elem(0);
  v=spinReconstructDir2Minus(u*hv);
  hv2  = u *hv;
  v2   = spinReconstructDir2Minus(hv2);
  diff_v = v2 - v;
  QDPIO::cout << "reconProdDir3Minus: || old - new || / || old || = " << sqrt( norm2(diff_v)/norm2(v2)) << endl;


  hv.elem() = hferm_lat.elem(0);
  v=spinReconstructDir2Plus(u*hv);
  hv2  = u *hv;
  v2   = spinReconstructDir2Plus(hv2);
  diff_v = v2 - v;
  QDPIO::cout << "reconProdDir3Plus: || old - new || / || old || = " << sqrt( norm2(diff_v)/norm2(v2)) << endl;



#if 0

  gaussian(hferm_lat);
  hferm2_lat = u_lat*hferm_lat;
  ferm_lat = spinReconstructDir0Minus(hferm2_lat);

  ferm2_lat= spinReconstructDir0Minus(u_lat*hferm_lat);
  diff_lat = ferm_lat - ferm2_lat;
  QDPIO::cout << "reconProdDir0Minus: || old - new || / || old || = " << sqrt( norm2(diff_lat)/norm2(ferm_lat)) << endl;


  hferm2_lat = u_lat*hferm_lat;
  ferm_lat = spinReconstructDir0Minus(hferm2_lat);

  ferm2_lat= spinReconstructDir0Minus(u_lat*hferm_lat);
  diff_lat = ferm_lat - ferm2_lat;
  QDPIO::cout << "reconProdDir0Minus: || old - new || / || old || = " << sqrt( norm2(diff_lat)/norm2(ferm_lat)) << endl;
  // Timed To bolt

#endif
  QDP_finalize();

  exit(0);
}
  
   
