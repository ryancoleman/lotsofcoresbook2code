// $Id: t_qdp.cc,v 1.24 2005-03-21 05:31:07 edwards Exp $
//
/*! \file
 *  \brief Silly little internal test code
 *
 */

#include <iostream>
#include <cstdio>

#include "qdp.h"
#include "examples.h"
#include "qdp_util.h"

using namespace QDP;


int main(int argc, char **argv)
{
  // Put the machine into a known state
  QDP_initialize(&argc, &argv);

  // Setup the layout
  const int foo[] = {2,2,2,2};
  multi1d<int> nrow(Nd);
  nrow = foo;  // Use only Nd elements
  Layout::setLattSize(nrow);
  Layout::create();

  // Open a file for some sample output
  XMLFileWriter xml("t_qdp.xml");
  push(xml,"t_qdp.xml");

  // Initialize the random number generator
  Seed seed;
  seed = 11;
  RNG::setrn(seed);
  QDPIO::cout << "After setrn" << std::endl;

  // Time to play...
  Real r1;
  r1 = 17.0;
  QDPIO::cout << "r1 after fill\n" << r1 << std::endl;

  for(int i=0; i < 10; ++i)
  {
    random(r1);
    QDPIO::cout << "r1 after random\n" << r1 << std::endl;
  }

  // Check the multi-dim arrays
  multi2d<Real> boo(2,3);
  boo = 0.0;
  QDPIO::cout << "Fill boo with 0\n";
  for(int j=0; j < 2; ++j)
    for(int i=0; i < 3; ++i)
      QDPIO::cout << boo(j,i) << std::endl;

  QDPIO::cout << "Fill boo with random\n";
  for(int j=0; j < 2; ++j)
    for(int i=0; i < 3; ++i)
      random(boo(j,i));

  QDPIO::cout << "Check boo filled with random\n";
  for(int j=0; j < 2; ++j)
    for(int i=0; i < 3; ++i)
      QDPIO::cout << boo(j,i) << std::endl;

  QDPIO::cout << "Test indexing of boo\n";
  for(int j=0; j < 2; ++j)
    for(int i=0; i < 3; ++i)
      QDPIO::cout << boo[j][i] << std::endl;

  // Check the multi-dim arrays
  multi3d<Real> goo(2,3,2);
  goo = 0.0;
  QDPIO::cout << "Fill goo with 0\n";
  for(int k=0; k < 2; ++k)
    for(int j=0; j < 3; ++j)
      for(int i=0; i < 2; ++i)
	QDPIO::cout << goo(k,j,i) << std::endl;

  QDPIO::cout << "Fill goo with random\n";
  for(int k=0; k < 2; ++k)
    for(int j=0; j < 3; ++j)
      for(int i=0; i < 2; ++i)
      {
	random(goo(k,j,i));
	QDPIO::cout << goo(k,j,i) << std::endl;
      }


  QDPIO::cout << "Check goo filled with random\n";
  for(int k=0; k < 2; ++k)
    for(int j=0; j < 3; ++j)
      for(int i=0; i < 2; ++i)
      {
	goo[k][j][i] = goo(k,j,i);
	QDPIO::cout << goo(k,j,i) << std::endl; 
      }
  
  QDPIO::cout << "Test indexing of goo\n";
  for(int k=0; k < 2; ++k)
    for(int j=0; j < 3; ++j)
      for(int i=0; i < 2; ++i)
	QDPIO::cout << goo[k][j][i] << std::endl;

  // Test out lattice fields
  LatticeColorMatrix b1,b2,b3;

  b1 = 1.0;
  QDPIO::cout << "b1 after fill\n" << std::endl;
  write(xml,"b1", b1);

  random(b1);
  QDPIO::cout << "b1 after random\n" << std::endl;
  write(xml,"b1", b1);

  random(b2);
  gaussian(b3);
  QDPIO::cout << "b3 after gaussian\n";
  push(xml,"test_stuff");
  write(xml,"b3",b3);
  write(xml,"b3",b3);
  pop(xml);
  
#if 0
  Double dsum;
  dsum = norm2(b1);
  QDPIO::cout << "dsum = " << dsum << std::endl;
  xml << "dsum = ";
  write(xml,"dsum", dsum);

  junk(xml,b3,b1,b2,all);
#endif

#if 1
  // Test comparisons and mask operations
  LatticeBoolean lbtmp1, lbtmp2;
  LatticeReal lftmp1, lftmp2;

  random(lftmp1);
  random(lftmp2);
  lbtmp1 = lftmp1 < lftmp2;
#endif

  //! Example of calling a plaquette routine
  /*! NOTE: the STL is *not* used to hold gauge fields */
  QDPIO::cout << "Initialize vector of latticegauge\n";
  multi1d<LatticeColorMatrix> u(Nd);
  QDPIO::cout << "After initialize vector of latticegauge\n";
  Double w_plaq, s_plaq, t_plaq, link;

  QDPIO::cout << "Start gaussian\n";
  for(int m=0; m < u.size(); ++m)
    gaussian(u[m]);

  MesPlq(u, w_plaq, s_plaq, t_plaq, link);
  QDPIO::cout << "w_plaq = " << w_plaq << std::endl;
  QDPIO::cout << "link = " << link << std::endl;

#if 1
  // Play with gamma matrices - they should be implemented...
  //! Test out propagators
  LatticePropagator quark_prop_1, quark_prop_2;
  gaussian(quark_prop_1);
  gaussian(quark_prop_2);

  int j_decay = Nd-1;
  int length = Layout::lattSize()[j_decay];
  multi1d< multi1d<Real> > meson_prop;
  multi1d<int> t_source(Nd);
  t_source = 0;

  mesons(quark_prop_1, quark_prop_2, meson_prop, t_source, j_decay);
  write(xml,"meson_prop",meson_prop);

#if 1
  multi1d< multi1d<Complex> > baryon_prop;
  write(xml,"quark_prop_1", quark_prop_1);
  baryon(quark_prop_1, baryon_prop, t_source, j_decay, 1);
  write(xml,"t_source", t_source);
  write(xml,"j_decay", j_decay);
  write(xml,"baryon_prop", baryon_prop);
#endif

#endif

  // More frivolity and gaiety.
  LatticeFermion psi, chi;
  random(psi);
  chi = zero;
  dslash(chi, u, psi, +1, 0);

  write(xml,"psi", psi);
  write(xml,"chi", chi);

  pop(xml);
  xml.close();

  // Time to bolt
  QDP_finalize();

  exit(0);
}
