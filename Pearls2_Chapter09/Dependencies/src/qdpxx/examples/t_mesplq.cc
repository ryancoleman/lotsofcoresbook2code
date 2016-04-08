// $Id: t_mesplq.cc,v 1.26 2008-02-20 21:27:58 bjoo Exp $

#include <iostream>
#include <cstdio>

#include "qdp.h"
#include "examples.h"

using namespace QDP;

int main(int argc, char *argv[])
{
  START_CODE();

  // Put the machine into a known state
  QDP_initialize(&argc, &argv);

  // Setup the layout
  const int foo[] = {4,2,2,4};
  multi1d<int> nrow(Nd);
  nrow = foo;  // Use only Nd elements
  Layout::setLattSize(nrow);
  Layout::create();

#ifdef QDP_USE_LIBXML2
  XMLFileWriter xml("t_mesplq.xml");
  push(xml,"t_mesplq");

  push(xml,"lattis");
  write(xml,"Nd",Nd);
  write(xml,"Nc",Nc);
  write(xml,"nrow",nrow);
  pop(xml);
#endif

  //! Example of calling a plaquette routine
  /*! NOTE: the STL is *not* used to hold gauge fields */
  multi1d<LatticeColorMatrix> u(Nd);
  Double w_plaq, s_plaq, t_plaq, link;

  QDPIO::cout << "Start gaussian\n";
  for(int m=0; m < u.size(); ++m)
    gaussian(u[m]);

  // Reunitarize the gauge field
  QDPIO::cout << "Start reunit\n";
  for(int m=0; m < u.size(); ++m)
    reunit(u[m]);

  // Try out the plaquette routine
  QDPIO::cout << "Start mesplq\n";
  MesPlq(u, w_plaq, s_plaq, t_plaq, link);
  QDPIO::cout << "w_plaq = " << w_plaq << std::endl;
  QDPIO::cout << "link = " << link << std::endl;

#ifdef QDP_USE_LIBXML2
  // Write out the results
  push(xml,"observables");
  write(xml,"w_plaq",w_plaq);
  write(xml,"link",link);
  pop(xml);

  pop(xml);
  xml.close();
#endif

  // 
  QDPIO::cout << "rb[0] has ordered rep=" << rb[0].hasOrderedRep() << std::endl;
  QDPIO::cout << "rb[1] has ordered rep=" << rb[1].hasOrderedRep() << std::endl;

  QDPIO::cout << "rb3[0] has ordered rep=" << rb3[0].hasOrderedRep() << std::endl;
  QDPIO::cout << "rb3[1] has ordered rep=" << rb3[1].hasOrderedRep() << std::endl;
  // Time to bolt
  QDP_finalize();

  END_CODE();
  exit(0);
}
