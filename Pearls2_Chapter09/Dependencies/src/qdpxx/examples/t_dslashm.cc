// $Id: t_dslashm.cc,v 1.20 2007-02-24 01:00:29 bjoo Exp $
/*! \file
 *  \brief Test the Wilson-Dirac operator (dslash)
 */

#include <iostream>
#include <cstdio>

#include "qdp.h"
#include "examples.h"

#include <sys/time.h>

using namespace QDP;


int main(int argc, char **argv)
{
  // Put the machine into a known state
  QDP_initialize(&argc, &argv);

  // Setup the layout
  const int foo[] = {4,2,2,2};
  multi1d<int> nrow(Nd);
  nrow = foo;  // Use only Nd elements
  Layout::setLattSize(nrow);
  Layout::create();

  //! Test out propagators
  multi1d<LatticeColorMatrix> u(Nd);
  for(int m=0; m < u.size(); ++m)
    gaussian(u[m]);

  LatticeFermion psi, chi;
  random(psi);
  chi = zero;

  int iter = 1000;

  {
    int isign = +1;
    int cb = 0;
    QDPIO::cout << "Applying D" << std::endl;
      
    clock_t myt1=clock();
    for(int i=0; i < iter; i++)
      dslash(chi, u, psi, isign, cb);
    clock_t myt2=clock();
      
    double mydt=(double)(myt2-myt1)/((double)(CLOCKS_PER_SEC));
    mydt=1.0e6*mydt/((double)(iter*(Layout::vol()/2)));
      
    QDPIO::cout << "cb = " << cb << " isign = " << isign << std::endl;
    QDPIO::cout << "The time per lattice point is "<< mydt << " micro sec" 
		<< " (" <<  (double)(1392.0f/mydt) << ") Mflops " << std::endl;
  }

  chi = zero;

  {
    int isign = +1;
    int cb = 0;
    QDPIO::cout << "Applying D" << std::endl;
      
    clock_t myt1=clock();
    for(int i=0; i < iter; i++)
      dslash2(chi, u, psi, isign, cb);
    clock_t myt2=clock();
      
    double mydt=(double)(myt2-myt1)/((double)(CLOCKS_PER_SEC));
    mydt=1.0e6*mydt/((double)(iter*(Layout::vol()/2)));
      
    QDPIO::cout << "cb = " << cb << " isign = " << isign << std::endl;
    QDPIO::cout << "The time per lattice point is "<< mydt << " micro sec" 
		<< " (" <<  (double)(1392.0f/mydt) << ") Mflops " << std::endl;
  }


#if 1
  XMLFileWriter xml("t_dslashm.xml");
  push(xml,"t_dslashm");
  write(xml,"Nd", Nd);
  write(xml,"Nc", Nc);
  write(xml,"Ns", Ns);
  write(xml,"nrow", nrow);
  write(xml,"psi", psi);
  write(xml,"chi", chi);
  pop(xml);
  xml.close();
#endif

  LatticeFermion chi2;
  random(psi);
  chi = zero;
  chi2 = zero;
  for(int isign=-1; isign < 2; isign+=2) { 
    for(int cb=0; cb<2; cb++) {
      int otherCB= cb == 0 ? 1 : 0;
      dslash(chi, u, psi, isign, cb);
      dslash2(chi2, u, psi, isign, cb);
      LatticeFermion diff;
      diff[rb[cb]]= chi2 - chi;
      QDPIO::cout << "isign="<<isign<<" cb=" << cb << " Diff = " << sqrt( norm2(diff,rb[cb]) / norm2(psi, rb[otherCB]))<< std::endl;
    }
  }

  // Time to bolt
  QDP_finalize();

  exit(0);
}
