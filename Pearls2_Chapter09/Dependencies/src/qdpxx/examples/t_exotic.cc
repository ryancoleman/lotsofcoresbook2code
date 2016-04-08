/*! \file
 *  \brief Test various exotic qdp routines
 */

#include <iostream>
#include <cstdio>

#include "qdp.h"
#include "examples.h"

using namespace QDP;

int main(int argc, char *argv[])
{
  // Put the machine into a known state
  QDP_initialize(&argc, &argv);

  // Setup the layout
  const int foo[] = {2,2,2,2};
  multi1d<int> nrow(Nd);
  nrow = foo;  // Use only Nd elements
  Layout::setLattSize(nrow);
  Layout::create();

  XMLFileWriter xml("t_exotic.xml");
  push(xml,"t_exotic");

  push(xml,"lattis");
  write(xml,"Nd", Nd);
  write(xml,"Nc", Nc);
  write(xml,"nrow", nrow);
  pop(xml);


  // Construct a random gauge transformation
  LatticeColorMatrix g;

  gaussian(g);
  taproj(g);
  expm12(g);
  reunit(g);

#if 0
  // Try out colorContract
  {
    LatticeColorMatrix a,b,c;
    gaussian(a);
    gaussian(b);
    gaussian(c);

    LatticeComplex lc1 = colorContract(a,b,c);

    push(xml,"color_contract_orig");
    write(xml,"lc1",lc1);
    pop(xml);
   
    // Do a random gauge transformation
    LatticeColorMatrix tmp;
    tmp = g * a * adj(g);  a = tmp;
    tmp = g * b * adj(g);  b = tmp;
    tmp = g * c * adj(g);  c = tmp;

    // Try colorcontract again
    LatticeComplex lc2 = colorContract(a,b,c);
    push(xml,"color_contract_gauge_transf");
    write(xml,"lc2", lc2);
    pop(xml);
  }
#endif

#if 1
  // Try out colorCrossProduct
  {
    LatticeColorVector a,b;
    gaussian(a);
    gaussian(b);

    LatticeColorVector lc1 = colorCrossProduct(a,b);

    push(xml,"color_cross_product_orig");
    write(xml,"lc1",lc1);
    pop(xml);
   
    // Do a random gauge transformation
    LatticeColorVector tmp;
    tmp = g * a;  a = tmp;
    tmp = g * b;  b = tmp;

    // Try colorCrossProduct again
    LatticeColorVector lc2 = adj(g) * colorCrossProduct(a,b);
    push(xml,"color_cross_product_gauge_transf");
    write(xml,"lc2", lc2);
    pop(xml);
  }
#endif

#if 1
  // Try out colorVectorContract
  {
    LatticeColorVector a,b;
    gaussian(a);
    gaussian(b);

    LatticeComplex lc1 = colorVectorContract(a,b);

    push(xml,"color_vector_contract_orig");
    write(xml,"lc1",lc1);
    pop(xml);
   
    // Do a random gauge transformation
    LatticeColorVector tmp;
    tmp = g * a;  a = tmp;
    tmp = g * b;  b = tmp;

    // Try colorVectorContract again
    LatticeComplex lc2 = colorVectorContract(a,b);
    push(xml,"color_vector_contract_gauge_transf");
    write(xml,"lc2", lc2);
    pop(xml);
  }
#endif

#if 0
  // Try out localInnerProduct
  // Not working.
  {
    LatticeColorVector a;
    LatticeFermion b;
    LatticeSpinVector c;
    gaussian(a);
    gaussian(b);

    LatticeSpinVector lc1 = localInnerProduct(a,b);

    push(xml,"local_inner_product_orig");
    write(xml,"lc1",lc1);
    pop(xml);
   
    // Do a random gauge transformation
    a = LatticeColorVector(g * a);
    b = LatticeFermion(g * b);

    // Try colorCrossProduct again
    LatticeSpinVector lc2 = localInnerProduct(a,b);
    push(xml,"local_inner_product_gauge_transf");
    write(xml,"lc2", lc2);
    pop(xml);
  }
#endif

#if 0
  // Try out chiralProject{Plus,Minus}
  {
    LatticeFermion psi, chi1, chi2;
    gaussian(psi);

    chi1 = 0.5*(psi + Gamma(Ns*Ns-1)*psi);
    chi2 = chiralProjectPlus(psi);
    QDPIO::cout << "|chi1|^2 = " << norm2(chi1) << endl
		<< "|chi2|^2 = " << norm2(chi2) << endl
		<< "|chi2 - chi1|^2 = " << norm2(chi2-chi1) << endl;

    chi1 = 0.5*(psi - Gamma(Ns*Ns-1)*psi);
    chi2 = chiralProjectMinus(psi);
    QDPIO::cout << "|chi1|^2 = " << norm2(chi1) << endl
		<< "|chi2|^2 = " << norm2(chi2) << endl
		<< "|chi2 - chi1|^2 = " << norm2(chi2-chi1) << endl;
  }
#endif

#if 1
  // Try out norm2 on arrays
  {
    int N = 5;
    multi1d<LatticeFermion> psi(N);
    for(int n=0; n < N; ++n)
      gaussian(psi[n]);

    Double dnorm1 = 0;
    for(int n=0; n < N; ++n)
      dnorm1 += norm2(psi[n],odd);

    Double dnorm2 = norm2(psi,odd);

    QDPIO::cout << "|dnorm1|^2 = " << dnorm1 << std::endl
		<< "|dnorm2|^2 = " << dnorm2 << std::endl;
  }
#endif

#if 1
  // Try out some other funky ops
  {
    LatticeFermion psi;
    gaussian(psi);

    // Should be okay
    bool bad1 = isnan(psi);
    QDPIO::cout << "Test isnan(LF) [should be false] = " << bad1 << std::endl;

    LatticeColorMatrix u;
    gaussian(u);

    // Should be okay
    bool bad2 = isnan(u + u);
    QDPIO::cout << "Test isnan(LCM) [should be false] = " << bad2 << std::endl;

    // Intentionally do bad things
    u = sqrt(Real(-1));
    bool bad3 = isnan(u);

    QDPIO::cout << "Test isnan(LCM) [should be true] = " << bad3 << std::endl;
  }
#endif

  pop(xml);
  xml.close();

  // Time to bolt
  QDP_finalize();

  exit(0);
}
