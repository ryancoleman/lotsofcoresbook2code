#include "qdp.h"
#include "testVaxpyDouble.h"
#include "unittest.h"

// #include "scalarsite_sse/sse_blas_vaxpy3_double.h"

using namespace QDP;
using namespace Assertions;

// Trinity of tests: Check SSE against Handrolled
//                   Check Handrolled against QDP++
//                   Check QDP++ against SSE


// Test 1. Check hand rolled against -- 'optimized'
void
testVaxpy4_1::run()
{
  Double a=Double(2.3);
  LatticeDiracFermionD3 x;
  LatticeDiracFermionD3 y;
  LatticeDiracFermionD3 z1;
  LatticeDiracFermionD3 z2;

  gaussian(x);
  gaussian(y);

  z1=y;
  z2=y;

  // Loop over sites
  for(int i=all.start(); i <= all.end(); i++) { 
    // Loop over spins
    for(int spin=0; spin < 4; spin++) {
      // Loop over colors (case 3 only) 
      for(int col=0; col < 3; col++) { 

	z1.elem(i).elem(spin).elem(col).real() = 
	  a.elem().elem().elem().elem()*x.elem(i).elem(spin).elem(col).real()
	  + z1.elem(i).elem(spin).elem(col).real();

	z1.elem(i).elem(spin).elem(col).imag() = 
	  a.elem().elem().elem().elem()*x.elem(i).elem(spin).elem(col).imag()
	  + z1.elem(i).elem(spin).elem(col).imag();
      }
    }
  }

  int n_4vec = (all.end() - all.start() + 1);
  REAL64* xptr = (REAL64 *)&(x.elem(all.start()).elem(0).elem(0).real());
  REAL64* zptr = &(z2.elem(all.start()).elem(0).elem(0).real());
  REAL64 ar = a.elem().elem().elem().elem();
  REAL64* aptr = &ar;

  vaxpy4(zptr, aptr, xptr, n_4vec);

  for(int i=all.start(); i <= all.end(); i++) { 
    // Loop over spins
    for(int spin=0; spin < 4; spin++) {
      // Loop over colors (case 3 only) 
      for(int col=0; col < 3; col++) { 

	double realdiff = z1.elem(i).elem(spin).elem(col).real()
	  - z2.elem(i).elem(spin).elem(col).real();

	assertion( fabs(realdiff) < 1.0e-14 );

	double imagdiff = z1.elem(i).elem(spin).elem(col).imag()
	  - z2.elem(i).elem(spin).elem(col).imag();

	assertion( fabs(imagdiff) < 1.0e-14 );

      }
    }
  }
}

// Test 2. Test Hand rolled against QDP++ expression
void
testVaxpy4_2::run()
{
  Double a=Double(2.3);
  LatticeDiracFermionD3 x;
  LatticeDiracFermionD3 y;
  LatticeDiracFermionD3 z1;
  LatticeDiracFermionD3 z2;

  gaussian(x);
  gaussian(y);

  z1=y;
  z2=y;

  // Loop over sites
  for(int i=all.start(); i <= all.end(); i++) { 
    // Loop over spins
    for(int spin=0; spin < 4; spin++) {
      // Loop over colors (case 3 only) 
      for(int col=0; col < 3; col++) { 

	z1.elem(i).elem(spin).elem(col).real() = 
	  a.elem().elem().elem().elem()*x.elem(i).elem(spin).elem(col).real()
	  + z1.elem(i).elem(spin).elem(col).real();

	z1.elem(i).elem(spin).elem(col).imag() = 
	  a.elem().elem().elem().elem()*x.elem(i).elem(spin).elem(col).imag()
	  + z1.elem(i).elem(spin).elem(col).imag();
      }
    }
  }

  z2 = a*x + z2;

  for(int i=all.start(); i <= all.end(); i++) { 
    // Loop over spins
    for(int spin=0; spin < 4; spin++) {
      // Loop over colors (case 3 only) 
      for(int col=0; col < 3; col++) { 

	double realdiff = z1.elem(i).elem(spin).elem(col).real()
	  - z2.elem(i).elem(spin).elem(col).real();

	assertion( fabs(realdiff) < 1.0e-14 );

	double imagdiff = z1.elem(i).elem(spin).elem(col).imag()
	  - z2.elem(i).elem(spin).elem(col).imag();

	assertion( fabs(imagdiff) < 1.0e-14 );

      }
    }
  }
  
}


// Test 3. Check QDP++ against 'optimized'
void
testVaxpy4_3::run()
{
  Double a=Double(2.3);
  LatticeDiracFermionD3 x;
  LatticeDiracFermionD3 y;
  LatticeDiracFermionD3 z1;
  LatticeDiracFermionD3 z2;

  gaussian(x);
  gaussian(y);
  z1=y;
  z2=y;
  // QDP++ 
  z1 = a*x + z1;

  // Optimized
  int n_4vec = (all.end() - all.start() + 1);
  REAL64* xptr = (REAL64 *)&(x.elem(all.start()).elem(0).elem(0).real());
  REAL64* zptr = &(z2.elem(all.start()).elem(0).elem(0).real());
  REAL64 ar = a.elem().elem().elem().elem();
  REAL64* aptr = &ar;

  vaxpy4(zptr, aptr, xptr, n_4vec);

  for(int i=all.start(); i <= all.end(); i++) { 
    // Loop over spins
    for(int spin=0; spin < 4; spin++) {
      // Loop over colors (case 3 only) 
      for(int col=0; col < 3; col++) { 

	double realdiff = z1.elem(i).elem(spin).elem(col).real()
	  - z2.elem(i).elem(spin).elem(col).real();

	assertion( fabs(realdiff) < 1.0e-14 );

	double imagdiff = z1.elem(i).elem(spin).elem(col).imag()
	  - z2.elem(i).elem(spin).elem(col).imag();

	assertion( fabs(imagdiff) < 1.0e-14 );

      }
    }
  }
}



// Test 4. Check hand rolled against -- 'optimized'
// On cb2 subset
void
testVaxpy4_RB0_1::run()
{
  Double a=Double(2.3);
  LatticeDiracFermionD3 x;
  LatticeDiracFermionD3 y;
  LatticeDiracFermionD3 z1;
  LatticeDiracFermionD3 z2;

  gaussian(x);
  gaussian(y);

  // Zero (the other cb)
  z1=zero;
  z2=zero;

  const int *tab = rb[0].siteTable().slice();

  // By hand do the loop over the sites 
  for(int j=0; j < rb[0].numSiteTable(); j++) {
    int i = tab[j];

    for(int spin=0; spin < 4; spin++) {
      // Loop over colors (case 3 only) 
      for(int col=0; col < 3; col++) { 

	z1.elem(i).elem(spin).elem(col).real() = 
	  a.elem().elem().elem().elem()*x.elem(i).elem(spin).elem(col).real()
	  + y.elem(i).elem(spin).elem(col).real();

	z1.elem(i).elem(spin).elem(col).imag() = 
	  a.elem().elem().elem().elem()*x.elem(i).elem(spin).elem(col).imag()
	  + y.elem(i).elem(spin).elem(col).imag();
      }
    }
  }

  if( rb[0].hasOrderedRep() ) { 
    int n_4vec = (rb[0].end() - rb[0].start() + 1);
    
    REAL64* xptr = (REAL64 *)&(x.elem(rb[0].start()).elem(0).elem(0).real());
    REAL64* yptr = &(y.elem(rb[0].start()).elem(0).elem(0).real());
    REAL64* zptr = &(z2.elem(rb[0].start()).elem(0).elem(0).real());
    REAL64 ar = a.elem().elem().elem().elem();
    REAL64* aptr = &ar;

    vaxpyz4(zptr, aptr, xptr, yptr, n_4vec);
  }
  else { 
    const int *tab = rb[0].siteTable().slice();
    for(int j=0; j < rb[0].numSiteTable(); j++) {
      int i = tab[j];
      REAL64* xptr = (REAL64 *)&(x.elem(i).elem(0).elem(0).real());
      REAL64* yptr = &(y.elem(i).elem(0).elem(0).real());
      REAL64* zptr = &(z2.elem(i).elem(0).elem(0).real());
      REAL64 ar = a.elem().elem().elem().elem();
      REAL64* aptr = &ar;
      vaxpyz4(zptr, aptr, xptr, yptr, 1);
    }
  }

  for(int i=all.start(); i <= all.end(); i++) { 
    // Loop over spins
    for(int spin=0; spin < 4; spin++) {
      // Loop over colors (case 3 only) 
      for(int col=0; col < 3; col++) { 

	double realdiff = z1.elem(i).elem(spin).elem(col).real()
	  - z2.elem(i).elem(spin).elem(col).real();

	assertion( fabs(realdiff) < 1.0e-14 );

	double imagdiff = z1.elem(i).elem(spin).elem(col).imag()
	  - z2.elem(i).elem(spin).elem(col).imag();

	assertion( fabs(imagdiff) < 1.0e-14 );

      }
    }
  }
}

// Test 5. Test Hand rolled against QDP++ expression
void
testVaxpy4_RB0_2::run()
{
  Double a=Double(2.3);
  LatticeDiracFermionD3 x;
  LatticeDiracFermionD3 y;
  LatticeDiracFermionD3 z1;
  LatticeDiracFermionD3 z2;

  gaussian(x);
  gaussian(y);

  // Zero (the other cb)
  z1=zero;
  z2=zero;

  const int *tab = rb[0].siteTable().slice();

  // By hand do the loop over the sites 
  for(int j=0; j < rb[0].numSiteTable(); j++) {
    int i = tab[j];

    for(int spin=0; spin < 4; spin++) {
      // Loop over colors (case 3 only) 
      for(int col=0; col < 3; col++) { 

	z1.elem(i).elem(spin).elem(col).real() = 
	  a.elem().elem().elem().elem()*x.elem(i).elem(spin).elem(col).real()
	  + y.elem(i).elem(spin).elem(col).real();

	z1.elem(i).elem(spin).elem(col).imag() = 
	  a.elem().elem().elem().elem()*x.elem(i).elem(spin).elem(col).imag()
	  + y.elem(i).elem(spin).elem(col).imag();
      }
    }
  }

  z2[rb[0]] = a*x + y;

  if( rb[0].hasOrderedRep() ) { 
    int n_4vec = (rb[0].end() - rb[0].start() + 1);
    
    REAL64* xptr = (REAL64 *)&(x.elem(rb[0].start()).elem(0).elem(0).real());
    REAL64* yptr = &(y.elem(rb[0].start()).elem(0).elem(0).real());
    REAL64* zptr = &(z2.elem(rb[0].start()).elem(0).elem(0).real());
    REAL64 ar = a.elem().elem().elem().elem();
    REAL64* aptr = &ar;

    vaxpyz4(zptr, aptr, xptr, yptr, n_4vec);
  }
  else { 
    const int *tab = rb[0].siteTable().slice();
    for(int j=0; j < rb[0].numSiteTable(); j++) {
      int i = tab[j];
      REAL64* xptr = (REAL64 *)&(x.elem(i).elem(0).elem(0).real());
      REAL64* yptr = &(y.elem(i).elem(0).elem(0).real());
      REAL64* zptr = &(z2.elem(i).elem(0).elem(0).real());
      REAL64 ar = a.elem().elem().elem().elem();
      REAL64* aptr = &ar;
      vaxpyz4(zptr, aptr, xptr, yptr, 1);
    }
  }

  for(int i=all.start(); i <= all.end(); i++) { 
    // Loop over spins
    for(int spin=0; spin < 4; spin++) {
      // Loop over colors (case 3 only) 
      for(int col=0; col < 3; col++) { 

	double realdiff = z1.elem(i).elem(spin).elem(col).real()
	  - z2.elem(i).elem(spin).elem(col).real();

	assertion( fabs(realdiff) < 1.0e-14 );

	double imagdiff = z1.elem(i).elem(spin).elem(col).imag()
	  - z2.elem(i).elem(spin).elem(col).imag();

	assertion( fabs(imagdiff) < 1.0e-14 );

      }
    }
  }
  
}


// Test 6. Check QDP++ against 'optimized'
void
testVaxpy4_RB0_3::run()
{
  Double a=Double(2.3);
  LatticeDiracFermionD3 x;
  LatticeDiracFermionD3 y;
  LatticeDiracFermionD3 z1;
  LatticeDiracFermionD3 z2;

  gaussian(x);
  gaussian(y);

  // QDP++ 
  z1[rb[0]] = a*x + y;

  // Optimized
 if( rb[0].hasOrderedRep() ) { 
    int n_4vec = (rb[0].end() - rb[0].start() + 1);
    
    REAL64* xptr = (REAL64 *)&(x.elem(rb[0].start()).elem(0).elem(0).real());
    REAL64* yptr = &(y.elem(rb[0].start()).elem(0).elem(0).real());
    REAL64* zptr = &(z2.elem(rb[0].start()).elem(0).elem(0).real());
    REAL64 ar = a.elem().elem().elem().elem();
    REAL64* aptr = &ar;

    vaxpyz4(zptr, aptr, xptr, yptr, n_4vec);
  }
  else { 
    const int *tab = rb[0].siteTable().slice();
    for(int j=0; j < rb[0].numSiteTable(); j++) {
      int i = tab[j];
      REAL64* xptr = (REAL64 *)&(x.elem(i).elem(0).elem(0).real());
      REAL64* yptr = &(y.elem(i).elem(0).elem(0).real());
      REAL64* zptr = &(z2.elem(i).elem(0).elem(0).real());
      REAL64 ar = a.elem().elem().elem().elem();
      REAL64* aptr = &ar;
      vaxpyz4(zptr, aptr, xptr, yptr, 1);
    }
  }

  for(int i=all.start(); i <= all.end(); i++) { 
    // Loop over spins
    for(int spin=0; spin < 4; spin++) {
      // Loop over colors (case 3 only) 
      for(int col=0; col < 3; col++) { 

	double realdiff = z1.elem(i).elem(spin).elem(col).real()
	  - z2.elem(i).elem(spin).elem(col).real();

	assertion( fabs(realdiff) < 1.0e-14 );

	double imagdiff = z1.elem(i).elem(spin).elem(col).imag()
	  - z2.elem(i).elem(spin).elem(col).imag();

	assertion( fabs(imagdiff) < 1.0e-14 );

      }
    }
  }
}

// Test 7. Check hand rolled against -- 'optimized'
// On cb2 subset
void
testVaxpy4_RB1_1::run()
{
  Double a=Double(2.3);
  LatticeDiracFermionD3 x;
  LatticeDiracFermionD3 y;
  LatticeDiracFermionD3 z1;
  LatticeDiracFermionD3 z2;

  gaussian(x);
  gaussian(y);

  // Zero (the other cb)
  z1=zero;
  z2=zero;

  const int *tab = rb[1].siteTable().slice();

  // By hand do the loop over the sites 
  for(int j=0; j < rb[1].numSiteTable(); j++) {
    int i = tab[j];

    for(int spin=0; spin < 4; spin++) {
      // Loop over colors (case 3 only) 
      for(int col=0; col < 3; col++) { 

	z1.elem(i).elem(spin).elem(col).real() = 
	  a.elem().elem().elem().elem()*x.elem(i).elem(spin).elem(col).real()
	  + y.elem(i).elem(spin).elem(col).real();

	z1.elem(i).elem(spin).elem(col).imag() = 
	  a.elem().elem().elem().elem()*x.elem(i).elem(spin).elem(col).imag()
	  + y.elem(i).elem(spin).elem(col).imag();
      }
    }
  }

  if( rb[1].hasOrderedRep() ) { 
    int n_4vec = (rb[1].end() - rb[1].start() + 1);
    
    REAL64* xptr = (REAL64 *)&(x.elem(rb[1].start()).elem(0).elem(0).real());
    REAL64* yptr = &(y.elem(rb[1].start()).elem(0).elem(0).real());
    REAL64* zptr = &(z2.elem(rb[1].start()).elem(0).elem(0).real());
    REAL64 ar = a.elem().elem().elem().elem();
    REAL64* aptr = &ar;

    vaxpyz4(zptr, aptr, xptr, yptr, n_4vec);
  }
  else { 
    const int *tab = rb[1].siteTable().slice();
    for(int j=0; j < rb[1].numSiteTable(); j++) {
      int i = tab[j];
      REAL64* xptr = (REAL64 *)&(x.elem(i).elem(0).elem(0).real());
      REAL64* yptr = &(y.elem(i).elem(0).elem(0).real());
      REAL64* zptr = &(z2.elem(i).elem(0).elem(0).real());
      REAL64 ar = a.elem().elem().elem().elem();
      REAL64* aptr = &ar;
      vaxpyz4(zptr, aptr, xptr, yptr, 1);
    }
  }

  for(int i=all.start(); i <= all.end(); i++) { 
    // Loop over spins
    for(int spin=0; spin < 4; spin++) {
      // Loop over colors (case 3 only) 
      for(int col=0; col < 3; col++) { 

	double realdiff = z1.elem(i).elem(spin).elem(col).real()
	  - z2.elem(i).elem(spin).elem(col).real();

	assertion( fabs(realdiff) < 1.0e-14 );

	double imagdiff = z1.elem(i).elem(spin).elem(col).imag()
	  - z2.elem(i).elem(spin).elem(col).imag();

	assertion( fabs(imagdiff) < 1.0e-14 );

      }
    }
  }
}

// Test 8. Test Hand rolled against QDP++ expression
void
testVaxpy4_RB1_2::run()
{
  Double a=Double(2.3);
  LatticeDiracFermionD3 x;
  LatticeDiracFermionD3 y;
  LatticeDiracFermionD3 z1;
  LatticeDiracFermionD3 z2;

  gaussian(x);
  gaussian(y);

  // Zero (the other cb)
  z1=zero;
  z2=zero;

  const int *tab = rb[1].siteTable().slice();

  // By hand do the loop over the sites 
  for(int j=0; j < rb[1].numSiteTable(); j++) {
    int i = tab[j];

    for(int spin=0; spin < 4; spin++) {
      // Loop over colors (case 3 only) 
      for(int col=0; col < 3; col++) { 

	z1.elem(i).elem(spin).elem(col).real() = 
	  a.elem().elem().elem().elem()*x.elem(i).elem(spin).elem(col).real()
	  + y.elem(i).elem(spin).elem(col).real();

	z1.elem(i).elem(spin).elem(col).imag() = 
	  a.elem().elem().elem().elem()*x.elem(i).elem(spin).elem(col).imag()
	  + y.elem(i).elem(spin).elem(col).imag();
      }
    }
  }

  z2[rb[1]] = a*x + y;

  if( rb[1].hasOrderedRep() ) { 
    int n_4vec = (rb[1].end() - rb[1].start() + 1);
    
    REAL64* xptr = (REAL64 *)&(x.elem(rb[1].start()).elem(0).elem(0).real());
    REAL64* yptr = &(y.elem(rb[1].start()).elem(0).elem(0).real());
    REAL64* zptr = &(z2.elem(rb[1].start()).elem(0).elem(0).real());
    REAL64 ar = a.elem().elem().elem().elem();
    REAL64* aptr = &ar;

    vaxpyz4(zptr, aptr, xptr, yptr, n_4vec);
  }
  else { 
    const int *tab = rb[1].siteTable().slice();
    for(int j=0; j < rb[1].numSiteTable(); j++) {
      int i = tab[j];
      REAL64* xptr = (REAL64 *)&(x.elem(i).elem(0).elem(0).real());
      REAL64* yptr = &(y.elem(i).elem(0).elem(0).real());
      REAL64* zptr = &(z2.elem(i).elem(0).elem(0).real());
      REAL64 ar = a.elem().elem().elem().elem();
      REAL64* aptr = &ar;
      vaxpyz4(zptr, aptr, xptr, yptr, 1);
    }
  }

  for(int i=all.start(); i <= all.end(); i++) { 
    // Loop over spins
    for(int spin=0; spin < 4; spin++) {
      // Loop over colors (case 3 only) 
      for(int col=0; col < 3; col++) { 

	double realdiff = z1.elem(i).elem(spin).elem(col).real()
	  - z2.elem(i).elem(spin).elem(col).real();

	assertion( fabs(realdiff) < 1.0e-14 );

	double imagdiff = z1.elem(i).elem(spin).elem(col).imag()
	  - z2.elem(i).elem(spin).elem(col).imag();

	assertion( fabs(imagdiff) < 1.0e-14 );

      }
    }
  }
  
}


// Test 9. Check QDP++ against 'optimized'
void
testVaxpy4_RB1_3::run()
{
  Double a=Double(2.3);
  LatticeDiracFermionD3 x;
  LatticeDiracFermionD3 y;
  LatticeDiracFermionD3 z1;
  LatticeDiracFermionD3 z2;

  gaussian(x);
  gaussian(y);

  // QDP++ 
  z1[rb[1]] = a*x + y;

  // Optimized
 if( rb[1].hasOrderedRep() ) { 
    int n_4vec = (rb[1].end() - rb[1].start() + 1);
    
    REAL64* xptr = (REAL64 *)&(x.elem(rb[1].start()).elem(0).elem(0).real());
    REAL64* yptr = &(y.elem(rb[1].start()).elem(0).elem(0).real());
    REAL64* zptr = &(z2.elem(rb[1].start()).elem(0).elem(0).real());
    REAL64 ar = a.elem().elem().elem().elem();
    REAL64* aptr = &ar;

    vaxpyz4(zptr, aptr, xptr, yptr, n_4vec);
  }
  else { 
    const int *tab = rb[1].siteTable().slice();
    for(int j=0; j < rb[1].numSiteTable(); j++) {
      int i = tab[j];
      REAL64* xptr = (REAL64 *)&(x.elem(i).elem(0).elem(0).real());
      REAL64* yptr = &(y.elem(i).elem(0).elem(0).real());
      REAL64* zptr = &(z2.elem(i).elem(0).elem(0).real());
      REAL64 ar = a.elem().elem().elem().elem();
      REAL64* aptr = &ar;
      vaxpyz4(zptr, aptr, xptr, yptr, 1);
    }
  }

  for(int i=all.start(); i <= all.end(); i++) { 
    // Loop over spins
    for(int spin=0; spin < 4; spin++) {
      // Loop over colors (case 3 only) 
      for(int col=0; col < 3; col++) { 

	double realdiff = z1.elem(i).elem(spin).elem(col).real()
	  - z2.elem(i).elem(spin).elem(col).real();

	assertion( fabs(realdiff) < 1.0e-14 );

	double imagdiff = z1.elem(i).elem(spin).elem(col).imag()
	  - z2.elem(i).elem(spin).elem(col).imag();

	assertion( fabs(imagdiff) < 1.0e-14 );

      }
    }
  }
}

// Test 10. Check hand rolled against -- 'optimized'
// On cb3d subset
void
testVaxpy4_RB31_1::run()
{
  Double a=Double(2.3);
  LatticeDiracFermionD3 x;
  LatticeDiracFermionD3 y;
  LatticeDiracFermionD3 z1;
  LatticeDiracFermionD3 z2;

  gaussian(x);
  gaussian(y);

  // Zero (the other cb)
  z1=zero;
  z2=zero;

  const int *tab = rb3[1].siteTable().slice();

  // By hand do the loop over the sites 
  for(int j=0; j < rb3[1].numSiteTable(); j++) {
    int i = tab[j];

    for(int spin=0; spin < 4; spin++) {
      // Loop over colors (case 3 only) 
      for(int col=0; col < 3; col++) { 

	z1.elem(i).elem(spin).elem(col).real() = 
	  a.elem().elem().elem().elem()*x.elem(i).elem(spin).elem(col).real()
	  + y.elem(i).elem(spin).elem(col).real();

	z1.elem(i).elem(spin).elem(col).imag() = 
	  a.elem().elem().elem().elem()*x.elem(i).elem(spin).elem(col).imag()
	  + y.elem(i).elem(spin).elem(col).imag();
      }
    }
  }

  if( rb3[1].hasOrderedRep() ) { 
    int n_4vec = (rb3[1].end() - rb3[1].start() + 1);
    
    REAL64* xptr = (REAL64 *)&(x.elem(rb3[1].start()).elem(0).elem(0).real());
    REAL64* yptr = &(y.elem(rb3[1].start()).elem(0).elem(0).real());
    REAL64* zptr = &(z2.elem(rb3[1].start()).elem(0).elem(0).real());
    REAL64 ar = a.elem().elem().elem().elem();
    REAL64* aptr = &ar;

    vaxpyz4(zptr, aptr, xptr, yptr, n_4vec);
  }
  else { 
    const int *tab = rb3[1].siteTable().slice();
    for(int j=0; j < rb3[1].numSiteTable(); j++) {
      int i = tab[j];
      REAL64* xptr = (REAL64 *)&(x.elem(i).elem(0).elem(0).real());
      REAL64* yptr = &(y.elem(i).elem(0).elem(0).real());
      REAL64* zptr = &(z2.elem(i).elem(0).elem(0).real());
      REAL64 ar = a.elem().elem().elem().elem();
      REAL64* aptr = &ar;
      vaxpyz4(zptr, aptr, xptr, yptr, 1);
    }
  }

  for(int i=all.start(); i <= all.end(); i++) { 
    // Loop over spins
    for(int spin=0; spin < 4; spin++) {
      // Loop over colors (case 3 only) 
      for(int col=0; col < 3; col++) { 

	double realdiff = z1.elem(i).elem(spin).elem(col).real()
	  - z2.elem(i).elem(spin).elem(col).real();

	assertion( fabs(realdiff) < 1.0e-14 );

	double imagdiff = z1.elem(i).elem(spin).elem(col).imag()
	  - z2.elem(i).elem(spin).elem(col).imag();

	assertion( fabs(imagdiff) < 1.0e-14 );

      }
    }
  }
}

// Test 11. Test Hand rolled against QDP++ expression
void
testVaxpy4_RB31_2::run()
{
  Double a=Double(2.3);
  LatticeDiracFermionD3 x;
  LatticeDiracFermionD3 y;
  LatticeDiracFermionD3 z1;
  LatticeDiracFermionD3 z2;

  gaussian(x);
  gaussian(y);

  // Zero (the other cb)
  z1=zero;
  z2=zero;

  const int *tab = rb3[1].siteTable().slice();

  // By hand do the loop over the sites 
  for(int j=0; j < rb3[1].numSiteTable(); j++) {
    int i = tab[j];

    for(int spin=0; spin < 4; spin++) {
      // Loop over colors (case 3 only) 
      for(int col=0; col < 3; col++) { 

	z1.elem(i).elem(spin).elem(col).real() = 
	  a.elem().elem().elem().elem()*x.elem(i).elem(spin).elem(col).real()
	  + y.elem(i).elem(spin).elem(col).real();

	z1.elem(i).elem(spin).elem(col).imag() = 
	  a.elem().elem().elem().elem()*x.elem(i).elem(spin).elem(col).imag()
	  + y.elem(i).elem(spin).elem(col).imag();
      }
    }
  }

  z2[rb3[1]] = a*x + y;

  if( rb3[1].hasOrderedRep() ) { 
    int n_4vec = (rb3[1].end() - rb3[1].start() + 1);
    
    REAL64* xptr = (REAL64 *)&(x.elem(rb3[1].start()).elem(0).elem(0).real());
    REAL64* yptr = &(y.elem(rb3[1].start()).elem(0).elem(0).real());
    REAL64* zptr = &(z2.elem(rb3[1].start()).elem(0).elem(0).real());
    REAL64 ar = a.elem().elem().elem().elem();
    REAL64* aptr = &ar;

    vaxpyz4(zptr, aptr, xptr, yptr, n_4vec);
  }
  else { 
    const int *tab = rb3[1].siteTable().slice();
    for(int j=0; j < rb3[1].numSiteTable(); j++) {
      int i = tab[j];
      REAL64* xptr = (REAL64 *)&(x.elem(i).elem(0).elem(0).real());
      REAL64* yptr = &(y.elem(i).elem(0).elem(0).real());
      REAL64* zptr = &(z2.elem(i).elem(0).elem(0).real());
      REAL64 ar = a.elem().elem().elem().elem();
      REAL64* aptr = &ar;
      vaxpyz4(zptr, aptr, xptr, yptr, 1);
    }
  }

  for(int i=all.start(); i <= all.end(); i++) { 
    // Loop over spins
    for(int spin=0; spin < 4; spin++) {
      // Loop over colors (case 3 only) 
      for(int col=0; col < 3; col++) { 

	double realdiff = z1.elem(i).elem(spin).elem(col).real()
	  - z2.elem(i).elem(spin).elem(col).real();

	assertion( fabs(realdiff) < 1.0e-14 );

	double imagdiff = z1.elem(i).elem(spin).elem(col).imag()
	  - z2.elem(i).elem(spin).elem(col).imag();

	assertion( fabs(imagdiff) < 1.0e-14 );

      }
    }
  }
  
}


// Test 12. Check QDP++ against 'optimized'
void
testVaxpy4_RB31_3::run()
{
  Double a=Double(2.3);
  LatticeDiracFermionD3 x;
  LatticeDiracFermionD3 y;
  LatticeDiracFermionD3 z1;
  LatticeDiracFermionD3 z2;

  gaussian(x);
  gaussian(y);

  // QDP++ 
  z1[rb3[1]] = a*x + y;

  // Optimized
 if( rb3[1].hasOrderedRep() ) { 
    int n_4vec = (rb3[1].end() - rb3[1].start() + 1);
    
    REAL64* xptr = (REAL64 *)&(x.elem(rb3[1].start()).elem(0).elem(0).real());
    REAL64* yptr = &(y.elem(rb3[1].start()).elem(0).elem(0).real());
    REAL64* zptr = &(z2.elem(rb3[1].start()).elem(0).elem(0).real());
    REAL64 ar = a.elem().elem().elem().elem();
    REAL64* aptr = &ar;

    vaxpyz4(zptr, aptr, xptr, yptr, n_4vec);
  }
  else { 
    const int *tab = rb3[1].siteTable().slice();
    for(int j=0; j < rb3[1].numSiteTable(); j++) {
      int i = tab[j];
      REAL64* xptr = (REAL64 *)&(x.elem(i).elem(0).elem(0).real());
      REAL64* yptr = &(y.elem(i).elem(0).elem(0).real());
      REAL64* zptr = &(z2.elem(i).elem(0).elem(0).real());
      REAL64 ar = a.elem().elem().elem().elem();
      REAL64* aptr = &ar;
      vaxpyz4(zptr, aptr, xptr, yptr, 1);
    }
  }

  for(int i=all.start(); i <= all.end(); i++) { 
    // Loop over spins
    for(int spin=0; spin < 4; spin++) {
      // Loop over colors (case 3 only) 
      for(int col=0; col < 3; col++) { 

	double realdiff = z1.elem(i).elem(spin).elem(col).real()
	  - z2.elem(i).elem(spin).elem(col).real();

	assertion( fabs(realdiff) < 1.0e-14 );

	double imagdiff = z1.elem(i).elem(spin).elem(col).imag()
	  - z2.elem(i).elem(spin).elem(col).imag();

	assertion( fabs(imagdiff) < 1.0e-14 );

      }
    }
  }
}


// Test 13. Check hand rolled against -- 'optimized'
// On cb3d subset
void
testVaxpy4_RB30_1::run()
{
  Double a=Double(2.3);
  LatticeDiracFermionD3 x;
  LatticeDiracFermionD3 y;
  LatticeDiracFermionD3 z1;
  LatticeDiracFermionD3 z2;

  gaussian(x);
  gaussian(y);

  // Zero (the other cb)
  z1=zero;
  z2=zero;

  const int *tab = rb3[0].siteTable().slice();

  // By hand do the loop over the sites 
  for(int j=0; j < rb3[0].numSiteTable(); j++) {
    int i = tab[j];

    for(int spin=0; spin < 4; spin++) {
      // Loop over colors (case 3 only) 
      for(int col=0; col < 3; col++) { 

	z1.elem(i).elem(spin).elem(col).real() = 
	  a.elem().elem().elem().elem()*x.elem(i).elem(spin).elem(col).real()
	  + y.elem(i).elem(spin).elem(col).real();

	z1.elem(i).elem(spin).elem(col).imag() = 
	  a.elem().elem().elem().elem()*x.elem(i).elem(spin).elem(col).imag()
	  + y.elem(i).elem(spin).elem(col).imag();
      }
    }
  }

  if( rb3[0].hasOrderedRep() ) { 
    int n_4vec = (rb3[0].end() - rb3[0].start() + 1);
    
    REAL64* xptr = (REAL64 *)&(x.elem(rb3[0].start()).elem(0).elem(0).real());
    REAL64* yptr = &(y.elem(rb3[0].start()).elem(0).elem(0).real());
    REAL64* zptr = &(z2.elem(rb3[0].start()).elem(0).elem(0).real());
    REAL64 ar = a.elem().elem().elem().elem();
    REAL64* aptr = &ar;

    vaxpyz4(zptr, aptr, xptr, yptr, n_4vec);
  }
  else { 
    const int *tab = rb3[0].siteTable().slice();
    for(int j=0; j < rb3[0].numSiteTable(); j++) {
      int i = tab[j];
      REAL64* xptr = (REAL64 *)&(x.elem(i).elem(0).elem(0).real());
      REAL64* yptr = &(y.elem(i).elem(0).elem(0).real());
      REAL64* zptr = &(z2.elem(i).elem(0).elem(0).real());
      REAL64 ar = a.elem().elem().elem().elem();
      REAL64* aptr = &ar;
      vaxpyz4(zptr, aptr, xptr, yptr, 1);
    }
  }

  for(int i=all.start(); i <= all.end(); i++) { 
    // Loop over spins
    for(int spin=0; spin < 4; spin++) {
      // Loop over colors (case 3 only) 
      for(int col=0; col < 3; col++) { 

	double realdiff = z1.elem(i).elem(spin).elem(col).real()
	  - z2.elem(i).elem(spin).elem(col).real();

	assertion( fabs(realdiff) < 1.0e-14 );

	double imagdiff = z1.elem(i).elem(spin).elem(col).imag()
	  - z2.elem(i).elem(spin).elem(col).imag();

	assertion( fabs(imagdiff) < 1.0e-14 );

      }
    }
  }
}

// Test 14. Test Hand rolled against QDP++ expression
void
testVaxpy4_RB30_2::run()
{
  Double a=Double(2.3);
  LatticeDiracFermionD3 x;
  LatticeDiracFermionD3 y;
  LatticeDiracFermionD3 z1;
  LatticeDiracFermionD3 z2;

  gaussian(x);
  gaussian(y);

  // Zero (the other cb)
  z1=zero;
  z2=zero;

  const int *tab = rb3[0].siteTable().slice();

  // By hand do the loop over the sites 
  for(int j=0; j < rb3[0].numSiteTable(); j++) {
    int i = tab[j];

    for(int spin=0; spin < 4; spin++) {
      // Loop over colors (case 3 only) 
      for(int col=0; col < 3; col++) { 

	z1.elem(i).elem(spin).elem(col).real() = 
	  a.elem().elem().elem().elem()*x.elem(i).elem(spin).elem(col).real()
	  + y.elem(i).elem(spin).elem(col).real();

	z1.elem(i).elem(spin).elem(col).imag() = 
	  a.elem().elem().elem().elem()*x.elem(i).elem(spin).elem(col).imag()
	  + y.elem(i).elem(spin).elem(col).imag();
      }
    }
  }

  z2[rb3[0]] = a*x + y;

  if( rb3[0].hasOrderedRep() ) { 
    int n_4vec = (rb3[0].end() - rb3[0].start() + 1);
    
    REAL64* xptr = (REAL64 *)&(x.elem(rb3[0].start()).elem(0).elem(0).real());
    REAL64* yptr = &(y.elem(rb3[0].start()).elem(0).elem(0).real());
    REAL64* zptr = &(z2.elem(rb3[0].start()).elem(0).elem(0).real());
    REAL64 ar = a.elem().elem().elem().elem();
    REAL64* aptr = &ar;

    vaxpyz4(zptr, aptr, xptr, yptr, n_4vec);
  }
  else { 
    const int *tab = rb3[0].siteTable().slice();
    for(int j=0; j < rb3[0].numSiteTable(); j++) {
      int i = tab[j];
      REAL64* xptr = (REAL64 *)&(x.elem(i).elem(0).elem(0).real());
      REAL64* yptr = &(y.elem(i).elem(0).elem(0).real());
      REAL64* zptr = &(z2.elem(i).elem(0).elem(0).real());
      REAL64 ar = a.elem().elem().elem().elem();
      REAL64* aptr = &ar;
      vaxpyz4(zptr, aptr, xptr, yptr, 1);
    }
  }

  for(int i=all.start(); i <= all.end(); i++) { 
    // Loop over spins
    for(int spin=0; spin < 4; spin++) {
      // Loop over colors (case 3 only) 
      for(int col=0; col < 3; col++) { 

	double realdiff = z1.elem(i).elem(spin).elem(col).real()
	  - z2.elem(i).elem(spin).elem(col).real();

	assertion( fabs(realdiff) < 1.0e-14 );

	double imagdiff = z1.elem(i).elem(spin).elem(col).imag()
	  - z2.elem(i).elem(spin).elem(col).imag();

	assertion( fabs(imagdiff) < 1.0e-14 );

      }
    }
  }
  
}


// Test 15. Check QDP++ against 'optimized'
void
testVaxpy4_RB30_3::run()
{
  Double a=Double(2.3);
  LatticeDiracFermionD3 x;
  LatticeDiracFermionD3 y;
  LatticeDiracFermionD3 z1;
  LatticeDiracFermionD3 z2;

  gaussian(x);
  gaussian(y);

  // QDP++ 
  z1[rb3[0]] = a*x + y;

  // Optimized
 if( rb3[0].hasOrderedRep() ) { 
    int n_4vec = (rb3[0].end() - rb3[0].start() + 1);
    
    REAL64* xptr = (REAL64 *)&(x.elem(rb3[0].start()).elem(0).elem(0).real());
    REAL64* yptr = &(y.elem(rb3[0].start()).elem(0).elem(0).real());
    REAL64* zptr = &(z2.elem(rb3[0].start()).elem(0).elem(0).real());
    REAL64 ar = a.elem().elem().elem().elem();
    REAL64* aptr = &ar;

    vaxpyz4(zptr, aptr, xptr, yptr, n_4vec);
  }
  else { 
    const int *tab = rb3[0].siteTable().slice();
    for(int j=0; j < rb3[0].numSiteTable(); j++) {
      int i = tab[j];
      REAL64* xptr = (REAL64 *)&(x.elem(i).elem(0).elem(0).real());
      REAL64* yptr = &(y.elem(i).elem(0).elem(0).real());
      REAL64* zptr = &(z2.elem(i).elem(0).elem(0).real());
      REAL64 ar = a.elem().elem().elem().elem();
      REAL64* aptr = &ar;
      vaxpyz4(zptr, aptr, xptr, yptr, 1);
    }
  }

  for(int i=all.start(); i <= all.end(); i++) { 
    // Loop over spins
    for(int spin=0; spin < 4; spin++) {
      // Loop over colors (case 3 only) 
      for(int col=0; col < 3; col++) { 

	double realdiff = z1.elem(i).elem(spin).elem(col).real()
	  - z2.elem(i).elem(spin).elem(col).real();

	assertion( fabs(realdiff) < 1.0e-14 );

	double imagdiff = z1.elem(i).elem(spin).elem(col).imag()
	  - z2.elem(i).elem(spin).elem(col).imag();

	assertion( fabs(imagdiff) < 1.0e-14 );

      }
    }
  }
}

// Test 16. Check QDP++ against handrolled
void
testVaxpy4_RB0_PEQ_1::run()
{
  Double a=Double(2.3);
  LatticeDiracFermionD3 x;
  LatticeDiracFermionD3 y;
  LatticeDiracFermionD3 z1;
  LatticeDiracFermionD3 z2;

  gaussian(x);
  gaussian(y);

  z1=zero;
  z2=zero;

  // QDP++
  z1[rb[0]]=y;
  z1[rb[0]] += a*x;

  // Optimized
  z2[rb[0]]=y;


  const int *tab = rb[0].siteTable().slice();
  for(int j=0; j < rb[0].numSiteTable(); j++) {
    int i = tab[j];
    for(int spin=0; spin < 4; spin++) { 
      for(int col=0; col < 3; col++) { 
	z2.elem(i).elem(spin).elem(col).real() +=
	  a.elem().elem().elem().elem()*x.elem(i).elem(spin).elem(col).real();
	z2.elem(i).elem(spin).elem(col).imag() +=
	  a.elem().elem().elem().elem()*x.elem(i).elem(spin).elem(col).imag();
      }
    }
  }
  

  for(int i=all.start(); i <= all.end(); i++) { 
    // Loop over spins
    for(int spin=0; spin < 4; spin++) {
      // Loop over colors (case 3 only) 
      for(int col=0; col < 3; col++) { 

	double realdiff = z1.elem(i).elem(spin).elem(col).real()
	  - z2.elem(i).elem(spin).elem(col).real();

	assertion( fabs(realdiff) < 1.0e-14 );

	double imagdiff = z1.elem(i).elem(spin).elem(col).imag()
	  - z2.elem(i).elem(spin).elem(col).imag();

	assertion( fabs(imagdiff) < 1.0e-14 );

      }
    }
  }
}

// Test 17. Check optimized against handrolled
void
testVaxpy4_RB0_PEQ_2::run()
{
  Double a=Double(2.3);
  LatticeDiracFermionD3 x;
  LatticeDiracFermionD3 y;
  LatticeDiracFermionD3 z1;
  LatticeDiracFermionD3 z2;

  gaussian(x);
  gaussian(y);

  z1=zero;
  z2=zero;

  // Handrolled
  z1[rb[0]]=y;
  const int *tab = rb[0].siteTable().slice();
  for(int j=0; j < rb[0].numSiteTable(); j++) {
    int i = tab[j];
    for(int spin=0; spin < 4; spin++) { 
      for(int col=0; col < 3; col++) { 
	z1.elem(i).elem(spin).elem(col).real() +=
	  a.elem().elem().elem().elem()*x.elem(i).elem(spin).elem(col).real();
	z1.elem(i).elem(spin).elem(col).imag() +=
	  a.elem().elem().elem().elem()*x.elem(i).elem(spin).elem(col).imag();
      }
    }
  }
  
  // Optimized
  z2[rb[0]]=y;

  if( rb[0].hasOrderedRep() ) { 
    int n_4vec = (rb[0].end() - rb[0].start() + 1);

    REAL64* xptr = (REAL64 *)&(x.elem(rb[0].start()).elem(0).elem(0).real());
    REAL64* zptr = &(z2.elem(rb[0].start()).elem(0).elem(0).real());
    REAL64 ar = a.elem().elem().elem().elem();
    REAL64* aptr = &ar;

    vaxpy4(zptr, aptr, xptr, n_4vec);
  }
  else { 
    const int *tab = rb[0].siteTable().slice();
    for(int j=0; j < rb[0].numSiteTable(); j++) {
      int i = tab[j];
      REAL64* xptr = (REAL64 *)&(x.elem(i).elem(0).elem(0).real());
      REAL64* zptr = &(z2.elem(i).elem(0).elem(0).real());
      REAL64 ar = a.elem().elem().elem().elem();
      REAL64* aptr = &ar;
      vaxpy4(zptr, aptr, xptr, 1);
    }
  }

  for(int i=all.start(); i <= all.end(); i++) { 
    // Loop over spins
    for(int spin=0; spin < 4; spin++) {
      // Loop over colors (case 3 only) 
      for(int col=0; col < 3; col++) { 

	double realdiff = z1.elem(i).elem(spin).elem(col).real()
	  - z2.elem(i).elem(spin).elem(col).real();

	assertion( fabs(realdiff) < 1.0e-14 );

	double imagdiff = z1.elem(i).elem(spin).elem(col).imag()
	  - z2.elem(i).elem(spin).elem(col).imag();

	assertion( fabs(imagdiff) < 1.0e-14 );

      }
    }
  }
}

// Test 18. Check QDP++ against optimized 
void
testVaxpy4_RB0_PEQ_3::run()
{
  Double a=Double(2.3);
  LatticeDiracFermionD3 x;
  LatticeDiracFermionD3 y;
  LatticeDiracFermionD3 z1;
  LatticeDiracFermionD3 z2;

  gaussian(x);
  gaussian(y);

  z1=zero;
  z2=zero;

  // Handrolled
  z1[rb[0]]=y;
  z1[rb[0]]+=a*x;
  
  // Optimized
  z2[rb[0]]=y;

  if( rb[0].hasOrderedRep() ) { 
    int n_4vec = (rb[0].end() - rb[0].start() + 1);

    REAL64* xptr = (REAL64 *)&(x.elem(rb[0].start()).elem(0).elem(0).real());
    REAL64* zptr = &(z2.elem(rb[0].start()).elem(0).elem(0).real());
    REAL64 ar = a.elem().elem().elem().elem();
    REAL64* aptr = &ar;

    vaxpy4(zptr, aptr, xptr, n_4vec);
  }
  else { 
    const int *tab = rb[0].siteTable().slice();
    for(int j=0; j < rb[0].numSiteTable(); j++) {
      int i = tab[j];
      REAL64* xptr = (REAL64 *)&(x.elem(i).elem(0).elem(0).real());
      REAL64* zptr = &(z2.elem(i).elem(0).elem(0).real());
      REAL64 ar = a.elem().elem().elem().elem();
      REAL64* aptr = &ar;
      vaxpy4(zptr, aptr, xptr, 1);
    }
  }

  for(int i=all.start(); i <= all.end(); i++) { 
    // Loop over spins
    for(int spin=0; spin < 4; spin++) {
      // Loop over colors (case 3 only) 
      for(int col=0; col < 3; col++) { 

	double realdiff = z1.elem(i).elem(spin).elem(col).real()
	  - z2.elem(i).elem(spin).elem(col).real();

	assertion( fabs(realdiff) < 1.0e-14 );

	double imagdiff = z1.elem(i).elem(spin).elem(col).imag()
	  - z2.elem(i).elem(spin).elem(col).imag();

	assertion( fabs(imagdiff) < 1.0e-14 );

      }
    }
  }
}



// Test 19. Check QDP++ against 'optimized'
void
testVaxpy4_RB30_PEQ::run()
{
  Double a=Double(2.3);
  LatticeDiracFermionD3 x;
  LatticeDiracFermionD3 y;
  LatticeDiracFermionD3 z1;
  LatticeDiracFermionD3 z2;

  gaussian(x);
  gaussian(y);

  z1=zero;
  z2=zero;

  // QDP++
  z1[rb3[0]]=y;
  z1[rb3[0]] += a*x;

  // Optimized
  z2[rb3[0]]=y;

 if( rb3[0].hasOrderedRep() ) { 
    int n_4vec = (rb3[0].end() - rb3[0].start() + 1);

    REAL64* xptr = (REAL64 *)&(x.elem(rb3[0].start()).elem(0).elem(0).real());
    REAL64* zptr = &(z2.elem(rb3[0].start()).elem(0).elem(0).real());
    REAL64 ar = a.elem().elem().elem().elem();
    REAL64* aptr = &ar;

    vaxpy4(zptr, aptr, xptr, n_4vec);
  }
  else { 
    const int *tab = rb3[0].siteTable().slice();
    for(int j=0; j < rb3[0].numSiteTable(); j++) {
      int i = tab[j];
      REAL64* xptr = (REAL64 *)&(x.elem(i).elem(0).elem(0).real());
      REAL64* zptr = &(z2.elem(i).elem(0).elem(0).real());
      REAL64 ar = a.elem().elem().elem().elem();
      REAL64* aptr = &ar;
      vaxpy4(zptr, aptr, xptr, 1);
    }
  }

  for(int i=all.start(); i <= all.end(); i++) { 
    // Loop over spins
    for(int spin=0; spin < 4; spin++) {
      // Loop over colors (case 3 only) 
      for(int col=0; col < 3; col++) { 

	double realdiff = z1.elem(i).elem(spin).elem(col).real()
	  - z2.elem(i).elem(spin).elem(col).real();

	assertion( fabs(realdiff) < 1.0e-14 );

	double imagdiff = z1.elem(i).elem(spin).elem(col).imag()
	  - z2.elem(i).elem(spin).elem(col).imag();

	assertion( fabs(imagdiff) < 1.0e-14 );

      }
    }
  }
}


