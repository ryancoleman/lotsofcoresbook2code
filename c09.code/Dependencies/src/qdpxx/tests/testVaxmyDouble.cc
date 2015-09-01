#include "qdp.h"
#include "testVaxmyDouble.h"
#include "unittest.h"
// #include <qmt.h>

// #include "scalarsite_sse/sse_blas_vaxpy3_double.h"

using namespace QDP;
using namespace Assertions;

// Trinity of tests: Check SSE against Handrolled
//                   Check Handrolled against QDP++
//                   Check QDP++ against SSE


// Test 1. Check hand rolled against -- 'optimized'
void
testVaxmyz4_1::run()
{
  Double a=Double(2.3);
  LatticeDiracFermionD3 x;
  LatticeDiracFermionD3 y;
  LatticeDiracFermionD3 z1;
  LatticeDiracFermionD3 z2;

  gaussian(x);
  gaussian(y);

  // Loop over sites
  for(int i=all.start(); i <= all.end(); i++) { 
    // Loop over spins
    for(int spin=0; spin < 4; spin++) {
      // Loop over colors (case 3 only) 
      for(int col=0; col < 3; col++) { 

	z1.elem(i).elem(spin).elem(col).real() = 
	  a.elem().elem().elem().elem()*x.elem(i).elem(spin).elem(col).real()
	  - y.elem(i).elem(spin).elem(col).real();

	z1.elem(i).elem(spin).elem(col).imag() = 
	  a.elem().elem().elem().elem()*x.elem(i).elem(spin).elem(col).imag()
	  - y.elem(i).elem(spin).elem(col).imag();
      }
    }
  }

  int n_4vec = (all.end() - all.start() + 1);
  REAL64* xptr = (REAL64 *)&(x.elem(all.start()).elem(0).elem(0).real());
  REAL64* yptr = &(y.elem(all.start()).elem(0).elem(0).real());
  REAL64* zptr = &(z2.elem(all.start()).elem(0).elem(0).real());
  REAL64 ar = a.elem().elem().elem().elem();
  REAL64* aptr = &ar;

  vaxmyz4(zptr, aptr, xptr, yptr, n_4vec);

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


// Test 2. Check hand QDP++ against -- 'optimized'
void
testVaxmyz4_2::run()
{
  Double a=Double(2.3);
  LatticeDiracFermionD3 x;
  LatticeDiracFermionD3 y;
  LatticeDiracFermionD3 z1;
  LatticeDiracFermionD3 z2;

  gaussian(x);
  gaussian(y);

  z1 = a*x - y;
 

  int n_4vec = (all.end() - all.start() + 1);
  REAL64* xptr = (REAL64 *)&(x.elem(all.start()).elem(0).elem(0).real());
  REAL64* yptr = (REAL64 *)&(y.elem(all.start()).elem(0).elem(0).real());
  REAL64* zptr = &(z2.elem(all.start()).elem(0).elem(0).real());
  REAL64 ar = a.elem().elem().elem().elem();
  REAL64* aptr = &ar;

  vaxmyz4(zptr, aptr, xptr, yptr, n_4vec);


  for(int i=all.start(); i <= all.end(); i++) { 
    // Loop over spins
    for(int spin=0; spin < 4; spin++) {
      // Loop over colors (case 3 only) 
      for(int col=0; col < 3; col++) { 

	double realdiff = z1.elem(i).elem(spin).elem(col).real()
	  - z2.elem(i).elem(spin).elem(col).real();

	//QDPIO::cout << "rediff=" << realdiff << endl;
	assertion( fabs(realdiff) < 1.0e-13 );


	double imagdiff = z1.elem(i).elem(spin).elem(col).imag()
	  - z2.elem(i).elem(spin).elem(col).imag();
	//QDPIO::cout << "imagdiff=" << imagdiff << endl;
	assertion( fabs(imagdiff) < 1.0e-14 );

      }
    }
  }

}


// Test 3. Check hand rolled against QDP++
void
testVaxmyz4_3::run()
{
  Double a=Double(2.3);
  LatticeDiracFermionD3 x;
  LatticeDiracFermionD3 y;
  LatticeDiracFermionD3 z1;
  LatticeDiracFermionD3 z2;

  gaussian(x);
  gaussian(y);

  // Loop over sites
  for(int i=all.start(); i <= all.end(); i++) { 
    // Loop over spins
    for(int spin=0; spin < 4; spin++) {
      // Loop over colors (case 3 only) 
      for(int col=0; col < 3; col++) { 

	z1.elem(i).elem(spin).elem(col).real() = 
	  a.elem().elem().elem().elem()*x.elem(i).elem(spin).elem(col).real()
	  - y.elem(i).elem(spin).elem(col).real();

	z1.elem(i).elem(spin).elem(col).imag() = 
	  a.elem().elem().elem().elem()*x.elem(i).elem(spin).elem(col).imag()
	  - y.elem(i).elem(spin).elem(col).imag();
      }
    }
  }

  z2 = a*x - y;

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


// Test 1. Check hand rolled against -- 'optimized'
void
testVaxmy4_1::run()
{
  Double a=Double(2.3);
  LatticeDiracFermionD3 x;
  LatticeDiracFermionD3 y;
  LatticeDiracFermionD3 z1;
  LatticeDiracFermionD3 z2;

  gaussian(x);
  gaussian(y);

  z1 = y;

  // Loop over sites
  for(int i=all.start(); i <= all.end(); i++) { 
    // Loop over spins
    for(int spin=0; spin < 4; spin++) {
      // Loop over colors (case 3 only) 
      for(int col=0; col < 3; col++) { 

	z1.elem(i).elem(spin).elem(col).real() = 
	  a.elem().elem().elem().elem()*x.elem(i).elem(spin).elem(col).real()
	  - z1.elem(i).elem(spin).elem(col).real();

	z1.elem(i).elem(spin).elem(col).imag() = 
	  a.elem().elem().elem().elem()*x.elem(i).elem(spin).elem(col).imag()
	  - z1.elem(i).elem(spin).elem(col).imag();
      }
    }
  }

  z2 = y;

  int n_4vec = (all.end() - all.start() + 1);
  REAL64* xptr = (REAL64 *)&(x.elem(all.start()).elem(0).elem(0).real());
  REAL64* zptr = &(z2.elem(all.start()).elem(0).elem(0).real());
  REAL64 ar = a.elem().elem().elem().elem();
  REAL64* aptr = &ar;

  vaxmy4(zptr, aptr, xptr, n_4vec);

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


// Test 2. Check hand QDP++ against -- 'optimized'
void
testVaxmy4_2::run()
{
  Double a=Double(2.3);
  LatticeDiracFermionD3 x;
  LatticeDiracFermionD3 y;
  LatticeDiracFermionD3 z1;
  LatticeDiracFermionD3 z2;

  gaussian(x);
  gaussian(y);

  z1=y;
  z1 = a*x - z1;

  z2 = y;

  int n_4vec = (all.end() - all.start() + 1);
  REAL64* xptr = (REAL64 *)&(x.elem(all.start()).elem(0).elem(0).real());
  REAL64* yptr = &(y.elem(all.start()).elem(0).elem(0).real());
  REAL64* zptr = &(z2.elem(all.start()).elem(0).elem(0).real());
  REAL64 ar = a.elem().elem().elem().elem();
  REAL64* aptr = &ar;

  vaxmy4(zptr, aptr, xptr, n_4vec);

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


// Test 3. Check hand rolled against QDP++
void
testVaxmy4_3::run()
{
  Double a=Double(2.3);
  LatticeDiracFermionD3 x;
  LatticeDiracFermionD3 y;
  LatticeDiracFermionD3 z1;
  LatticeDiracFermionD3 z2;

  gaussian(x);
  gaussian(y);

  z1 = y;

  // Loop over sites
  for(int i=all.start(); i <= all.end(); i++) { 
    // Loop over spins
    for(int spin=0; spin < 4; spin++) {
      // Loop over colors (case 3 only) 
      for(int col=0; col < 3; col++) { 

	z1.elem(i).elem(spin).elem(col).real() = 
	  a.elem().elem().elem().elem()*x.elem(i).elem(spin).elem(col).real()
	  - z1.elem(i).elem(spin).elem(col).real();

	z1.elem(i).elem(spin).elem(col).imag() = 
	  a.elem().elem().elem().elem()*x.elem(i).elem(spin).elem(col).imag()
	  - z1.elem(i).elem(spin).elem(col).imag();
      }
    }
  }

  z2 = y;

  z2 = a*x - z2;

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




