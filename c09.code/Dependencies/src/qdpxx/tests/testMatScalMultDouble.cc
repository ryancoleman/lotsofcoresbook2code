#include "qdp.h"
#include "testMatScalMultDouble.h"
#include "unittest.h"

// #include "scalarsite_sse/sse_blas_vaxpy3_double.h"

using namespace QDP;
using namespace Assertions;

// Trinity of tests: Check SSE against Handrolled
//                   Check Handrolled against QDP++
//                   Check QDP++ against SSE


// Test 1. Check hand rolled against optimized
void
testScalMult1_1::run()
{
  Double a=Double(2.3);
  LatticeColorMatrixD3 x; 
  LatticeColorMatrixD3 z1, z2;

  gaussian(x);
  
  for(int site=all.start(); site <= all.end(); site++) { 
    for(int col1=0; col1 < 3; col1++) { 
      for(int col2=0; col2 < 3; col2++) { 
	z1.elem(site).elem().elem(col1,col2).real() = 
	  a.elem().elem().elem().elem()*x.elem(site).elem().elem(col1,col2).real();

	z1.elem(site).elem().elem(col1,col2).imag() = 
	  a.elem().elem().elem().elem()*x.elem(site).elem().elem(col1,col2).imag();
      }
    }
  }


  // Now the optimized one
  REAL64 *aptr = (REAL64 *)&(a.elem().elem().elem().elem());
  REAL64 *m1ptr = (REAL64 *)&(x.elem(0).elem().elem(0,0).real());
  REAL64 *m2ptr = (REAL64 *)&(z2.elem(0).elem().elem(0,0).real());
  int n_mat=all.end() - all.start() + 1;
  ssed_m_eq_scal_m(m2ptr, aptr, m1ptr, n_mat);

  for(int i=all.start(); i <= all.end(); i++) { 

    // Loop over spins
    for(int col1=0; col1< 3; col1++) {
      // Loop over colors (case 3 only) 
      for(int col2=0; col2 < 3; col2++) { 

	double realdiff = z1.elem(i).elem().elem(col1,col2).real()
	  - z2.elem(i).elem().elem(col1,col2).real();

	assertion( fabs(realdiff) < 1.0e-14 );

	double imagdiff = z1.elem(i).elem().elem(col1,col2).imag()
	  - z2.elem(i).elem().elem(col1,col2).imag();

	assertion( fabs(imagdiff) < 1.0e-14 );

      }
    }
  }
}

// Test 2. Check hand rolled against qdp++
void
testScalMult1_2::run()
{
  Double a=Double(2.3);
  LatticeColorMatrixD3 x; 
  LatticeColorMatrixD3 z1, z2;

  gaussian(x);
  
  for(int site=all.start(); site <= all.end(); site++) { 
    for(int col1=0; col1 < 3; col1++) { 
      for(int col2=0; col2 < 3; col2++) { 
	z1.elem(site).elem().elem(col1,col2).real() = 
	  a.elem().elem().elem().elem()*x.elem(site).elem().elem(col1,col2).real();

	z1.elem(site).elem().elem(col1,col2).imag() = 
	  a.elem().elem().elem().elem()*x.elem(site).elem().elem(col1,col2).imag();
      }
    }
  }


  z2 = a*x;

  for(int i=all.start(); i <= all.end(); i++) { 

    // Loop over spins
    for(int col1=0; col1< 3; col1++) {
      // Loop over colors (case 3 only) 
      for(int col2=0; col2 < 3; col2++) { 

	double realdiff = z1.elem(i).elem().elem(col1,col2).real()
	  - z2.elem(i).elem().elem(col1,col2).real();

	assertion( fabs(realdiff) < 1.0e-14 );

	double imagdiff = z1.elem(i).elem().elem(col1,col2).imag()
	  - z2.elem(i).elem().elem(col1,col2).imag();

	assertion( fabs(imagdiff) < 1.0e-14 );

      }
    }
  }

}

// Test 2. Check QDP++ against OPT
void
testScalMult1_3::run()
{
  Double a=Double(2.3);
  LatticeColorMatrixD3 x; 
  LatticeColorMatrixD3 z1, z2;

  gaussian(x);
  z1 = a*x;
  
  // Now the optimized one
  REAL64 *aptr = (REAL64 *)&(a.elem().elem().elem().elem());
  REAL64 *m1ptr = (REAL64 *)&(x.elem(0).elem().elem(0,0).real());
  REAL64 *m2ptr = (REAL64 *)&(z2.elem(0).elem().elem(0,0).real());
  int n_mat=all.end() - all.start() + 1;
  ssed_m_eq_scal_m(m2ptr, aptr, m1ptr, n_mat);

  for(int i=all.start(); i <= all.end(); i++) { 

    // Loop over spins
    for(int col1=0; col1< 3; col1++) {
      // Loop over colors (case 3 only) 
      for(int col2=0; col2 < 3; col2++) { 

	double realdiff = z1.elem(i).elem().elem(col1,col2).real()
	  - z2.elem(i).elem().elem(col1,col2).real();

	assertion( fabs(realdiff) < 1.0e-14 );

	double imagdiff = z1.elem(i).elem().elem(col1,col2).imag()
	  - z2.elem(i).elem().elem(col1,col2).imag();

	assertion( fabs(imagdiff) < 1.0e-14 );

      }
    }
  }

}

// Test 1. Check hand rolled against optimized
void
testScalMult2_1::run()
{
  Double a=Double(2.3);
  LatticeColorMatrixD3 x; 
  LatticeColorMatrixD3 z1, z2;

  gaussian(x);
  z1 = x;

  for(int site=all.start(); site <= all.end(); site++) { 
    for(int col1=0; col1 < 3; col1++) { 
      for(int col2=0; col2 < 3; col2++) { 
	z1.elem(site).elem().elem(col1,col2).real() *= 
	  a.elem().elem().elem().elem();

	z1.elem(site).elem().elem(col1,col2).imag() *= 
	  a.elem().elem().elem().elem();

      }
    }
  }

  z2 = x;
  // Now the optimized one
  REAL64 *aptr = (REAL64 *)&(a.elem().elem().elem().elem());
  REAL64 *m2ptr = (REAL64 *)&(z2.elem(0).elem().elem(0,0).real());
  int n_mat=all.end() - all.start() + 1;
  ssed_m_muleq_scal(m2ptr, aptr, n_mat);

  for(int i=all.start(); i <= all.end(); i++) { 

    // Loop over spins
    for(int col1=0; col1< 3; col1++) {
      // Loop over colors (case 3 only) 
      for(int col2=0; col2 < 3; col2++) { 

	double realdiff = z1.elem(i).elem().elem(col1,col2).real()
	  - z2.elem(i).elem().elem(col1,col2).real();

	assertion( fabs(realdiff) < 1.0e-14 );

	double imagdiff = z1.elem(i).elem().elem(col1,col2).imag()
	  - z2.elem(i).elem().elem(col1,col2).imag();

	assertion( fabs(imagdiff) < 1.0e-14 );

      }
    }
  }
}

// Test 2. Check hand rolled against qdp++
void
testScalMult2_2::run()
{
  Double a=Double(2.3);
  LatticeColorMatrixD3 x; 
  LatticeColorMatrixD3 z1, z2;

  gaussian(x);
  z1 = x;  
  for(int site=all.start(); site <= all.end(); site++) { 
    for(int col1=0; col1 < 3; col1++) { 
      for(int col2=0; col2 < 3; col2++) { 
	z1.elem(site).elem().elem(col1,col2).real() *= 
	  a.elem().elem().elem().elem();

	z1.elem(site).elem().elem(col1,col2).imag() *= 
	  a.elem().elem().elem().elem();
      }
    }
  }


  z2 = x;
  z2 *= a;

  for(int i=all.start(); i <= all.end(); i++) { 

    // Loop over spins
    for(int col1=0; col1< 3; col1++) {
      // Loop over colors (case 3 only) 
      for(int col2=0; col2 < 3; col2++) { 

	double realdiff = z1.elem(i).elem().elem(col1,col2).real()
	  - z2.elem(i).elem().elem(col1,col2).real();

	assertion( fabs(realdiff) < 1.0e-14 );

	double imagdiff = z1.elem(i).elem().elem(col1,col2).imag()
	  - z2.elem(i).elem().elem(col1,col2).imag();

	assertion( fabs(imagdiff) < 1.0e-14 );

      }
    }
  }

}

// Test 2. Check QDP++ against OPT
void
testScalMult2_3::run()
{
  Double a=Double(2.3);
  LatticeColorMatrixD3 x; 
  LatticeColorMatrixD3 z1, z2;

  gaussian(x);
  z1 =x;
  z1 *= a;

  z2 = x;

  // Now the optimized one
  REAL64 *aptr = (REAL64 *)&(a.elem().elem().elem().elem());
  REAL64 *m2ptr = (REAL64 *)&(z2.elem(0).elem().elem(0,0).real());
  int n_mat=all.end() - all.start() + 1;
  ssed_m_muleq_scal(m2ptr, aptr, n_mat);

  for(int i=all.start(); i <= all.end(); i++) { 

    // Loop over spins
    for(int col1=0; col1< 3; col1++) {
      // Loop over colors (case 3 only) 
      for(int col2=0; col2 < 3; col2++) { 

	double realdiff = z1.elem(i).elem().elem(col1,col2).real()
	  - z2.elem(i).elem().elem(col1,col2).real();

	assertion( fabs(realdiff) < 1.0e-14 );

	double imagdiff = z1.elem(i).elem().elem(col1,col2).imag()
	  - z2.elem(i).elem().elem(col1,col2).imag();

	assertion( fabs(imagdiff) < 1.0e-14 );

      }
    }
  }

}





