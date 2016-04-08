#include "qdp.h"
#include "testMatPeqMatDouble.h"
#include "unittest.h"

// #include "scalarsite_sse/sse_blas_vaxpy3_double.h"

using namespace QDP;
using namespace Assertions;

// Trinity of tests: Check SSE against Handrolled
//                   Check Handrolled against QDP++
//                   Check QDP++ against SSE


// Test 1. Check hand rolled against optimized
void
testMPeqM_1::run()
{

  LatticeColorMatrixD3 x; 
  LatticeColorMatrixD3 z1, z2;

  gaussian(x);
  gaussian(z1);
  z2 = z1;

  for(int site=all.start(); site <= all.end(); site++) { 
    for(int col1=0; col1 < 3; col1++) { 
      for(int col2=0; col2 < 3; col2++) { 
	z1.elem(site).elem().elem(col1,col2).real() += 
	  x.elem(site).elem().elem(col1,col2).real();

	z1.elem(site).elem().elem(col1,col2).imag() += 
	  x.elem(site).elem().elem(col1,col2).imag();

      }
    }
  }


  // Now the optimized one
  REAL64 *m1ptr = (REAL64 *)&(x.elem(0).elem().elem(0,0).real());
  REAL64 *m2ptr = (REAL64 *)&(z2.elem(0).elem().elem(0,0).real());

  int n_mat=all.end() - all.start() + 1;
  ssed_m_peq_m(m2ptr, m1ptr, n_mat);

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
testMPeqM_2::run()
{
  LatticeColorMatrixD3 x; 
  LatticeColorMatrixD3 z1, z2;

  gaussian(x);
  gaussian(z1);
  z2 = z1;
  for(int site=all.start(); site <= all.end(); site++) { 
    for(int col1=0; col1 < 3; col1++) { 
      for(int col2=0; col2 < 3; col2++) { 
	z1.elem(site).elem().elem(col1,col2).real() += 
	  x.elem(site).elem().elem(col1,col2).real();

	z1.elem(site).elem().elem(col1,col2).imag() += 
	  x.elem(site).elem().elem(col1,col2).imag();

      }
    }
  }
  z2 += x;

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
testMPeqM_3::run()
{
  LatticeColorMatrixD3 x; 
  LatticeColorMatrixD3 z1, z2;

  gaussian(x);
  gaussian(z1);
  z2 = z1;

  z1 += x;

  // Now the optimized one
  REAL64 *m1ptr = (REAL64 *)&(x.elem(0).elem().elem(0,0).real());
  REAL64 *m2ptr = (REAL64 *)&(z2.elem(0).elem().elem(0,0).real());
  int n_mat=all.end() - all.start() + 1;
  ssed_m_peq_m(m2ptr,  m1ptr, n_mat);

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



void
testMMeqM_1::run()
{

  LatticeColorMatrixD3 x; 
  LatticeColorMatrixD3 z1, z2;

  gaussian(x);
  gaussian(z1);
  z2 = z1;

  for(int site=all.start(); site <= all.end(); site++) { 
    for(int col1=0; col1 < 3; col1++) { 
      for(int col2=0; col2 < 3; col2++) { 
	z1.elem(site).elem().elem(col1,col2).real() -= 
	  x.elem(site).elem().elem(col1,col2).real();

	z1.elem(site).elem().elem(col1,col2).imag() -= 
	  x.elem(site).elem().elem(col1,col2).imag();

      }
    }
  }


  // Now the optimized one
  REAL64 *m1ptr = (REAL64 *)&(x.elem(0).elem().elem(0,0).real());
  REAL64 *m2ptr = (REAL64 *)&(z2.elem(0).elem().elem(0,0).real());

  int n_mat=all.end() - all.start() + 1;
  ssed_m_meq_m(m2ptr, m1ptr, n_mat);

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
testMMeqM_2::run()
{
  LatticeColorMatrixD3 x; 
  LatticeColorMatrixD3 z1, z2;

  gaussian(x);
  gaussian(z1);
  z2 = z1;
  for(int site=all.start(); site <= all.end(); site++) { 
    for(int col1=0; col1 < 3; col1++) { 
      for(int col2=0; col2 < 3; col2++) { 
	z1.elem(site).elem().elem(col1,col2).real() -= 
	  x.elem(site).elem().elem(col1,col2).real();

	z1.elem(site).elem().elem(col1,col2).imag() -= 
	  x.elem(site).elem().elem(col1,col2).imag();

      }
    }
  }
  z2 -= x;

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
testMMeqM_3::run()
{
  LatticeColorMatrixD3 x; 
  LatticeColorMatrixD3 z1, z2;

  gaussian(x);
  gaussian(z1);
  z2 = z1;

  z1 -= x;

  // Now the optimized one
  REAL64 *m1ptr = (REAL64 *)&(x.elem(0).elem().elem(0,0).real());
  REAL64 *m2ptr = (REAL64 *)&(z2.elem(0).elem().elem(0,0).real());
  int n_mat=all.end() - all.start() + 1;
  ssed_m_meq_m(m2ptr,  m1ptr, n_mat);

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
testMPeqH_1::run()
{

  LatticeColorMatrixD3 x; 
  LatticeColorMatrixD3 z1, z2;

  gaussian(x);
  gaussian(z1);
  z2 = z1;

  for(int site=all.start(); site <= all.end(); site++) { 
    for(int col1=0; col1 < 3; col1++) { 
      for(int col2=0; col2 < 3; col2++) { 
	z1.elem(site).elem().elem(col1,col2).real() += 
	  x.elem(site).elem().elem(col2,col1).real();

	z1.elem(site).elem().elem(col1,col2).imag() -= 
	  x.elem(site).elem().elem(col2,col1).imag();

      }
    }
  }


  // Now the optimized one
  REAL64 *m1ptr = (REAL64 *)&(x.elem(0).elem().elem(0,0).real());
  REAL64 *m2ptr = (REAL64 *)&(z2.elem(0).elem().elem(0,0).real());

  int n_mat=all.end() - all.start() + 1;
  ssed_m_peq_h(m2ptr, m1ptr, n_mat);

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
testMPeqH_2::run()
{
  LatticeColorMatrixD3 x; 
  LatticeColorMatrixD3 z1, z2;

  gaussian(x);
  gaussian(z1);
  z2 = z1;
  for(int site=all.start(); site <= all.end(); site++) { 
    for(int col1=0; col1 < 3; col1++) { 
      for(int col2=0; col2 < 3; col2++) { 
	z1.elem(site).elem().elem(col1,col2).real() += 
	  x.elem(site).elem().elem(col2,col1).real();

	z1.elem(site).elem().elem(col1,col2).imag() -= 
	  x.elem(site).elem().elem(col2,col1).imag();

      }
    }
  }
  z2 += adj(x);

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
testMPeqH_3::run()
{
  LatticeColorMatrixD3 x; 
  LatticeColorMatrixD3 z1, z2;

  gaussian(x);
  gaussian(z1);
  z2 = z1;

  z1 += adj(x);

  // Now the optimized one
  REAL64 *m1ptr = (REAL64 *)&(x.elem(0).elem().elem(0,0).real());
  REAL64 *m2ptr = (REAL64 *)&(z2.elem(0).elem().elem(0,0).real());
  int n_mat=all.end() - all.start() + 1;
  ssed_m_peq_h(m2ptr,  m1ptr, n_mat);

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
testMMeqH_1::run()
{

  LatticeColorMatrixD3 x; 
  LatticeColorMatrixD3 z1, z2;

  gaussian(x);
  gaussian(z1);
  z2 = z1;

  for(int site=all.start(); site <= all.end(); site++) { 
    for(int col1=0; col1 < 3; col1++) { 
      for(int col2=0; col2 < 3; col2++) { 
	z1.elem(site).elem().elem(col1,col2).real() -= 
	  x.elem(site).elem().elem(col2,col1).real();

	z1.elem(site).elem().elem(col1,col2).imag() += 
	  x.elem(site).elem().elem(col2,col1).imag();

      }
    }
  }


  // Now the optimized one
  REAL64 *m1ptr = (REAL64 *)&(x.elem(0).elem().elem(0,0).real());
  REAL64 *m2ptr = (REAL64 *)&(z2.elem(0).elem().elem(0,0).real());

  int n_mat=all.end() - all.start() + 1;
  ssed_m_meq_h(m2ptr, m1ptr, n_mat);

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
testMMeqH_2::run()
{
  LatticeColorMatrixD3 x; 
  LatticeColorMatrixD3 z1, z2;

  gaussian(x);
  gaussian(z1);
  z2 = z1;
  for(int site=all.start(); site <= all.end(); site++) { 
    for(int col1=0; col1 < 3; col1++) { 
      for(int col2=0; col2 < 3; col2++) { 
	z1.elem(site).elem().elem(col1,col2).real() -= 
	  x.elem(site).elem().elem(col2,col1).real();

	z1.elem(site).elem().elem(col1,col2).imag() += 
	  x.elem(site).elem().elem(col2,col1).imag();

      }
    }
  }
  z2 -= adj(x);

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
testMMeqH_3::run()
{
  LatticeColorMatrixD3 x; 
  LatticeColorMatrixD3 z1, z2;

  gaussian(x);
  gaussian(z1);
  z2 = z1;

  z1 -= adj(x);

  // Now the optimized one
  REAL64 *m1ptr = (REAL64 *)&(x.elem(0).elem().elem(0,0).real());
  REAL64 *m2ptr = (REAL64 *)&(z2.elem(0).elem().elem(0,0).real());
  int n_mat=all.end() - all.start() + 1;
  ssed_m_meq_h(m2ptr,  m1ptr, n_mat);

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










