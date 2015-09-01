#include "qdp.h"
#include "testMatEqMatHermDouble.h"
#include "unittest.h"

// #include "scalarsite_sse/sse_blas_vaxpy3_double.h"

using namespace QDP;
using namespace Assertions;

// Trinity of tests: Check SSE against Handrolled
//                   Check Handrolled against QDP++
//                   Check QDP++ against SSE


// Test 1. Check hand rolled against optimized
void
testMeqMH_1::run()
{

  LatticeColorMatrixD3 x,y; 
  LatticeColorMatrixD3 z1, z2;

  gaussian(x);
  gaussian(y);

  for(int site=all.start(); site <= all.end(); site++) { 
    for(int i=0; i < 3; i++) { 
      for(int j=0; j < 3; j++) { 
	z1.elem(site).elem().elem(i,j).real() = (REAL64)0;
	z1.elem(site).elem().elem(i,j).imag() = (REAL64)0;

	for(int k=0; k < 3; k++) { 
	  z1.elem(site).elem().elem(i,j).real() += 
	    x.elem(site).elem().elem(i,k).real()
	    * y.elem(site).elem().elem(j,k).real();

	  z1.elem(site).elem().elem(i,j).real() +=
	     x.elem(site).elem().elem(i,k).imag()
	    * y.elem(site).elem().elem(j,k).imag();

	  z1.elem(site).elem().elem(i,j).imag() -=
	     x.elem(site).elem().elem(i,k).real()
	    * y.elem(site).elem().elem(j,k).imag();

	  z1.elem(site).elem().elem(i,j).imag() += 
	    x.elem(site).elem().elem(i,k).imag()
	    * y.elem(site).elem().elem(j,k).real();

	}

      }
    }
  }


  // Now the optimized one
  REAL64 *xptr = (REAL64 *)&(x.elem(0).elem().elem(0,0).real());
  REAL64 *yptr = (REAL64 *)&(y.elem(0).elem().elem(0,0).real());

  REAL64 *zptr = (REAL64 *)&(z2.elem(0).elem().elem(0,0).real());

  int n_mat=all.end() - all.start() + 1;
  ssed_m_eq_mh(zptr, xptr, yptr, n_mat);

  
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



// Test 1. Check hand rolled against QDP++
void
testMeqMH_2::run()
{

  LatticeColorMatrixD3 x,y; 
  LatticeColorMatrixD3 z1, z2;

  gaussian(x);
  gaussian(y);

  for(int site=all.start(); site <= all.end(); site++) { 
    for(int i=0; i < 3; i++) { 
      for(int j=0; j < 3; j++) { 
	z1.elem(site).elem().elem(i,j).real() = (REAL64)0;
	z1.elem(site).elem().elem(i,j).imag() = (REAL64)0;

	for(int k=0; k < 3; k++) { 
	  z1.elem(site).elem().elem(i,j).real() += 
	    x.elem(site).elem().elem(i,k).real()
	    * y.elem(site).elem().elem(j,k).real();

	  z1.elem(site).elem().elem(i,j).real() +=
	     x.elem(site).elem().elem(i,k).imag()
	    * y.elem(site).elem().elem(j,k).imag();

	  z1.elem(site).elem().elem(i,j).imag() -=
	     x.elem(site).elem().elem(i,k).real()
	    * y.elem(site).elem().elem(j,k).imag();

	  z1.elem(site).elem().elem(i,j).imag() += 
	    x.elem(site).elem().elem(i,k).imag()
	    * y.elem(site).elem().elem(j,k).real();

	}

      }
    }
  }

  z2 = x*adj(y);

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

// Check hand QDP++ against optimized
void
testMeqMH_3::run()
{

  LatticeColorMatrixD3 x,y; 
  LatticeColorMatrixD3 z1, z2;

  gaussian(x);
  gaussian(y);
  z1= x*adj(y);

  // Now the optimized one
  REAL64 *xptr = (REAL64 *)&(x.elem(0).elem().elem(0,0).real());
  REAL64 *yptr = (REAL64 *)&(y.elem(0).elem().elem(0,0).real());

  REAL64 *zptr = (REAL64 *)&(z2.elem(0).elem().elem(0,0).real());

  int n_mat=all.end() - all.start() + 1;
  ssed_m_eq_mh(zptr, xptr, yptr, n_mat);

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
testMPeqaMH_1::run()
{

  LatticeColorMatrixD3 x,y; 
  LatticeColorMatrixD3 z1, z2;
  Double a(0.6);

  gaussian(x);
  gaussian(y);
  gaussian(z1);
  z2=z1;

  for(int site=all.start(); site <= all.end(); site++) { 
    for(int i=0; i < 3; i++) { 
      for(int j=0; j < 3; j++) { 
	
	for(int k=0; k < 3; k++) { 
	  z1.elem(site).elem().elem(i,j).real() += 
	    a.elem().elem().elem().elem()
	    *x.elem(site).elem().elem(i,k).real()
	    * y.elem(site).elem().elem(j,k).real();

	  z1.elem(site).elem().elem(i,j).real() +=
	    a.elem().elem().elem().elem()
	    *x.elem(site).elem().elem(i,k).imag()
	    * y.elem(site).elem().elem(j,k).imag();

	  z1.elem(site).elem().elem(i,j).imag() += 
	    a.elem().elem().elem().elem()
	    *x.elem(site).elem().elem(i,k).imag()
	    * y.elem(site).elem().elem(j,k).real();

	  z1.elem(site).elem().elem(i,j).imag() -=
	    a.elem().elem().elem().elem()
	    *x.elem(site).elem().elem(i,k).real()
	    * y.elem(site).elem().elem(j,k).imag();
	}

      }
    }
  }


  // Now the optimized one
  REAL64 *xptr = (REAL64 *)&(x.elem(0).elem().elem(0,0).real());
  REAL64 *yptr = (REAL64 *)&(y.elem(0).elem().elem(0,0).real());
  REAL64 *aptr = (REAL64 *)&(a.elem().elem().elem().elem());
  REAL64 *zptr = (REAL64 *)&(z2.elem(0).elem().elem(0,0).real());

  int n_mat=all.end() - all.start() + 1;
  ssed_m_peq_amh(zptr, aptr, xptr, yptr, n_mat);

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

// Test 1. Check hand rolled against QDP++
void
testMPeqaMH_2::run()
{

  LatticeColorMatrixD3 x,y; 
  LatticeColorMatrixD3 z1, z2;
  Double a(-0.4);

  gaussian(x);
  gaussian(y);
  gaussian(z1);
  z2 = z1;

  for(int site=all.start(); site <= all.end(); site++) { 
    for(int i=0; i < 3; i++) { 
      for(int j=0; j < 3; j++) { 
	for(int k=0; k < 3; k++) { 
	  z1.elem(site).elem().elem(i,j).real() +=
	    a.elem().elem().elem().elem() 
	    *x.elem(site).elem().elem(i,k).real()
	    * y.elem(site).elem().elem(j,k).real();

	  z1.elem(site).elem().elem(i,j).real() +=
	    a.elem().elem().elem().elem()
	    *x.elem(site).elem().elem(i,k).imag()
	    * y.elem(site).elem().elem(j,k).imag();

	  z1.elem(site).elem().elem(i,j).imag() += 
	    a.elem().elem().elem().elem()
	    *x.elem(site).elem().elem(i,k).imag()
	    * y.elem(site).elem().elem(j,k).real();

	  z1.elem(site).elem().elem(i,j).imag() -=
	    a.elem().elem().elem().elem()
	    *x.elem(site).elem().elem(i,k).real()
	    * y.elem(site).elem().elem(j,k).imag();
	}

      }
    }
  }

  LatticeColorMatrixD3 t=a*x;
  z2 += t*adj(y);

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
testMPeqaMH_3::run()
{

  LatticeColorMatrixD3 x,y; 
  LatticeColorMatrixD3 z1, z2;
  Double a(0.6);

  gaussian(x);
  gaussian(y);
  gaussian(z1);
  z2=z1;

  LatticeColorMatrixD3 t=a*x;
  z1 += t*adj(y);

  // Now the optimized one
  REAL64 *xptr = (REAL64 *)&(x.elem(0).elem().elem(0,0).real());
  REAL64 *yptr = (REAL64 *)&(y.elem(0).elem().elem(0,0).real());
  REAL64 *aptr = (REAL64 *)&(a.elem().elem().elem().elem());
  REAL64 *zptr = (REAL64 *)&(z2.elem(0).elem().elem(0,0).real());

  int n_mat=all.end() - all.start() + 1;
  ssed_m_peq_amh(zptr, aptr, xptr, yptr, n_mat);

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




