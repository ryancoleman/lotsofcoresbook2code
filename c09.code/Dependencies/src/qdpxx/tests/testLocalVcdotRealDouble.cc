#include "qdp.h"
#include "testLocalVcdotRealDouble.h"
#include "unittest.h"

// #include "scalarsite_sse/sse_blas_vaxpy3_double.h"

using namespace QDP;
using namespace Assertions;

// Trinity of tests: Check SSE against Handrolled
//                   Check Handrolled against QDP++
//                   Check QDP++ against SSE


// Test 1. Check hand rolled against -- 'optimized'
void
testLocalVcdotReal4_1::run()
{

  Double lsum_hand=Double(0);

  LatticeDiracFermionD3 x;
  LatticeDiracFermionD3 y;

  gaussian(x);
  gaussian(y);


  // Loop over sites
  for(int i=all.start(); i <= all.end(); i++) { 

    // Loop over spins
    for(int spin=0; spin < 4; spin++) {
      // Loop over colors (case 3 only) 
      for(int col=0; col < 3; col++) { 

	lsum_hand +=(x.elem(i).elem(spin).elem(col).real()
		     * y.elem(i).elem(spin).elem(col).real());

	lsum_hand +=(x.elem(i).elem(spin).elem(col).imag()
		     * y.elem(i).elem(spin).elem(col).imag());

      }
    }
  }


  Double lsum_opt=Double(0);
  REAL64* sumptr=&(lsum_opt.elem().elem().elem().elem());
  REAL64* xptr = (REAL64 *)&(x.elem(all.start()).elem(0).elem(0).real());
  REAL64* yptr = (REAL64 *)&(y.elem(all.start()).elem(0).elem(0).real());
  int n_4vec=all.end()-all.start()+1;
  local_vcdot_real4(sumptr, yptr, xptr, n_4vec);

  Double diff = fabs(lsum_opt - lsum_hand);
  Double dof  = Double(Layout::vol()*4*3*2);

  QDPIO::cout << std::endl << "\tDiff = " << toDouble(diff) << std::endl;
  QDPIO::cout << "\tDiff/dof = " << toDouble(diff/dof) << std::endl;
  
  assertion( toBool( diff/dof  < 1.0e-14 ) );

}


// Test 2. Check hand QDP++ against -- 'optimized'
void
testLocalVcdotReal4_2::run()
{
  LatticeDiracFermionD3 x;
  LatticeDiracFermionD3 y;
  gaussian(x);
  gaussian(y);
  Double lsum_qdp = innerProductReal(y,x);

  Double lsum_opt=Double(0);
  REAL64* sumptr=&(lsum_opt.elem().elem().elem().elem());
  REAL64* xptr = (REAL64 *)&(x.elem(all.start()).elem(0).elem(0).real());
  REAL64* yptr = (REAL64 *)&(y.elem(all.start()).elem(0).elem(0).real());
  int n_4vec=all.end()-all.start()+1;
  local_vcdot_real4(sumptr, yptr, xptr, n_4vec);
  QDPInternal::globalSum(lsum_opt);

  Double diff = fabs(lsum_opt - lsum_qdp);
  Double dof=Double(Layout::vol()*4*3*2);
  
  QDPIO::cout << std::endl <<  "\tDiff = " << diff << std::endl;
  QDPIO::cout <<  "\tDiff/d.o.f = " << diff/dof << std::endl;


  assertion( toBool(diff/dof < 1.0e-14) );

}


// Test 3. Check hand rolled against QDP++
void
testLocalVcdotReal4_3::run()
{

  LatticeDiracFermionD3 x;
  LatticeDiracFermionD3 y;
  gaussian(x);
  gaussian(y);

  Double lsum_hand=Double(0);

  // Loop over sites
  for(int i=all.start(); i <= all.end(); i++) { 
    // Loop over spins
    for(int spin=0; spin < 4; spin++) {
      // Loop over colors (case 3 only) 
      for(int col=0; col < 3; col++) { 

	lsum_hand +=x.elem(i).elem(spin).elem(col).real()
	  * y.elem(i).elem(spin).elem(col).real();

	lsum_hand +=x.elem(i).elem(spin).elem(col).imag()
	  * y.elem(i).elem(spin).elem(col).imag();

      }
    }
  }
  QDPInternal::globalSum(lsum_hand);

  Double lsum_qdp = innerProductReal(y,x);
  Double diff = fabs(lsum_qdp - lsum_hand);
  Double dof=Double(Layout::vol()*4*3*2);
  QDPIO::cout << std::endl << "\tDiff = " << diff << std::endl;
  QDPIO::cout << "\tDiff/d.o.f=" << diff/dof << std::endl;


  assertion( toBool(diff/dof < 1.0e-14) );
}

