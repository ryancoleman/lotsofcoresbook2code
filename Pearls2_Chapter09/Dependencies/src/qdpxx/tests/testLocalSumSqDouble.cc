#include "qdp.h"
#include "testLocalSumSqDouble.h"
#include "unittest.h"

// #include "scalarsite_sse/sse_blas_vaxpy3_double.h"

using namespace QDP;
using namespace Assertions;

// Trinity of tests: Check SSE against Handrolled
//                   Check Handrolled against QDP++
//                   Check QDP++ against SSE


// Test 1. Check hand rolled against -- 'optimized'
void
testLocalSumSq4_1::run()
{

  Double lsum_hand=Double(0);
  LatticeDiracFermionD3 x;
  gaussian(x);


  // Loop over sites
  // This is NOT THREADED but it is test code. That is OK.
  for(int i=all.start(); i <= all.end(); i++) { 
    // Loop over spins
    for(int spin=0; spin < 4; spin++) {
      // Loop over colors (case 3 only) 
      for(int col=0; col < 3; col++) { 

	lsum_hand +=(x.elem(i).elem(spin).elem(col).real()
		     * x.elem(i).elem(spin).elem(col).real());

	lsum_hand +=(x.elem(i).elem(spin).elem(col).imag()
		     * x.elem(i).elem(spin).elem(col).imag());

      }
    }
  }
  QDPInternal::globalSum(lsum_hand);

  Double lsum_opt=Double(0);
  REAL64* sumptr=&(lsum_opt.elem().elem().elem().elem());
  REAL64* xptr = (REAL64 *)&(x.elem(all.start()).elem(0).elem(0).real());
  int n_4vec=all.end()-all.start()+1;

  // This may eventually be threaded
  local_sumsq4(sumptr, xptr, n_4vec);
  QDPInternal::globalSum(lsum_opt);

  Double diff = fabs(lsum_opt - lsum_hand);
  Double dof  = Double(Layout::vol()*4*3*2);

  QDPIO::cout << std::endl << "\tDiff = " << toDouble(diff) << std::endl;
  QDPIO::cout << "\tDiff/dof = " << toDouble(diff/dof) << std::endl;
  
  assertion( toBool( diff/dof  < 1.0e-13 ) );

}


// Test 2. Check hand QDP++ against -- 'optimized'
void
testLocalSumSq4_2::run()
{
  
  LatticeDiracFermionD3 x;

  // These may eventually be threaded
  gaussian(x);
  Double lsum_qdp = norm2(x);

  Double lsum_opt=Double(0);
  REAL64* sumptr=&(lsum_opt.elem().elem().elem().elem());
  REAL64* xptr = (REAL64 *)&(x.elem(all.start()).elem(0).elem(0).real());
  int n_4vec=all.end()-all.start()+1;
  // This may be threaded under the hood.

  local_sumsq4(sumptr, xptr, n_4vec);
  QDPInternal::globalSum(lsum_opt);

  Double diff = fabs(lsum_opt - lsum_qdp);
  Double dof=Double(Layout::vol()*4*3*2);
  
  QDPIO::cout << std::endl <<  "\tDiff = " << diff << std::endl;
  QDPIO::cout <<  "\tDiff/d.o.f = " << diff/dof << std::endl;


  assertion( toBool(diff/dof < 1.0e-14) );

}


// Test 3. Check hand rolled against QDP++
void
testLocalSumSq4_3::run()
{

  LatticeDiracFermionD3 x;
  gaussian(x);
  Double lsum_hand=Double(0);

  // Loop over sites
  // Not threaded but that is OK
  for(int i=all.start(); i <= all.end(); i++) { 
    // Loop over spins
    for(int spin=0; spin < 4; spin++) {
      // Loop over colors (case 3 only) 
      for(int col=0; col < 3; col++) { 

	lsum_hand +=x.elem(i).elem(spin).elem(col).real()
	  * x.elem(i).elem(spin).elem(col).real();

	lsum_hand +=x.elem(i).elem(spin).elem(col).imag()
	  * x.elem(i).elem(spin).elem(col).imag();

      }
    }
  }
  QDPInternal::globalSum(lsum_hand);

  // This should be threaded...
  Double lsum_qdp = norm2(x);


  Double diff = fabs(lsum_qdp - lsum_hand);
  Double dof=Double(Layout::vol()*4*3*2);
  QDPIO::cout << std::endl << "\tDiff = " << diff << std::endl;
  QDPIO::cout << "\tDiff/d.o.f=" << diff/dof << std::endl;


  assertion( toBool(diff/dof < 1.0e-14) );
}

