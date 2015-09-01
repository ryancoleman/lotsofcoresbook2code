#include "qdp.h"
#include "testLocalVcdotDouble.h"
#include "unittest.h"

// #include "scalarsite_sse/sse_blas_vaxpy3_double.h"

using namespace QDP;
using namespace Assertions;

// Trinity of tests: Check SSE against Handrolled
//                   Check Handrolled against QDP++
//                   Check QDP++ against SSE


// Test 1. Check hand rolled against -- 'optimized'
void
testLocalVcdot4_1::run()
{

  REAL64 lsum_hand_re=(REAL64)(0);
  REAL64 lsum_hand_im=(REAL64)(0);
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

	lsum_hand_re +=(x.elem(i).elem(spin).elem(col).real()
		        * y.elem(i).elem(spin).elem(col).real());

	lsum_hand_re +=(x.elem(i).elem(spin).elem(col).imag()
		      * y.elem(i).elem(spin).elem(col).imag());


	lsum_hand_im -=(x.elem(i).elem(spin).elem(col).real()
			* y.elem(i).elem(spin).elem(col).imag());

	lsum_hand_im +=(x.elem(i).elem(spin).elem(col).imag()
			* y.elem(i).elem(spin).elem(col).real());

      }
    }
  }

  QDPInternal::globalSum(lsum_hand_re);
  QDPInternal::globalSum(lsum_hand_im);

  DComplex lsum_opt=zero;
  REAL64* sumptr=&(lsum_opt.elem().elem().elem().real());
  REAL64* xptr = (REAL64 *)&(x.elem(all.start()).elem(0).elem(0).real());
  REAL64* yptr = (REAL64 *)&(y.elem(all.start()).elem(0).elem(0).real());
  int n_4vec=all.end()-all.start()+1;
  local_vcdot4(sumptr, yptr, xptr, n_4vec);
  QDPInternal::globalSum(lsum_opt);



  REAL64 lsum_opt_re = lsum_opt.elem().elem().elem().real();
  REAL64 lsum_opt_im = lsum_opt.elem().elem().elem().imag();
  REAL64 diff = fabs(lsum_opt_re - lsum_hand_re);
  REAL64 dof  = (REAL64)(Layout::vol()*4*3*2);
  QDPIO::cout << std::endl << "\tDiff = " << diff << std::endl;
  QDPIO::cout << "\tDiff/dof = " << diff/dof << std::endl;
  assertion( diff/dof  < 1.0e-14 );

  diff = fabs(lsum_opt_im - lsum_hand_im);
  QDPIO::cout << std::endl << "\tDiff = " << diff << std::endl;
  QDPIO::cout << "\tDiff/dof = " << diff/dof << std::endl;
  assertion( diff/dof  < 1.0e-14 );

}


// Test 2. Check hand QDP++ against -- 'optimized'
void
testLocalVcdot4_2::run()
{
  LatticeDiracFermionD3 x;
  LatticeDiracFermionD3 y;
  gaussian(x);
  gaussian(y);
  DComplex lsum_qdp = innerProduct(y,x);
  DComplex lsum_opt=zero;

  REAL64* sumptr=&(lsum_opt.elem().elem().elem().real());
  REAL64* xptr = (REAL64 *)&(x.elem(all.start()).elem(0).elem(0).real());
  REAL64* yptr = (REAL64 *)&(y.elem(all.start()).elem(0).elem(0).real());
  int n_4vec=all.end()-all.start()+1;
  local_vcdot4(sumptr, yptr, xptr, n_4vec);
  QDPInternal::globalSum(lsum_opt);

  REAL64 lsum_qdp_re = lsum_qdp.elem().elem().elem().real();
  REAL64 lsum_qdp_im = lsum_qdp.elem().elem().elem().imag();

  REAL64 lsum_opt_re = lsum_opt.elem().elem().elem().real();
  REAL64 lsum_opt_im = lsum_opt.elem().elem().elem().imag();

  REAL64 diff = fabs(lsum_opt_re - lsum_qdp_re);
  REAL64 dof  = (REAL64)(Layout::vol()*4*3*2);
  QDPIO::cout << std::endl << "\tDiff = " << diff << std::endl;
  QDPIO::cout << "\tDiff/dof = " << diff/dof << std::endl;
  assertion( diff/dof  < 1.0e-14 );

  diff = fabs(lsum_opt_im - lsum_qdp_im);
  QDPIO::cout << std::endl << "\tDiff = " << diff << std::endl;
  QDPIO::cout << "\tDiff/dof = " << diff/dof << std::endl;
  assertion( diff/dof  < 1.0e-14  );

}


// Test 3. Check hand rolled against QDP++
void
testLocalVcdot4_3::run()
{

  LatticeDiracFermionD3 x;
  LatticeDiracFermionD3 y;
  gaussian(x);
  gaussian(y);

  REAL64 lsum_hand_re=(REAL64)(0);
  REAL64 lsum_hand_im=(REAL64)(0);

  // Loop over sites
  for(int i=all.start(); i <= all.end(); i++) { 
    // Loop over spins
    for(int spin=0; spin < 4; spin++) {
      // Loop over colors (case 3 only) 
      for(int col=0; col < 3; col++) { 

	lsum_hand_re +=x.elem(i).elem(spin).elem(col).real()
	  * y.elem(i).elem(spin).elem(col).real();

	lsum_hand_re +=x.elem(i).elem(spin).elem(col).imag()
	  * y.elem(i).elem(spin).elem(col).imag();

	lsum_hand_im +=x.elem(i).elem(spin).elem(col).imag()
	  * y.elem(i).elem(spin).elem(col).real();

	lsum_hand_im -=x.elem(i).elem(spin).elem(col).real()
	  * y.elem(i).elem(spin).elem(col).imag();

      }
    }
  }
  QDPInternal::globalSum(lsum_hand_re);
  QDPInternal::globalSum(lsum_hand_im);

  DComplex lsum_qdp = innerProduct(y,x);
  REAL64 lsum_qdp_re = lsum_qdp.elem().elem().elem().real();
  REAL64 lsum_qdp_im = lsum_qdp.elem().elem().elem().imag();

  REAL64 diff = fabs(lsum_hand_re - lsum_qdp_re);
  REAL64 dof  = (REAL64)(Layout::vol()*4*3*2);
  QDPIO::cout << std::endl << "\tDiff = " << toDouble(diff) << std::endl;
  QDPIO::cout << "\tDiff/dof = " << toDouble(diff/dof) << std::endl;
  assertion( toBool( diff/dof  < 1.0e-14 ) );

  diff = fabs(lsum_hand_im - lsum_qdp_im);
  QDPIO::cout << std::endl << "\tDiff = " << toDouble(diff) << std::endl;
  QDPIO::cout << "\tDiff/dof = " << toDouble(diff/dof) << std::endl;
  assertion( toBool( diff/dof  < 1.0e-14 ) );


}

