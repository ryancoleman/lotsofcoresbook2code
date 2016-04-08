#include "qdp.h"
#include "unittest.h"
#include "testvol.h"


#include "testMatScalMultDouble.h"
#include "testMatPeqMatDouble.h"
#include "testMatEqMatMatDouble.h"
#include "testMatEqMatHermDouble.h"
#include "testMatEqHermMatDouble.h"
#include "testMatEqHermHermDouble.h"

using namespace QDP;

int main(int argc, char **argv)
{
  // Initialize UnitTest jig
  TestRunner  tests(&argc, &argv, nrow_in);
  QDPIO::cout << "Volume= { " << Layout::lattSize()[0]
	      << " , " << Layout::lattSize()[1]
	      << " , " << Layout::lattSize()[2]
	      << " , " << Layout::lattSize()[3] << " } " << std::endl;

  
  tests.addTest(new testScalMult1_1(), "testScalMult1_1" );
  tests.addTest(new testScalMult1_2(), "testScalMult1_2" );
  tests.addTest(new testScalMult1_3(), "testScalMult1_3" );

  tests.addTest(new testScalMult2_1(), "testScalMult2_1" );
  tests.addTest(new testScalMult2_2(), "testScalMult2_2" );
  tests.addTest(new testScalMult2_3(), "testScalMult2_3" );

  tests.addTest(new testMPeqM_1(), "testMPeqM_1" );
  tests.addTest(new testMPeqM_2(), "testMPeqM_2" );
  tests.addTest(new testMPeqM_3(), "testMPeqM_3" );

  tests.addTest(new testMMeqM_1(), "testMMeqM_1" );
  tests.addTest(new testMMeqM_2(), "testMMeqM_2" );
  tests.addTest(new testMMeqM_3(), "testMMeqM_3" );

  tests.addTest(new testMPeqH_1(), "testMPeqH_1" );
  tests.addTest(new testMPeqH_2(), "testMPeqH_2" );
  tests.addTest(new testMPeqH_3(), "testMPeqH_3" );

  tests.addTest(new testMMeqH_1(), "testMMeqH_1" );
  tests.addTest(new testMMeqH_2(), "testMMeqH_2" );
  tests.addTest(new testMMeqH_3(), "testMMeqH_3" );
  
  tests.addTest(new testMeqMM_1(), "testMeqMM_1" );
  tests.addTest(new testMeqMM_2(), "testMeqMM_2" );
  tests.addTest(new testMeqMM_3(), "testMeqMM_3" );

  tests.addTest(new testMeqMH_1(), "testMeqMH_1" );
  tests.addTest(new testMeqMH_2(), "testMeqMH_2" );
  tests.addTest(new testMeqMH_3(), "testMeqMH_3" );

  tests.addTest(new testMeqHM_1(), "testMeqHM_1" );
  tests.addTest(new testMeqHM_2(), "testMeqHM_2" );
  tests.addTest(new testMeqHM_3(), "testMeqHM_3" );

  tests.addTest(new testMeqHH_1(), "testMeqHH_1" );
  tests.addTest(new testMeqHH_2(), "testMeqHH_2" );
  tests.addTest(new testMeqHH_3(), "testMeqHH_3" );

  tests.addTest(new testMPeqaMM_1(), "testMPeqaMM_1" );
  tests.addTest(new testMPeqaMM_2(), "testMPeqaMM_2" );
  tests.addTest(new testMPeqaMM_3(), "testMPeqaMM_3" );

  tests.addTest(new testMPeqaMH_1(), "testMPeqaMH_1" );
  tests.addTest(new testMPeqaMH_2(), "testMPeqaMH_2" );
  tests.addTest(new testMPeqaMH_3(), "testMPeqaMH_3" );

  tests.addTest(new testMPeqaHM_1(), "testMPeqaHM_1" );
  tests.addTest(new testMPeqaHM_2(), "testMPeqaHM_2" );
  tests.addTest(new testMPeqaHM_3(), "testMPeqaHM_3" );

  tests.addTest(new testMPeqaHH_1(), "testMPeqaHH_1" );
  tests.addTest(new testMPeqaHH_2(), "testMPeqaHH_2" );
  tests.addTest(new testMPeqaHH_3(), "testMPeqaHH_3" );

  // Run all tests
  tests.run();

  // Testjig is destroyed
  tests.summary();
}

