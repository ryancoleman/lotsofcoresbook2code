#include "qdp.h"
#include "unittest.h"
#include "testvol.h"


#include "timeMatEqMatMatDouble.h"
#include "timeMatEqMatHermDouble.h"
#include "timeMatEqHermMatDouble.h"
#include "timeMatEqHermHermDouble.h"

using namespace QDP;

int main(int argc, char **argv)
{
  // Initialize UnitTest jig
  TestRunner  tests(&argc, &argv, nrow_in);
  QDPIO::cout << "Volume= { " << Layout::lattSize()[0]
	      << " , " << Layout::lattSize()[1]
	      << " , " << Layout::lattSize()[2]
	      << " , " << Layout::lattSize()[3] << " } " << std::endl;

  
  tests.addTest(new timeMeqMM_QDP(), "timeMeqMM_QDP" );
  tests.addTest(new timeMeqMM(), "timeMeqMM" );

  tests.addTest(new timeMPeqaMM_QDP(), "timeMPeqaMM_QDP" );
  tests.addTest(new timeMPeqaMM(), "timeMPeqaMM" );

  tests.addTest(new timeMeqMH_QDP(), "timeMeqMH_QDP" );
  tests.addTest(new timeMeqMH(), "timeMeqMH" );

  tests.addTest(new timeMPeqaMH_QDP(), "timeMPeqaMH_QDP" );
  tests.addTest(new timeMPeqaMH(), "timeMPeqaMH" );

  tests.addTest(new timeMeqHM_QDP(), "timeMeqHM_QDP" );
  tests.addTest(new timeMeqHM(), "timeMeqHM" );

  tests.addTest(new timeMPeqaHM_QDP(), "timeMPeqaHM_QDP" );
  tests.addTest(new timeMPeqaHM(), "timeMPeqaHM" );

  tests.addTest(new timeMeqHH_QDP(), "timeMeqHH_QDP" );
  tests.addTest(new timeMeqHH(), "timeMeqHH" );

  tests.addTest(new timeMPeqaHH_QDP(), "timeMPeqaHH_QDP" );
  tests.addTest(new timeMPeqaHH(), "timeMPeqaHH" );

  // Run all tests
  tests.run();

  // Testjig is destroyed
  tests.summary();
}

