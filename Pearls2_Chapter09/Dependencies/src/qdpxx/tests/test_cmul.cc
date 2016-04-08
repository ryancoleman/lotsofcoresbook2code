#include "qdp.h"
#include "unittest.h"
#include "testvol.h"


#include "testCMulDouble.h"


using namespace QDP;

int main(int argc, char **argv)
{
  // Initialize UnitTest jig
  TestRunner  tests(&argc, &argv, nrow_in);
  QDPIO::cout << "Volume= { " << Layout::lattSize()[0]
	      << " , " << Layout::lattSize()[1]
	      << " , " << Layout::lattSize()[2]
	      << " , " << Layout::lattSize()[3] << " } " << std::endl;

  
  tests.addTest(new testCMul(), "testCMul" );
  tests.addTest(new testCMadd(), "testCMadd" );
  tests.addTest(new testConjMul(), "testConjMul" );
  tests.addTest(new testConjMadd(), "testConjMadd" );
  tests.addTest(new testCCMul(), "testCCMul" );
  tests.addTest(new testCCMadd(), "testCCMadd" );

  // Run all tests
  tests.run();

  // Testjig is destroyed
  tests.summary();
}

