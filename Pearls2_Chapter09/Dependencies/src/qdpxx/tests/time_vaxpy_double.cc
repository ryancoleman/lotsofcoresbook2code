#include "unittest.h"
#include "testvol.h"

#include <string>
#include "timeVaxpyDouble.h"

using namespace QDP;

int main(int argc, char **argv)
{
  // Initialize UnitTest jig
  TestRunner  tests(&argc, &argv, nrow_in);
  QDPIO::cout << "Volume= { " << Layout::lattSize()[0]
	      << " , " << Layout::lattSize()[1]
	      << " , " << Layout::lattSize()[2]
	      << " , " << Layout::lattSize()[3] << " } " << std::endl;


  // This behaves as expected
  tests.addTest(new time_VAXPBYZ_double(), "time_AXPBYZ_double" );
  tests.addTest(new time_VAXPBYZ_float(), "time_AXPBYZ_float" );

 
  // tests.addTest(new time_VAXPYZ(), "time_AXPYZ" );
  // tests.addTest(new time_VAXPY(), "time_AXPY" );
  // tests.addTest(new time_VAXMY(), "time_AXMY" );


  //tests.addTest(new time_LOCAL_SUMSQ(), "time_LOCAL_SUMSQ");
  // tests.addTest(new time_LOCAL_VCDOT(), "time_LOCAL_VCDOT");
  // tests.addTest(new time_LOCAL_VCDOT_REAL(), "time_LOCAL_VCDOT_REAL");
  // tests.addTest(new time_VAXMYZ(), "time_AXMYZ" );
  // tests.addTest(new time_VSCAL(), "time_VSCAL");
  // tests.addTest(new time_SUMSQ(), "time_SUMSQ");
  // tests.addTest(new time_VCDOT(), "time_VCDOT");
  // tests.addTest(new time_VCDOT_REAL(), "time_VCDOT_REAL");
  // tests.addTest(new time_QDP_PEQ(), "time_QDP_PEQ" );
  // tests.addTest(new time_QDP_AXPYZ(), "time_QDP_AXPYZ" );

 
  tests.run();
  // Testjig is destroyed
  tests.summary();
}

