#include "qdp.h"
#include "unittest.h"
#include "testvol.h"

#include "testVaxpyDouble.h"
#include "testVaypxDouble.h"
#include "testVaxmyDouble.h"
#include "testVaxpbyDouble.h"
#include "testVScalDouble.h"
#include "testLocalSumSqDouble.h"
#include "testLocalVcdotRealDouble.h"
#include "testLocalVcdotDouble.h"

using namespace QDP;

int main(int argc, char **argv)
{
  // Initialize UnitTest jig
  TestRunner  tests(&argc, &argv, nrow_in);
  QDPIO::cout << "Volume= { " << Layout::lattSize()[0]
	      << " , " << Layout::lattSize()[1]
	      << " , " << Layout::lattSize()[2]
	      << " , " << Layout::lattSize()[3] << " } " << std::endl;

  tests.addTest(new testVaxpy4_1(), "testVaxpy4_1" );
  tests.addTest(new testVaxpy4_2(), "testVaxpy4_2" );
  tests.addTest(new testVaxpy4_3(), "testVaxpy4_3" );

  tests.addTest(new testVaypx4_1(), "testVaypx4_1" );
  tests.addTest(new testVaypx4_2(), "testVaypx4_2" );
  tests.addTest(new testVaypx4_3(), "testVaypx4_3" );

  tests.addTest(new testVaxpy4_RB0_1(), "testVaxpy4_RB0_1" );
  tests.addTest(new testVaxpy4_RB0_2(), "testVaxpy4_RB0_2" );
  tests.addTest(new testVaxpy4_RB0_3(), "testVaxpy4_RB0_3" );

  tests.addTest(new testVaxpy4_RB1_1(), "testVaxpy4_RB1_1" );
  tests.addTest(new testVaxpy4_RB1_2(), "testVaxpy4_RB1_2" );
  tests.addTest(new testVaxpy4_RB1_3(), "testVaxpy4_RB1_3" );

  tests.addTest(new testVaxpy4_RB30_1(), "testVaxpy4_RB30_1" );
  tests.addTest(new testVaxpy4_RB30_2(), "testVaxpy4_RB30_2" );
  tests.addTest(new testVaxpy4_RB30_3(), "testVaxpy4_RB30_3" );

  tests.addTest(new testVaxpy4_RB31_1(), "testVaxpy4_RB31_1" );
  tests.addTest(new testVaxpy4_RB31_2(), "testVaxpy4_RB31_2" );
  tests.addTest(new testVaxpy4_RB31_3(), "testVaxpy4_RB31_3" );

  tests.addTest(new testVaxpy4_RB0_PEQ_1(), "testVaxpy4_RB0_PEQ_1" );
  tests.addTest(new testVaxpy4_RB0_PEQ_2(), "testVaxpy4_RB0_PEQ_2" );
  tests.addTest(new testVaxpy4_RB0_PEQ_3(), "testVaxpy4_RB0_PEQ_3" );

  tests.addTest(new testVaxpy4_RB30_PEQ(), "testVaxpy4_RB30_PEQ" );


  tests.addTest(new testVaxmyz4_1(), "testVaxmyz4_1" );
  tests.addTest(new testVaxmyz4_2(), "testVaxmyz4_2" );
  tests.addTest(new testVaxmyz4_3(), "testVaxmyz4_3" );

  tests.addTest(new testVaxmy4_1(), "testVaxmy4_1" );
  tests.addTest(new testVaxmy4_2(), "testVaxmy4_2" );
  tests.addTest(new testVaxmy4_3(), "testVaxmy4_3" );

  tests.addTest(new testVaxpbyz4_1(), "testVaxpbyz4_1" );
  //  tests.addTest(new testVaxpbyz4_2(), "testVaxpbyz4_2" );
  tests.addTest(new testVaxpbyz4_3(), "testVaxpbyz4_3" );

  tests.addTest(new testVaxpby4_1(), "testVaxpby4_1" );
  tests.addTest(new testVaxpby4_2(), "testVaxpby4_2" );
  tests.addTest(new testVaxpby4_3(), "testVaxpby4_3" );


  tests.addTest(new testVaxmbyz4_1(), "testVaxmbyz4_1" );
  tests.addTest(new testVaxmbyz4_2(), "testVaxmbyz4_2" );
  tests.addTest(new testVaxmbyz4_3(), "testVaxmbyz4_3" );

  tests.addTest(new testVaxmby4_1(), "testVaxmby4_1" );
  tests.addTest(new testVaxmby4_2(), "testVaxmby4_2" );
  tests.addTest(new testVaxmby4_3(), "testVaxmby4_3" );

  tests.addTest(new testVScal4_1(), "testVScal4_1" );
  tests.addTest(new testVScal4_2(), "testVScal4_2" );
  tests.addTest(new testVScal4_3(), "testVScal4_3" );

  tests.addTest(new testLocalSumSq4_1(), "testLocalSumSq4_1" );
  tests.addTest(new testLocalSumSq4_2(), "testLocalSumSq4_2" );
  tests.addTest(new testLocalSumSq4_3(), "testLocalSumSq4_3" );

  tests.addTest(new testLocalVcdotReal4_1(), "testLocalVcdotReal4_1" );
  tests.addTest(new testLocalVcdotReal4_2(), "testLocalVcdotReal4_2" );
  tests.addTest(new testLocalVcdotReal4_3(), "testLocalVcdotReal4_3" );

  tests.addTest(new testLocalVcdot4_1(), "testLocalVcdot4_1" );
  tests.addTest(new testLocalVcdot4_2(), "testLocalVcdot4_2" );
  tests.addTest(new testLocalVcdot4_3(), "testLocalVcdot4_3" );



  // Run all tests
  tests.run();

  // Testjig is destroyed
  tests.summary();
}

