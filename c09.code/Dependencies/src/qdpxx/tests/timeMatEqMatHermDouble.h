#ifndef TIME_MAT_EQ_MAT_HERM_DOUBLE
#define TIME_MAT_EQ_MAT_HERM_DOUBLE


#ifndef UNITTEST_H
#include "unittest.h"
#endif

#include "scalarsite_sse/sse_linalg_mm_su3_double.h"

class timeMeqMH_QDP  : public TestFixture { public: void run(void); };
class timeMeqMH      : public TestFixture { public: void run(void); };

class timeMPeqaMH_QDP  : public TestFixture { public: void run(void); };
class timeMPeqaMH      : public TestFixture { public: void run(void); };


#endif
