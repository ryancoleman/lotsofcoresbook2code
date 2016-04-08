#ifndef TIME_MAT_EQ_MAT_MAT_DOUBLE
#define TIME_MAT_EQ_MAT_MAT_DOUBLE


#ifndef UNITTEST_H
#include "unittest.h"
#endif

#include "scalarsite_sse/sse_linalg_mm_su3_double.h"

class timeMeqMM_QDP  : public TestFixture { public: void run(void); };
class timeMeqMM      : public TestFixture { public: void run(void); };

class timeMPeqaMM_QDP  : public TestFixture { public: void run(void); };
class timeMPeqaMM      : public TestFixture { public: void run(void); };


#endif
