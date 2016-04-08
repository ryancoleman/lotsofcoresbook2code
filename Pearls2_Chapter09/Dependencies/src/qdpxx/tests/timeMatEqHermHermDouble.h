#ifndef TIME_MAT_EQ_HERM_HERM_DOUBLE
#define TIME_MAT_EQ_HERM_HERM_DOUBLE


#ifndef UNITTEST_H
#include "unittest.h"
#endif

#include "scalarsite_sse/sse_linalg_mm_su3_double.h"

class timeMeqHH_QDP  : public TestFixture { public: void run(void); };
class timeMeqHH      : public TestFixture { public: void run(void); };

class timeMPeqaHH_QDP  : public TestFixture { public: void run(void); };
class timeMPeqaHH      : public TestFixture { public: void run(void); };


#endif
