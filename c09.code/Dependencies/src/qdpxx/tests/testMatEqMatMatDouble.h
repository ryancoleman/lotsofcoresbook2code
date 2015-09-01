#ifndef TEST_MAT_EQ_MAT_MAT_DOUBLE
#define TEST_MAT_EQ_MAT_MAT_DOUBLE


#ifndef UNITTEST_H
#include "unittest.h"
#endif

#include "scalarsite_sse/sse_linalg_mm_su3_double.h"

class testMeqMM_1  : public TestFixture { public: void run(void); };
class testMeqMM_2  : public TestFixture { public: void run(void); };
class testMeqMM_3  : public TestFixture { public: void run(void); };

class testMPeqaMM_1  : public TestFixture { public: void run(void); };
class testMPeqaMM_2  : public TestFixture { public: void run(void); };
class testMPeqaMM_3  : public TestFixture { public: void run(void); };

#endif
