#ifndef TEST_MAT_PEQ_MAT_DOUBLE
#define TEST_MAT_PEQ_MAT_DOUBLE


#ifndef UNITTEST_H
#include "unittest.h"
#endif

#include "scalarsite_sse/sse_linalg_mm_su3_double.h"

class testMPeqM_1  : public TestFixture { public: void run(void); };
class testMPeqM_2  : public TestFixture { public: void run(void); };
class testMPeqM_3  : public TestFixture { public: void run(void); };

class testMMeqM_1  : public TestFixture { public: void run(void); };
class testMMeqM_2  : public TestFixture { public: void run(void); };
class testMMeqM_3  : public TestFixture { public: void run(void); };

class testMPeqH_1  : public TestFixture { public: void run(void); };
class testMPeqH_2  : public TestFixture { public: void run(void); };
class testMPeqH_3  : public TestFixture { public: void run(void); };

class testMMeqH_1  : public TestFixture { public: void run(void); };
class testMMeqH_2  : public TestFixture { public: void run(void); };
class testMMeqH_3  : public TestFixture { public: void run(void); };


#endif
