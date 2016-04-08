#ifndef TEST_CMUL_H
#define TEST_CMUL_H

#ifndef UNITTEST_H
#include "unittest.h"
#endif

class testCMul   : public TestFixture { public: void run(void); };
class testCMadd  : public TestFixture { public: void run(void); };

class testConjMul   : public TestFixture { public: void run(void); };
class testConjMadd  : public TestFixture { public: void run(void); };

class testCCMul   : public TestFixture { public: void run(void); };
class testCCMadd  : public TestFixture { public: void run(void); };


#endif
