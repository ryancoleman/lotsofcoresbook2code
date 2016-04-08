#ifndef TEST_CLOVDSLASH_FULL
#define TEST_CLOVDSLASH_FULL


#ifndef UNITTEST_H
#include "unittest.h"
#endif

enum Prec { FLOAT_PREC=0, HALF_PREC, DOUBLE_PREC };


class testClovDslashFull : public TestFixture { 
public: 
 testClovDslashFull(int By_, int Bz_, int NCores_, int Sy_, int Sz_, int PadXY_, int PadXYZ_, int MinCt_, bool c12, Prec precision_) : By(By_), Bz(Bz_), NCores(NCores_), Sy(Sy_), Sz(Sz_), PadXY(PadXY_), PadXYZ(PadXYZ_), MinCt(MinCt_), N_simt(Sy_*Sz_), compress12(c12), precision(precision_) {}

  // Toplevel test function wrapper
  void run(void); 
 private:
  // Templated test function
  template<typename FT, int V, int S, bool compress, typename U, typename Phi>
    void runTest(void);

  const int By;
  const int Bz;
  const int NCores;
  const int Sy;
  const int Sz;
  const int PadXY;
  const int PadXYZ;
  const int MinCt;
  const int N_simt;
  const bool compress12;
  const Prec precision;
};

#endif
