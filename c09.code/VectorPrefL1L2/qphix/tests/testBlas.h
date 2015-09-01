#ifndef TEST_BLAS
#define TEST_BLAS

#include "qphix/print_utils.h"
class testBlas { 
public: 
 testBlas(int NCores_, int Sy_, int Sz_, int PadXY_, int PadXYZ_, int niters_) : NCores(NCores_), Sy(Sy_), Sz(Sz_), PadXY(PadXY_), PadXYZ(PadXYZ_), iters(niters_) {}
  void run(const int lattSize[], const int qmp_geom[]); 
 private:
  const int NCores;
  const int Sy;
  const int Sz;
  const int PadXY;
  const int PadXYZ;
  const int iters;
};

#endif
