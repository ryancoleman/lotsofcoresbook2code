#ifndef TEST_DSLASH_FULL
#define TEST_DSLASH_FULL


#ifndef UNITTEST_H
#include "unittest.h"
#endif

enum Prec { FLOAT_PREC=0, HALF_PREC, DOUBLE_PREC };

class testDslashFull : public TestFixture { 
public: 
 testDslashFull(int By_, int Bz_, int NCores_, int Sy_, int Sz_, int PadXY_, int PadXYZ_, int MinCt_, bool c12, Prec precision_, const int soalen_) : By(By_), Bz(Bz_), NCores(NCores_), Sy(Sy_), Sz(Sz_), PadXY(PadXY_), PadXYZ(PadXYZ_), MinCt(MinCt_), N_simt(Sy_*Sz_), compress12(c12), precision(precision_), soalen(soalen_) {}
  void run(void); 
 private:
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
  const int soalen;
  Seed rng_seed;

  template<typename T, int V, int S, typename U, typename Phi>
  void testDslashWrapper(const U& u) { 
    for(int t_bc=+1; t_bc >= -1; t_bc-=2) {
      if( compress12  ) { 
	testDslash<T,V,S,true,U,Phi>(u,t_bc);
      }
      else {
	testDslash<T,V,S,false,U,Phi>(u,t_bc);
      }
    }
  }
  
  template<typename T, int V, int S, typename U, typename Phi>
  void testDslashAChiMBDPsiWrapper(const U& u) {
    for(int t_bc=1; t_bc >= -1; t_bc-=2) {
      if( compress12  ) { 
	testDslashAChiMBDPsi<T,V,S,true,U,Phi>(u,t_bc);
      }
      else {
	testDslashAChiMBDPsi<T,V,S,false,U,Phi>(u,t_bc);
      }
    }
  }
  
  template<typename T, int V, int S, typename U, typename Phi>
  void testMWrapper(const U& u) {
    for(int t_bc=-1; t_bc <= +1; t_bc+=2) {
      if( compress12  ) { 
	testM<T,V,S,true,U,Phi>(u,t_bc);
      }
      else {
	testM<T,V,S,false,U,Phi>(u,t_bc);
      }
    }
  }
  
  template<typename T, int V, int S, typename U, typename Phi>
  void testCGWrapper(const U& u) 
  {
    for(int t_bc=-1; t_bc <= +1; t_bc+=2) {
      if( compress12  ) { 
	testCG<T,V,S,true,U,Phi>(u,t_bc);
      }
      else {
	testCG<T,V,S,false,U,Phi>(u,t_bc);
      }
    }
  }

  template<typename T, int V, int S, typename U, typename Phi>
  void testBiCGStabWrapper(const U& u)
  {
    //for(int t_bc=-1; t_bc <= +1; t_bc+=2) {
    int t_bc=-1;
      if( compress12  ) { 
	testBiCGStab<T,V,S,true,U,Phi>(u,t_bc);
      }
      else {
	testBiCGStab<T,V,S,false,U,Phi>(u,t_bc);
      }
      //}
  }
  
  template<typename T, int V, int S, bool compress, typename U, typename Phi>
    void testDslash(const U& u, int t_bc);

  template<typename T, int V, int S, bool compress, typename U, typename Phi>
    void testDslashAChiMBDPsi(const U& u, int t_bc);

  template<typename T, int V, int S, bool compress, typename U, typename Phi>
    void testM(const U& u, int t_bc);

  template<typename T, int V, int S, bool compress, typename U, typename Phi>
    void testCG(const U& u, int t_bc);

  template<typename T, int V, int S, bool compress, typename U, typename Phi>
    void testBiCGStab(const U& u, int t_bc);


  template<typename T1, int VEC1, int SOA1, bool compress, typename T2, int VEC2, int SOA2, typename U, typename Phi>
    void testRichardson(const U& u, int t_bc);
};

#endif
