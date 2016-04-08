#ifndef TIME_DSLASH_NOQDP
#define TIME_DSLASH_NOQDP

enum Prec { FLOAT_PREC=0, HALF_PREC, DOUBLE_PREC };

class timeDslashNoQDP { 
public: 
 timeDslashNoQDP(int By_, int Bz_, int NCores_, int Sy_, int Sz_, int PadXY_, int PadXYZ_, int MinCt_, int niters_, bool c12, Prec precision_) : By(By_), Bz(Bz_), NCores(NCores_), Sy(Sy_), Sz(Sz_), PadXY(PadXY_), PadXYZ(PadXYZ_), MinCt(MinCt_), N_simt(Sy_*Sz_),  compress12(c12), iters(niters_), precision(precision_) {}


  void run(const int lattSize[], const int qmp_geom[]); 
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
  const int iters;

  template<typename FT, int V, int S, bool compress>
    void runTest(const int lattSize[], const int qmp_geom[]);

};

#endif
