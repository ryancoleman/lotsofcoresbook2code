
#ifndef _DATA_TYPES_H_
#define _DATA_TYPES_H_

#include "instructions.h"

#if PRECISION == 1
#define FPTYPE float
#define LPTYPE short
#elif PRECISION == 2
#define FPTYPE double
#define LPTYPE float
#else
#error "Invalid FP PRECISION"
#endif

#ifndef USE_LP_SPINOR
typedef FPTYPE SpinorBaseType;
#define SpinorType	0
#else
typedef LPTYPE SpinorBaseType;
#define SpinorType	1
#endif

#ifndef USE_LP_GAUGE
typedef FPTYPE GaugeBaseType;
#define GaugeType	0
#else
typedef LPTYPE GaugeBaseType;
#define GaugeType	1
#endif

#ifndef USE_LP_CLOVER
typedef FPTYPE CloverBaseType;
#define CloverType	0
#else
typedef LPTYPE CloverBaseType;
#define CloverType	1
#endif

template<typename FT, int VL, int SL, bool compress12>
class data_types
{
public:
    typedef FT Spinor[4][3][2][SL];
#ifdef USE_PACKED_GAUGES
    typedef FT Gauge[8][(compress12? 2 : 3)][3][2][VL];
#else
    typedef FT Gauge[8][(compress12? 2 : 3)][3][2][SL];
#endif
#ifdef USE_PACKED_CLOVER
    typedef struct {
        FT diag1[6][VL];
        FT off_diag1[15][2][VL];
        FT diag2[6][VL];
        FT off_diag2[15][2][VL];
    } Clover;
#else
    typedef struct {
        FT diag1[6][SL];
        FT off_diag1[15][2][SL];
        FT diag2[6][SL];
        FT off_diag2[15][2][SL];
    } Clover;
#endif
};

void LoadSpinorElement(InstVector& ivector, const FVec& ret, string& base, string& offsets, int spin, int col, int reim, bool isFace, string mask, int dir = -1);
void LoadFullSpinor(InstVector& ivector, const FVec ret[4][3][2], string& base, string& offsets, string mask);
void StoreFullSpinor(InstVector& ivector, const FVec ret[4][3][2], string& base, string& offsets, int isStreaming = 0);
void StreamFullSpinor(InstVector& ivector, const FVec ret[4][3][2], string& base, string& offsets);
void PackHalfSpinor(InstVector& ivector, const FVec ret[2][3][2], string& base, int dir, string mask);
void UnpackHalfSpinor(InstVector& ivector, const FVec ret[2][3][2], string& base, int dir, string mask="");
void LoadFullGaugeDir(InstVector& ivector, const FVec ret[3][3][2], string& base, string& offsets, int dir, bool compress12);
void LoadFullCloverBlock(InstVector& ivector, const FVec diag[6], const FVec off_diag[15][2], string& base, string& offsets, int block);
void PrefetchL1FullSpinorDirIn(InstVector& ivector, const string& base, const string& off, int dir = -1, int type = 0);
void PrefetchL1FullSpinorOut(InstVector& ivector, const string& base, const string& off);
void PrefetchL2FullSpinorDirIn(InstVector& ivector, const string& base, const string& off, const string& pref_dist, int dir=-1);
void PrefetchL2FullSpinorOut(InstVector& ivector, const string& base, const string& off, const string& pref_dist);
void PrefetchL1FullGaugeDirIn(InstVector& ivector, const string& base, const string& off, int dir, bool compress12);
void PrefetchL2FullGaugeDirIn(InstVector& ivector, const string& base, const string& off, int dir, const string& pref_dist, bool compress12);
void PrefetchL2FullGaugeIn(InstVector& ivector, const string& base, const string& off, const string& pref_dist, bool compress12);
void PrefetchL1FullCloverBlockIn(InstVector& ivector, const string& base, const string& off, int block);
void PrefetchL2FullCloverIn(InstVector& ivector, const string& base, const string& off, const string& pref_dist);
void PrefetchL1HalfSpinorDir(InstVector& ivector, const string& base, int dir, bool isPrefforWrite, int type);
void PrefetchL2HalfSpinorDir(InstVector& ivector, const string& base, const string& pref_dist, int dir, bool isPrefforWrite, int type);
#endif // _DATA_TYPES_H_