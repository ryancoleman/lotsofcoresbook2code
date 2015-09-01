
#include "data_types.h"
#include "instructions.h"

typedef SpinorBaseType Spinor[4][3][2][SOALEN];
#ifdef USE_PACKED_GAUGES
typedef GaugeBaseType Gauge[8][3][3][2][VECLEN];
#else
typedef GaugeBaseType Gauge[8][3][3][2][SOALEN];
#endif
#ifdef USE_PACKED_CLOVER
typedef struct {
    CloverBaseType diag1[6][VECLEN];
    CloverBaseType off_diag1[15][2][VECLEN];
    CloverBaseType diag2[6][VECLEN];
    CloverBaseType off_diag2[15][2][VECLEN];
} Clover;
#else
typedef struct {
    CloverBaseType diag1[6][SOALEN];
    CloverBaseType off_diag1[15][2][SOALEN];
    CloverBaseType diag2[6][SOALEN];
    CloverBaseType off_diag2[15][2][SOALEN];
} Clover;
#endif

string serialize_data_types(bool compress12)
{
    std::ostringstream buf;
#ifdef USE_PACKED_GAUGES
    buf << "#define USE_PACKED_GAUGES" << endl;
#else
    buf << "// USE_PACKED_GAUGES not defined" << endl;
#endif

#ifdef USE_PACKED_CLOVER
    buf << "#define USE_PACKED_CLOVER" << endl;
#else
    buf << "// USE_PACKED_CLOVER not defined" << endl;
#endif

    buf << "#define PRECISION	" << PRECISION << endl;
    buf << "#define SOALEN	" << SOALEN << endl;
    buf << "#define VECLEN	" << VECLEN << endl;

    return buf.str();
}

void readFVecSpecialized(InstVector& ivector, const FVec& ret, GatherAddress *a, string mask)
{
#if (SOALEN == VECLEN) && defined(USE_SHUFFLES)
    loadFVec(ivector, ret, a->getAddr(0), mask);
#else
#ifndef USE_LDUNPK
    gatherFVec(ivector, ret, a, mask);
#else
    initFVec(ivector, ret);

    for(int i = 0; i < VECLEN; i += SOALEN) {
        loadSOAFVec(ivector, ret, a->getAddr(i), i/SOALEN, SOALEN, mask);
    }

#endif
#endif
}

void writeFVecSpecialized(InstVector& ivector, const FVec& ret, GatherAddress *a, int isStreaming)
{
#if (SOALEN == VECLEN) && defined(USE_SHUFFLES)
    storeFVec(ivector, ret, a->getAddr(0), isStreaming);
#else
#ifndef USE_PKST
    scatterFVec(ivector, ret, a);
#else

    for(int i = 0; i < VECLEN; i += SOALEN) {
        storeSOAFVec(ivector, ret, a->getAddr(i), i/SOALEN, SOALEN);
    }

#endif
#endif
}

void readUFVec(InstVector& ivector, const FVec& ret, GatherAddress *a, int forward, string &mask)
{
    ivector.push_back(new InitFVec(ret));

    for(int i = 0; i < VECLEN; i += SOALEN) {
        loadSplitSOAFVec(ivector, ret, a->getAddr(i), a->getAddr(i+(forward ? SOALEN-1 : 1)), i/SOALEN, SOALEN, forward, mask);
    }
}

void readFVecSpinor(InstVector& ivector, const FVec& ret, string& base, string& offset, int spin, int color, int reim, string& mask)
{
    readFVecSpecialized(ivector, ret, new GatherAddress(new SpinorAddress(base,spin,color,reim,SpinorType), offset), mask);
}

void readFVecSpinorXB(InstVector& ivector, const FVec& ret, string& base, string& offset, int spin, int color, int reim, string& mask)
{
#if defined(USE_LDUNPK) || ((SOALEN == VECLEN) && defined(USE_SHUFFLES))
    readUFVec(ivector, ret, new GatherAddress(new SpinorAddress(base,spin,color,reim,SpinorType), offset ), 0, mask);
#else
    gatherFVec(ivector, ret, new GatherAddress(new SpinorAddress(base,spin,color,reim,SpinorType), offset), mask);
#endif
}

void readFVecSpinorXF(InstVector& ivector, const FVec& ret, string& base, string& offset, int spin, int color, int reim, string& mask)
{
#if defined(USE_LDUNPK) || ((SOALEN == VECLEN) && defined(USE_SHUFFLES))
    readUFVec(ivector, ret, new GatherAddress(new SpinorAddress(base,spin,color,reim,SpinorType), offset ), 1, mask);
#else
    gatherFVec(ivector, ret, new GatherAddress(new SpinorAddress(base,spin,color,reim,SpinorType), offset), mask);
#endif
}

void writeFVecSpinor(InstVector& ivector, const FVec& ret, string& base, string& offset, int spin, int color, int reim, int isStreaming)
{
    writeFVecSpecialized(ivector, ret, new GatherAddress(new SpinorAddress(base,spin,color,reim,SpinorType), offset), isStreaming);
}

void readFVecGauge(InstVector& ivector, const FVec& ret, string& base, string& offset, int dir, int c1, int c2, int reim)
{
#ifndef USE_PACKED_GAUGES
    readFVecSpecialized(ivector, ret, new GatherAddress(new GaugeAddress(base,dir,c1,c2,reim,GaugeType), offset), string(""));
#else
    loadFVec(ivector, ret, new GaugeAddress(base,dir,c1,c2,reim,GaugeType), string(""));
#endif
}

void readFVecClovDiag(InstVector& ivector, const FVec& ret, string& base, string& offset, int block, int c)
{
#ifndef USE_PACKED_CLOVER
    readFVecSpecialized(ivector, ret, new GatherAddress(new ClovDiagAddress(base,block,c,CloverType), offset), string(""));
#else
    loadFVec(ivector, ret, new ClovDiagAddress(base,block,c,CloverType), string(""));
#endif
}

void readFVecClovOffDiag(InstVector& ivector, const FVec& ret, string& base, string& offset, int block, int c, int reim)
{
#ifndef USE_PACKED_CLOVER
    readFVecSpecialized(ivector, ret, new GatherAddress(new ClovOffDiagAddress(base,block,c,reim,CloverType), offset), string(""));
#else
    loadFVec(ivector, ret, new ClovOffDiagAddress(base,block,c,reim,CloverType), string(""));
#endif
}


void LoadSpinorElement(InstVector& ivector, const FVec& ret, string& base, string& offsets, int spin, int col, int reim, bool isFace, string mask, int dir)
{
    void (*readFunc)(InstVector& ivector, const FVec& ret, string& base, string& offset, int spin, int color, int reim, string& mask);

    readFunc = readFVecSpinor;

    if(dir == 0) {
        readFunc = readFVecSpinorXB;
    }

    if(dir == 1) {
        readFunc = readFVecSpinorXF;
    }

    readFunc(ivector, ret, base, offsets, spin, col, reim, mask);
}

void LoadFullSpinor(InstVector& ivector, const FVec ret[4][3][2], string& base, string& offsets, string mask)
{
    for(int col=0; col < 3; col++) {
        for(int spin=0; spin < 4; spin++) {
            LoadSpinorElement(ivector, ret[spin][col][RE], base, offsets, spin, col, RE, false, mask, -1);
            LoadSpinorElement(ivector, ret[spin][col][IM], base, offsets, spin, col, IM, false, mask, -1);
        }
    }
}

void StoreFullSpinor(InstVector& ivector, const FVec ret[4][3][2], string& base, string& offsets, int isStreaming)
{
    int useShuffles = 0;

#ifndef ENABLE_STREAMING_STORES
    isStreaming = 0;
#endif

#if defined(USE_SHUFFLES)
    useShuffles = 1;
#endif

    if(isStreaming) {
        useShuffles = 1;
    }

#if SOALEN == VECLEN
    useShuffles = 0;
#endif

    extern FVec tmp[];
    FVec in[24] = {
        ret[0][0][0], ret[0][0][1], ret[1][0][0], ret[1][0][1], ret[2][0][0], ret[2][0][1], ret[3][0][0], ret[3][0][1],
        ret[0][1][0], ret[0][1][1], ret[1][1][0], ret[1][1][1], ret[2][1][0], ret[2][1][1], ret[3][1][0], ret[3][1][1],
        ret[0][2][0], ret[0][2][1], ret[1][2][0], ret[1][2][1], ret[2][2][0], ret[2][2][1], ret[3][2][0], ret[3][2][1]
    };

    if(useShuffles) {
        GatherAddress *ga = new GatherAddress(new SpinorAddress(base,0,0,0,SpinorType), offsets );

        for(int i = 0; i < ((24*SOALEN)/VECLEN); i++) {
            transpose(ivector, &tmp[0], &in[i*(VECLEN/SOALEN)], SOALEN);

            for(int j = 0; j < (VECLEN/SOALEN); j++) {
                storeFVec(ivector, tmp[j], new AddressImm(ga->getAddr(j*SOALEN), i*VECLEN), isStreaming);
            }
        }
    }
    else {
        for(int spin=0; spin < 4; spin++) {
            for(int col=0; col < 3; col++) {
                writeFVecSpinor(ivector,  ret[spin][col][RE], base, offsets, spin,col,RE,isStreaming);
                writeFVecSpinor(ivector,  ret[spin][col][IM], base, offsets, spin,col,IM,isStreaming);
            }
        }
    }
}

void StreamFullSpinor(InstVector& ivector, const FVec ret[4][3][2], string& base, string& offsets)
{
    StoreFullSpinor(ivector, ret, base, offsets, 1);
}

void PackHalfSpinorElement(InstVector& ivector, const FVec& ret, string& base, int packoffset, int dir, string mask)
{
    if(dir >= 4 || (dir >= 2 && SOALEN == VECLEN)) {
        storeFVec(ivector, ret, new AddressImm(new GenericAddress(base, SpinorType), packoffset), 1);
    }
    else if(dir == 2) {
        storeSOAFVec(ivector, ret, new AddressImm(new GenericAddress(base, SpinorType), packoffset), 0, SOALEN);
    }
    else if(dir == 3) {
        storeSOAFVec(ivector, ret, new AddressImm(new GenericAddress(base, SpinorType), packoffset), (VECLEN/SOALEN)-1, SOALEN);
    }
    else {
        int possibleMsk = 0;
        int soaMsk = (dir == 0 ? 1 : (1 << (SOALEN-1)));

        for(int i = 0; i < VECLEN; i += SOALEN) {
            possibleMsk |= (soaMsk << i);
        }

        packFVec(ivector, ret, new AddressImm(new GenericAddress(base, SpinorType), packoffset), mask, possibleMsk);
    }
}

void PackHalfSpinor(InstVector& ivector, const FVec ret[2][3][2], string& base, int dir, string mask)
{
    int nActiveLanes = VECLEN;
    int ind = 0;

    if(dir == 0 || dir == 1) {
        nActiveLanes = (VECLEN/SOALEN+1)/2;    // +1 to avoid making it 0 when SOALEN=VECLEN
    }

    if(dir == 2 || dir == 3) {
        nActiveLanes = SOALEN;
    }

    for(int spin=0; spin < 2; spin++) {
        for(int col=0; col < 3; col++) {
            PackHalfSpinorElement(ivector, ret[spin][col][RE], base, (ind+0)*nActiveLanes, dir, mask);
            PackHalfSpinorElement(ivector, ret[spin][col][IM], base, (ind+1)*nActiveLanes, dir, mask);
            ind += 2;
        }
    }
}

void UnpackHalfSpinorElement(InstVector& ivector, const FVec& ret, string& base, int unpackoffset, int dir, string mask)
{
    if(dir >= 4 || (dir >= 2 && SOALEN == VECLEN)) {
        loadFVec(ivector, ret, new AddressImm(new GenericAddress(base, SpinorType), unpackoffset), string(""));
    }
    else if(dir == 2) {
        loadSOAFVec(ivector, ret, new AddressImm(new GenericAddress(base, SpinorType), unpackoffset), 0, SOALEN, string(""));
    }
    else if(dir == 3) {
        loadSOAFVec(ivector, ret, new AddressImm(new GenericAddress(base, SpinorType), unpackoffset), (VECLEN/SOALEN)-1, SOALEN, string(""));
    }
    else {
        int possibleMsk = 0;
        int soaMsk = (dir == 0 ? 1 : (1 << (SOALEN-1)));

        for(int i = 0; i < VECLEN; i += SOALEN) {
            possibleMsk |= (soaMsk << i);
        }

        unpackFVec(ivector, ret, new AddressImm(new GenericAddress(base, SpinorType), unpackoffset), mask, possibleMsk);
    }
}

void UnpackHalfSpinor(InstVector& ivector, const FVec ret[2][3][2], string& base, int dir, string mask)
{
    int nActiveLanes = VECLEN;
    int ind = 0;

    if(dir == 0 || dir == 1) {
        nActiveLanes = (VECLEN/SOALEN+1)/2;    // +1 to avoid making it 0 when SOALEN=VECLEN
    }

    if(dir == 2 || dir == 3) {
        nActiveLanes = SOALEN;
    }

    for(int spin=0; spin < 2; spin++) {
        for(int col=0; col < 3; col++) {
            UnpackHalfSpinorElement(ivector, ret[spin][col][RE], base, (ind+0)*nActiveLanes, dir, mask);
            UnpackHalfSpinorElement(ivector, ret[spin][col][IM], base, (ind+1)*nActiveLanes, dir, mask);
            ind += 2;
        }
    }
}

void LoadFullGaugeDir(InstVector& ivector, const FVec ret[3][3][2], string& base, string& offsets, int dir, bool compress12)
{
    int nrows = 3;

    if(compress12 == true) {
        nrows = 2;
    }

    for(int c1=0; c1 < nrows; c1++) {
        for(int c2=0; c2 < 3; c2++) {
            readFVecGauge(ivector, ret[c1][c2][RE], base, offsets, dir, c1, c2, RE);
            readFVecGauge(ivector, ret[c1][c2][IM], base, offsets, dir, c1, c2, IM);
        }
    }
}

void LoadFullCloverBlock(InstVector& ivector, const FVec diag[6], const FVec off_diag[15][2], string& base, string& offsets, int block)
{
    for(int c=0; c < 6; c++) {
        readFVecClovDiag(ivector, diag[c], base, offsets, block, c);
    }

    for(int c=0; c < 15; c++) {
        readFVecClovOffDiag(ivector, off_diag[c][RE], base, offsets, block, c, RE);
        readFVecClovOffDiag(ivector, off_diag[c][IM], base, offsets, block, c, IM);
    }
}


// Prefetches

void PrefetchL1Specialized(InstVector& ivector, GatherAddress *a, int type = 0, int dir = -1)
//void PrefetchL1Specialized(InstVector& ivector, GatherAddress *a, int type, int dir)
{
#ifndef NO_GPREF_L1
    gatherPrefetchL1(ivector, a, type);
#else

    for(int i = 0; i < VECLEN; i += SOALEN) {
        prefetchL1(ivector, a->getAddr(i), type);
    }

#endif
}

void PrefetchL2Specialized(InstVector& ivector, GatherAddress *a, int type = 0, int dir = -1)
//void PrefetchL2Specialized(InstVector& ivector, GatherAddress *a, int type, int dir)
{
#ifndef NO_GPREF_L2
    gatherPrefetchL2(ivector, a, type);
#else

    for(int i = 0; i < VECLEN; i += SOALEN) {
        prefetchL2(ivector, a->getAddr(i), type);
    }

#endif
}

void prefetchL1SpinorIn(InstVector& ivector, string base, string off, int imm, int dir, int type)
{
#ifdef PREF_L1_SPINOR_IN
    PrefetchL1Specialized(ivector, new GatherAddress(new AddressImm(new SpinorAddress(base,0,0,0,SpinorType), imm), off ), type, dir);
#endif
}

void prefetchL2SpinorIn(InstVector& ivector, string base, string off, const string& pref_dist, int imm, int dir)
{
#ifdef PREF_L2_SPINOR_IN
    PrefetchL2Specialized(ivector, new GatherAddress(new AddressImm(new AddressOffset(new SpinorAddress(base,0,0,0,SpinorType), pref_dist), imm), off ));
#endif
}

void prefetchL1SpinorOut(InstVector& ivector, string base, string off, int imm)
{
#ifdef PREF_L1_SPINOR_OUT
    PrefetchL1Specialized(ivector, new GatherAddress(new AddressImm(new SpinorAddress(base,0,0,0,SpinorType), imm), off ), 3 /*NT & EX */);
#endif
}

void prefetchL2SpinorOut(InstVector& ivector, string base, string off, const string& pref_dist, int imm)
{
#ifdef PREF_L2_SPINOR_OUT
    PrefetchL2Specialized(ivector, new GatherAddress(new AddressImm(new AddressOffset(new SpinorAddress(base,0,0,0,SpinorType), pref_dist), imm), off ), 3 /*NT & EX */);
#endif
}

void prefetchL1GuageDirIn(InstVector& ivector, string base, string off, int dir, int imm)
{
#ifdef PREF_L1_GAUGE
#ifndef USE_PACKED_GAUGES
    PrefetchL1Specialized(ivector, new GatherAddress(new AddressImm(new GaugeAddress(base,dir,0,0,0,GaugeType), imm), off ), 0);
#else
    prefetchL1(ivector, new AddressImm(new GaugeAddress(base,dir,0,0,0,GaugeType), imm), 0);
#endif
#endif
}

void prefetchL2GuageDirIn(InstVector& ivector, string base, string off, int dir, const string& pref_dist, int imm)
{
#ifdef PREF_L2_GAUGE
#ifndef USE_PACKED_GAUGES
    PrefetchL2Specialized(ivector, new GatherAddress(new AddressImm(new AddressOffset(new GaugeAddress(base,dir,0,0,0,GaugeType), pref_dist), imm), off ), 0);
#else
    prefetchL2(ivector, new AddressImm(new AddressOffset(new GaugeAddress(base,dir,0,0,0,GaugeType), pref_dist), imm), 0);
#endif
#endif
}

void prefetchL1CloverBlockIn(InstVector& ivector, string base, string off, int block, int imm)
{
#ifdef PREF_L1_CLOVER
#ifndef USE_PACKED_CLOVER
    PrefetchL1Specialized(ivector, new GatherAddress(new AddressImm(new ClovDiagAddress(base,block,0,CloverType), imm), off ), 0);
#else
    prefetchL1(ivector, new AddressImm(new ClovDiagAddress(base,block,0,CloverType), imm), 0);
#endif
#endif
}

void prefetchL2CloverBlockIn(InstVector& ivector, string base, string off, int block, const string& pref_dist, int imm)
{
#ifdef PREF_L2_CLOVER
#ifndef USE_PACKED_CLOVER
    PrefetchL2Specialized(ivector, new GatherAddress(new AddressImm(new AddressOffset(new ClovDiagAddress(base,block,0,CloverType), pref_dist), imm), off ), 0);
#else
    prefetchL2(ivector, new AddressImm(new AddressOffset(new ClovDiagAddress(base,block,0,CloverType), pref_dist), imm), 0);
#endif
#endif
}

void PrefetchL1FullSpinorDirIn(InstVector& ivector, const string& base, const string& off, int dir, int type)
{
    // for now we ignore direction but it can be used for specialization
    for(int i = 0; i < (24*SOALEN*sizeof(SpinorBaseType)+63)/64; i++) {
        prefetchL1SpinorIn(ivector, base, off, i*(64/sizeof(SpinorBaseType)), dir, type);
    }
}

void PrefetchL1FullSpinorOut(InstVector& ivector, const string& base, const string& off)
{
    for(int i = 0; i < (24*SOALEN*sizeof(SpinorBaseType)+63)/64; i++) {
        prefetchL1SpinorOut(ivector, base, off, i*(64/sizeof(SpinorBaseType)));
    }
}

void PrefetchL2FullSpinorDirIn(InstVector& ivector, const string& base, const string& off, const string& pref_dist, int dir)
{
    // for now we ignore direction but itcan be used for specialization
    for(int i = 0; i < (24*SOALEN*sizeof(SpinorBaseType)+63)/64; i++) {
        prefetchL2SpinorIn(ivector, base, off, pref_dist, i*(64/sizeof(SpinorBaseType)), dir);
    }
}

void PrefetchL2FullSpinorOut(InstVector& ivector, const string& base, const string& off, const string& pref_dist)
{
    for(int i = 0; i < (24*SOALEN*sizeof(SpinorBaseType)+63)/64; i++) {
        prefetchL2SpinorOut(ivector, base, off, pref_dist, i*(64/sizeof(SpinorBaseType)));
    }
}

void PrefetchL1FullGaugeDirIn(InstVector& ivector, const string& base, const string& off, int dir, bool compress12)
{
#ifndef USE_PACKED_GAUGES
    int nSites = SOALEN;
#else
    int nSites = VECLEN;
#endif
    int g_size=0;

    if( compress12 ) {
        g_size=2*3*2;
    }
    else {
        g_size=3*3*2;
    }

    for(int i = 0; i < ((g_size*nSites*sizeof(GaugeBaseType)+63)/64); i++) {
        prefetchL1GuageDirIn(ivector, base, off, dir, i*(64/sizeof(GaugeBaseType)));
    }
}

void PrefetchL2FullGaugeDirIn(InstVector& ivector, const string& base, const string& off, int dir, const string& pref_dist, bool compress12)
{
#ifndef USE_PACKED_GAUGES
    int nSites = SOALEN;
#else
    int nSites = VECLEN;
#endif
    int g_size=0;

    if( compress12 ) {
        g_size=2*3*2;
    }
    else {
        g_size=3*3*2;
    }

    for(int i = 0; i < ((g_size*nSites*sizeof(GaugeBaseType)+63)/64); i++) {
        prefetchL2GuageDirIn(ivector, base, off, dir, pref_dist, i*(64/sizeof(GaugeBaseType)));
    }
}

void PrefetchL2FullGaugeIn(InstVector& ivector, const string& base, const string& off, const string& pref_dist, bool compress12)
{
    for(int dir = 0; dir < 8; dir++) {
        PrefetchL2FullGaugeDirIn(ivector, base, off, dir, pref_dist, compress12);
    }
}

void PrefetchL1FullCloverBlockIn(InstVector& ivector, const string& base, const string& off, int block)
{
#ifndef USE_PACKED_CLOVER
    int nSites = SOALEN;
#else
    int nSites = VECLEN;
#endif

    for(int i = 0; i < ((36*nSites*sizeof(CloverBaseType)+63)/64); i++) {
        prefetchL1CloverBlockIn(ivector, base, off, block, i*(64/sizeof(CloverBaseType)));
    }
}

void PrefetchL2FullCloverIn(InstVector& ivector, const string& base, const string& off, const string& pref_dist)
{
#ifndef USE_PACKED_CLOVER
    int nSites = SOALEN;
#else
    int nSites = VECLEN;
#endif

    for(int i = 0; i < ((2*36*nSites*sizeof(CloverBaseType)+63)/64); i++) {
        prefetchL2CloverBlockIn(ivector, base, off, 0, pref_dist, i*(64/sizeof(CloverBaseType)));
    }
}

void PrefetchL1HalfSpinorDir(InstVector& ivector, const string& base, int dir, bool isPrefforWrite, int type)
{
    int nActiveLanes = VECLEN;
    int ind = 0;

    if(dir == 0 || dir == 1) {
        nActiveLanes = (VECLEN/SOALEN+1)/2;    // +1 to avoid making it 0 when SOALEN=VECLEN
    }

    if(dir == 2 || dir == 3) {
        nActiveLanes = SOALEN;
    }

#ifdef ENABLE_STREAMING_STORES

    if(nActiveLanes == VECLEN && isPrefforWrite) {
        return;
    }

#endif

    for(int i = 0; i < (12*nActiveLanes*sizeof(SpinorBaseType)+63)/64; i++) {
        prefetchL1(ivector, new AddressImm(new GenericAddress(base, SpinorType), i*(64/sizeof(SpinorBaseType))), type);
    }
}

void PrefetchL2HalfSpinorDir(InstVector& ivector, const string& base, const string& pref_dist, int dir, bool isPrefforWrite, int type)
{
    int nActiveLanes = VECLEN;
    int ind = 0;

    if(dir == 0 || dir == 1) {
        nActiveLanes = (VECLEN/SOALEN+1)/2;    // +1 to avoid making it 0 when SOALEN=VECLEN
    }

    if(dir == 2 || dir == 3) {
        nActiveLanes = SOALEN;
    }

#ifdef ENABLE_STREAMING_STORES

    if(nActiveLanes == VECLEN && isPrefforWrite) {
        return;
    }

#endif

    for(int i = 0; i < (12*nActiveLanes*sizeof(SpinorBaseType)+63)/64; i++) {
        prefetchL2(ivector, new AddressImm(new AddressOffset(new GenericAddress(base, SpinorType), pref_dist), i*(64/sizeof(SpinorBaseType))), type);
    }
}

