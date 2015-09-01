#include <stdlib.h>
#include <stdio.h>
#include "instructions.h"

#if PRECISION == 1 && VECLEN == 4
#pragma message "Using SSE Single Precision"

#define FVECTYPE "__m128"

const string fullIntMask("0xF");
const string fullMask("0xF");

FVec::FVec(const string& name_) : name(name_), type(FVECTYPE) {}

string DeclareFVec::serialize() const
{
#if 0
    return v.getType()+" "+v.getName()+ " = _mm_undefined_ps();";
#else
    return v.getType()+" "+v.getName()+ " = _mm_setzero_ps();";
#endif
}

string InitFVec::serialize() const
{
#if 0
    return v.getName()+ " = _mm_undefined_ps();";
#else
    return v.getName()+ " = _mm_setzero_ps();";
#endif
}

string DeclareMask::serialize() const
{
    ostringstream outbuf;

    if(value.empty()) {
        outbuf << "__m128 " << name << ";" << endl;
    }
    else {
        outbuf << "__m128 " << name << " = " << value << ";" << endl;
    }

    return outbuf.str();
}

string IntToMask::serialize() const
{
    ostringstream outbuf;
    outbuf << mask << " = _mm_int2mask_ps(" << value << ");" << endl;
    return outbuf.str();
}

string DeclareOffsets::serialize() const
{
    ostringstream outbuf;
    outbuf << "__m128i " << vname << " = _mm_load_si128((__m128i const *)" << pname << ");" << endl;
    return outbuf.str();
}

string IfAllOneCond::serialize() const
{
    return " if ((" + condition + " & " + fullIntMask + ") == " + fullIntMask + ") { ";
}

string LoadFVec::serialize() const
{
    std::ostringstream buf;

    if(mask.empty()) {
        if(!a->isHalfType()) {
            buf << v.getName() << " = _mm_load_ps(" << a->serialize() << ");" <<endl;
        }
        else {
            printf("Error: Half type not supported for SSE\n");
            exit(1);
        }
    }
    else {
        if(!a->isHalfType()) {
            printf("Error: Masked load not supported for SSE\n");
            exit(1);
        }
        else {
            printf("Error: Half type not supported for SSE\n");
            exit(1);
        }
    }

    return buf.str();

}

string StoreFVec::serialize() const
{
    ostringstream buf;
    int streaming = isStreaming;

    if(streaming) {
        if(!a->isHalfType()) {
            buf << "_mm_stream_ps(" << a->serialize() << ", " << v.getName() <<  ");" <<endl;
        }
        else {
            printf("Error: Half type not supported for SSE\n");
            exit(1);
        }
    }
    else {
        if(!a->isHalfType()) {
            buf << "_mm_store_ps(" << a->serialize() << ", " << v.getName() <<  ");" <<endl;
        }
        else {
            printf("Error: Half type not supported for SSE\n");
            exit(1);
        }
    }

    return buf.str();
}

string GatherFVec::serialize() const
{
    std::ostringstream buf;

    printf("Error: Gather not supported for SSE\n");
    exit(1);
    return buf.str();

}

string ScatterFVec::serialize() const
{
    std::ostringstream buf;
    printf("Error: Scatter not supported for SSE\n");
    exit(1);
    return buf.str();
}

string LoadBroadcast::serialize() const
{
    std::ostringstream buf;

    if(!a->isHalfType()) {
        buf << v.getName() << " = _mm_set1_ps(*" << a->serialize() << ");" << endl;
    }
    else {
        printf("Error: Half type not supported for SSE\n");
        exit(1);
    }

    return buf.str();
}

PrefetchL1::PrefetchL1( const Address* a_, int type) : a(a_)
{
    // Type: 0 - none, 1 - NT, 2 - Ex, 3 - NT+Ex
    switch (type) {
    case 0:
        hint = "_MM_HINT_T0";
        break;

    case 1:
        hint = "_MM_HINT_NTA";
        break;

    case 2:
        hint = "_MM_HINT_T0";
        break;

    case 3:
        hint = "_MM_HINT_NTA";
        break;
    }
}

PrefetchL2::PrefetchL2( const Address* a_, int type) : a(a_)
{
    // Type: 0 - none, 1 - NT, 2 - Ex, 3 - NT+Ex
    switch(type) {
    case 0:
        hint = "_MM_HINT_T1";
        break;

    case 1:
        hint = "_MM_HINT_T2";
        break;

    case 2:
        hint = "_MM_HINT_T1";
        break;

    case 3:
        hint = "_MM_HINT_T2";
        break;
    }
}

GatherPrefetchL1::GatherPrefetchL1( const GatherAddress* a_, int type) : a(a_)
{
    // Type: 0 - none, 1 - NT, 2 - Ex, 3 - NT+Ex
    switch (type) {
    case 0:
        hint = "_MM_HINT_T0";
        break;

    case 1:
        hint = "_MM_HINT_NTA";
        break;
//		case 2: hint = "_MM_HINT_ET0"; break;
//		case 3: hint = "_MM_HINT_ENTA"; break;
    }
}
string GatherPrefetchL1::serialize() const
{
    std::ostringstream buf;
    buf << "#error \"Gather Prefetch is not supported in SSE\"" <<endl;
    printf("Gather Prefetch is not supported in SSE\n");
    exit(1);
    return buf.str();
}

GatherPrefetchL2::GatherPrefetchL2( const GatherAddress* a_, int type) : a(a_)
{
    // Type: 0 - none, 1 - NT, 2 - Ex, 3 - NT+Ex
    switch(type) {
    case 0:
        hint = "_MM_HINT_T1";
        break;

    case 1:
        hint = "_MM_HINT_T2";
        break;
//		case 2: hint = "_MM_HINT_ET1"; break;
//		case 3: hint = "_MM_HINT_ET2"; break;
    }
}
string GatherPrefetchL2::serialize() const
{
    std::ostringstream buf;
    buf << "#error \"Gather Prefetch is not supported in SSE\"" <<endl;
    printf("Gather Prefetch is not supported in SSE\n");
    exit(1);
    return buf.str();

}

string SetZero::serialize() const
{
    return  ret.getName()+" = _mm_setzero_ps(); ";
}

string Mul::serialize() const
{
    if(mask.empty()) {
        return  ret.getName()+" = _mm_mul_ps( "+a.getName()+" , "+b.getName()+" );" ;
    }
    else {
        return  ret.getName()+ " = _mm_blendv_ps(" + ret.getName() + ", _mm_mul_ps( "+a.getName()+" , "+b.getName()+"), " + mask + ");" ;
    }
}

string FnMAdd::serialize() const
{
    if(mask.empty()) {
        return  ret.getName()+" = _mm_sub_ps("+c.getName()+", _mm_mul_ps("+a.getName()+" , "+b.getName()+"));" ;
    }
    else {
        return  ret.getName()+" = _mm_blendv_ps(" + ret.getName() + ", _mm_sub_ps("+c.getName()+", _mm_mul_ps("+a.getName()+" , "+b.getName()+")), " + mask + ");" ;
    }
}

string FMAdd::serialize() const
{
    if(mask.empty()) {
        return  ret.getName()+" = _mm_add_ps(_mm_mul_ps("+a.getName()+", "+b.getName()+"), "+c.getName()+" );" ;
    }
    else {
        return  ret.getName()+" = _mm_blendv_ps(" + ret.getName() + ", _mm_add_ps(_mm_mul_ps("+a.getName()+", "+b.getName()+"), "+c.getName()+"), " + mask + ");" ;
    }
}

string Add::serialize() const
{
    if(mask.empty()) {
        return  ret.getName()+" = _mm_add_ps( "+a.getName()+" , "+b.getName()+" );" ;
    }
    else {
        return  ret.getName()+" = _mm_blendv_ps(" + ret.getName() + ", _mm_add_ps( "+a.getName()+" , "+b.getName()+"), " + mask + ");" ;
    }
}

string Sub::serialize() const
{
    if(mask.empty()) {
        return  ret.getName()+" = _mm_sub_ps( "+a.getName()+" , "+b.getName()+" );" ;
    }
    else {
        return  ret.getName()+" = _mm_blendv_ps(" + ret.getName() + ", _mm_sub_ps( "+a.getName()+" , "+b.getName()+"), " + mask + ");" ;
    }
}

string MovFVec::serialize() const
{
    if(mask.empty()) {
        return  ret.getName()+" = " + a.getName()+";" ;
    }
    else {
        return ret.getName()+" = _mm_blendv_ps(" + ret.getName() + ", "+a.getName()+", " + mask + ");" ;
    }
}

class LoadSplitSOAFVec : public MemRefInstruction
{
public:
    LoadSplitSOAFVec( const FVec& v_, const Address* a1_, const Address* a2_, const int soanum_, const int soalen_, int forward_) : v(v_), a1(a1_), a2(a2_), soanum(soanum_), soalen(soalen_), forward(forward_) {}
    string serialize() const
    {
        std::ostringstream buf;

        if(!a1->isHalfType()) {
            if(forward) {
                if(soalen == 4) {
                    buf << v.getName() << " =  _mm_blend_ps(_mm_loadu_ps(" << a1->serialize() << "), _mm_set1_ps(*" <<
                        a2->serialize() << "), " << (1 << (soalen-1)) << ");" << endl;
                }
            }
            else {
                if(soalen == 4) {
                    buf << v.getName() << " =  _mm_blend_ps(_mm_loadu_ps((" << a2->serialize() << ")-1), _mm_set1_ps(*" << a1->serialize() << "), 1);" << endl;
                }
            }
        }
        else {
            printf("Error: Half type not supported for SSE\n");
            exit(1);
        }

        return buf.str();
    }
    const Address* getAddress() const
    {
        return a1;
    }
    MemRefType getType() const
    {
        return LOAD_MASKED_VEC;
    }
private:
    const FVec v;
    const Address* a1;
    const Address* a2;
    const int soalen, soanum;
    const int forward;
};

class CondInsertFVecElement : public MemRefInstruction
{
public:
    CondInsertFVecElement( const FVec& v_, const Address* a_, const string mask_, int pos_, bool skipCond_) : v(v_), a(a_), mask(mask_), pos(pos_), skipCond(skipCond_) {}
    string serialize() const
    {
        std::ostringstream buf;

        if(!skipCond) {
            buf << "if(" << mask << " & " << (1 << pos) << ") ";
        }

        if(!a->isHalfType()) {
            buf << v.getName() << " = _mm_blend_ps(" << v.getName() << ", _mm_set1_ps(*" << a->serialize() << "), " << (1<<pos) << ");" << endl;
        }
        else {
            printf("Error: Half type not supported for SSE\n");
            exit(1);
        }

        return buf.str();
    }
    const Address* getAddress() const
    {
        return a;
    }
    MemRefType getType() const
    {
        return STORE_MASKED_VEC;
    }
private:
    const FVec v;
    const Address* a;
    const string mask;
    const int pos;
    const bool skipCond;
};

class CondExtractFVecElement : public MemRefInstruction
{
public:
    CondExtractFVecElement( const FVec& v_, const Address* a_, const string mask_, int pos_, bool skipCond_) : v(v_), a(a_), mask(mask_), pos(pos_), skipCond(skipCond_) {}
    string serialize() const
    {
        std::ostringstream buf;

        if(!skipCond) {
            buf << "if(" << mask << " & " << (1 << pos) << ") ";
        }

        if(!a->isHalfType()) {
            buf << "((int*)" << a->serialize() << ")[0] = _mm_extract_ps(" << v.getName() << ", " << (pos % 4) << ");" << endl;
        }
        else {
            printf("Error: Half type not supported for SSE\n");
            exit(1);
        }

        return buf.str();
    }
    const Address* getAddress() const
    {
        return a;
    }
    MemRefType getType() const
    {
        return STORE_MASKED_VEC;
    }
private:
    const FVec v;
    const Address* a;
    const string mask;
    const int pos;
    const bool skipCond;
};

void loadSOAFVec(InstVector& ivector, const FVec& ret, const Address *a, int soanum, int soalen, string mask)
{
    if(soalen == 4) {
        ivector.push_back( new LoadFVec(ret, a, string("")));
    }
    else {
        printf("SOALEN = %d not supported\n", soalen);
        exit(1);
    }
}

void storeSOAFVec(InstVector& ivector, const FVec& ret, const Address *a, int soanum, int soalen)
{
    if(soalen == 4) {
        ivector.push_back( new StoreFVec(ret, a, 0));
    }
    else {
        printf("SOALEN = %d not supported\n", soalen);
        exit(1);
    }
}

void loadSplitSOAFVec(InstVector& ivector, const FVec& ret, const Address *a1, const Address *a2, int soanum, int soalen, int forward, string mask)
{
    ivector.push_back( new LoadSplitSOAFVec(ret, a1, a2, soanum, soalen, forward));
}

void unpackFVec(InstVector& ivector, const FVec& ret, Address *a, string mask, int possibleMask)
{
    int pos = 0, nBits = 0;

    for(int i = 0; i < 4; i++) if(possibleMask & (1 << i)) {
            nBits++;
        }

    for(int i = 0; i < 4; i++)
        if(possibleMask & (1 << i)) {
            ivector.push_back( new CondInsertFVecElement(ret, new AddressImm(a, pos), mask, i, nBits==1));
            //pos++;
        }
}

void packFVec(InstVector& ivector, const FVec& ret, Address *a, string mask, int possibleMask)
{
    int pos = 0, nBits = 0;

    for(int i = 0; i < 4; i++) if(possibleMask & (1 << i)) {
            nBits++;
        }

    for(int i = 0; i < 4; i++)
        if(possibleMask & (1 << i)) {
            ivector.push_back( new CondExtractFVecElement(ret, new AddressImm(a, pos), mask, i, nBits==1));
            //pos++;
        }
}

void gatherFVec(InstVector& ivector, const FVec& ret, GatherAddress *a, string mask)
{
    ivector.push_back( new GatherFVec(ret, a, mask));
}

void scatterFVec(InstVector& ivector, const FVec& ret, GatherAddress *a)
{
    ivector.push_back( new ScatterFVec(ret, a));
}

void gatherPrefetchL1(InstVector& ivector, GatherAddress *a, int type)
{
    ivector.push_back( new GatherPrefetchL1(a, type));
}

void gatherPrefetchL2(InstVector& ivector, GatherAddress *a, int type)
{
    ivector.push_back( new GatherPrefetchL2(a, type));
}

void fperm32x4(InstVector& ivector, const FVec& ret, const FVec& a, const FVec& b, const int imm)
{
    //ivector.push_back(new Perm32x4(ret, a, b, imm));
    printf("Perm32x4() Not supported\n");
    exit(1);
}

void transpose1x1(InstVector& ivector, const FVec r[1], const FVec f[1])
{
    movFVec(ivector, r[0], f[0], string(""));
}

void transpose(InstVector& ivector, const FVec r[], const FVec f[], int soalen)
{
    switch (soalen) {
    case 4:
        transpose1x1(ivector, r, f);
        break;

    default:
        printf("SOALEN = %d Not Supported (only SOALEN = 4 supported)\n", soalen);
    }
}

#endif // PRECISION == 1
