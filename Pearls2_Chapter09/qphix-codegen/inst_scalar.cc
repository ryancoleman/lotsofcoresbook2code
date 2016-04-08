#include <stdlib.h>
#include <stdio.h>
#include "instructions.h"

#if VECLEN == 1
#if PRECISION == 1
#pragma message "Using Scalar Single Precision"
#define FVECTYPE "float"
#define ZERO "0.0f"
#else
#pragma message "Using Scalar Double Precision"
#define FVECTYPE "double"
#define ZERO "0.0"
#endif

const string fullIntMask("0x1");
const string fullMask("0x1");

FVec::FVec(const string& name_) : name(name_), type(FVECTYPE) {}

string DeclareFVec::serialize() const
{

#if 0
    return v.getType()+" "+v.getName()+ ";";
#else
    return v.getType()+" "+v.getName()+ " = " ZERO ";";
#endif
}

string InitFVec::serialize() const
{
#if 1
    return "";
#else
    return v.getName()+ " = " ZERO ";";
#endif
}

string DeclareMask::serialize() const
{
    ostringstream outbuf;

    if(value.empty()) {
        outbuf << "int " << name << ";" << endl;
    }
    else {
        outbuf << "int " << name << " = " << value << ";" << endl;
    }

    return outbuf.str();
}

string IntToMask::serialize() const
{
    ostringstream outbuf;
    outbuf << mask << " = " << value << ";" << endl;
    return outbuf.str();
}

string DeclareOffsets::serialize() const
{
    ostringstream outbuf;
    outbuf << "int " << vname << " = (*(" << pname << "));" << endl;
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
            buf << v.getName() << " = *(" << a->serialize() << ");" <<endl;
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
            buf << "*(" << a->serialize() << ") = " << v.getName() <<  ";" <<endl;
        }
        else {
            printf("Error: Half type not supported for SSE\n");
            exit(1);
        }
    }
    else {
        if(!a->isHalfType()) {
            buf << "*(" << a->serialize() << ") = " << v.getName() <<  ";" <<endl;
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
        buf << v.getName() << " = (*" << a->serialize() << ");" << endl;
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
    return  ret.getName()+" = " ZERO ";";
}

string Mul::serialize() const
{
    if(mask.empty()) {
        return  ret.getName()+" = "+a.getName()+" * "+b.getName()+";" ;
    }
    else {
        return  ret.getName()+ " = _mm_blendv_ps(" + ret.getName() + ", _mm_mul_ps( "+a.getName()+" , "+b.getName()+"), " + mask + ");" ;
    }
}

string FnMAdd::serialize() const
{
    if(mask.empty()) {
        return  ret.getName()+" = "+c.getName()+" - ("+a.getName()+" * "+b.getName()+");" ;
    }
    else {
        return  ret.getName()+" = _mm_blendv_ps(" + ret.getName() + ", _mm_sub_ps("+c.getName()+", _mm_mul_ps("+a.getName()+" , "+b.getName()+")), " + mask + ");" ;
    }
}

string FMAdd::serialize() const
{
    if(mask.empty()) {
        return  ret.getName()+" = (("+a.getName()+" * "+b.getName()+") + "+c.getName()+");" ;
    }
    else {
        return  ret.getName()+" = _mm_blendv_ps(" + ret.getName() + ", _mm_add_ps(_mm_mul_ps("+a.getName()+", "+b.getName()+"), "+c.getName()+"), " + mask + ");" ;
    }
}

string Add::serialize() const
{
    if(mask.empty()) {
        return  ret.getName()+" = ("+a.getName()+" + "+b.getName()+");" ;
    }
    else {
        return  ret.getName()+" = _mm_blendv_ps(" + ret.getName() + ", _mm_add_ps( "+a.getName()+" , "+b.getName()+"), " + mask + ");" ;
    }
}

string Sub::serialize() const
{
    if(mask.empty()) {
        return  ret.getName()+" = ("+a.getName()+" - "+b.getName()+");" ;
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
            buf << v.getName() << " =  (*(" << a1->serialize() << "));" << endl;
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

void loadSOAFVec(InstVector& ivector, const FVec& ret, const Address *a, int soanum, int soalen, string mask)
{
    if(soalen == 1) {
        ivector.push_back( new LoadFVec(ret, a, string("")));
    }
    else {
        printf("SOALEN = %d not supported\n", soalen);
        exit(1);
    }
}

void storeSOAFVec(InstVector& ivector, const FVec& ret, const Address *a, int soanum, int soalen)
{
    if(soalen == 1) {
        ivector.push_back( new StoreFVec(ret, a, 0));
    }
    else {
        printf("SOALEN = %d not supported\n", soalen);
        exit(1);
    }
}

void loadSplitSOAFVec(InstVector& ivector, const FVec& ret, const Address *a1, const Address *a2, int soanum, int soalen, int forward, string mask)
{
    ivector.push_back( new LoadFVec(ret, a1, string("")));
}

void unpackFVec(InstVector& ivector, const FVec& ret, Address *a, string mask, int possibleMask)
{
    ivector.push_back( new LoadFVec(ret, a, string("")));
}

void packFVec(InstVector& ivector, const FVec& ret, Address *a, string mask, int possibleMask)
{
    ivector.push_back( new StoreFVec(ret, a, 0));
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

void transpose1x1(InstVector& ivector, const FVec r[1], const FVec f[1])
{
    movFVec(ivector, r[0], f[0], string(""));
}

void transpose(InstVector& ivector, const FVec r[], const FVec f[], int soalen)
{
    switch (soalen) {
    case 1:
        transpose1x1(ivector, r, f);
        break;

    default:
        printf("SOALEN = %d Not Supported (only SOALEN = 1 supported)\n", soalen);
    }
}

#endif // PRECISION == 1
