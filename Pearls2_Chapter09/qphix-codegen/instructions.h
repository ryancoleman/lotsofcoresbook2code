#ifndef __INSTRUCTIONS_H__
#define __INSTRUCTIONS_H__

#include <string>
#include <sstream>
#include <vector>

#ifndef PRECISION
#define PRECISION 1
#endif

#include "address_types.h"

using namespace std;

#define RE 0
#define IM 1

class FVec
{
public:
    FVec(const string& name_);
    FVec(const FVec& v_) : name(v_.getName()), type(v_.getType()) {}
    const string& getName() const
    {
        return name;
    }
    const string& getType() const
    {
        return type;
    }

private:
    const string name;
    const string type;
};

class Instruction
{
public:
    // string class return empty string
    virtual string serialize() const
    {
        return string("");
    }
    virtual bool hasAddress() const
    {
        return false;
    }
    virtual int numArithmeticInst() const
    {
        return 0;
    }
    virtual int numDeclarations() const
    {
        return 0;
    }
    virtual int numScopes() const
    {
        return 0;
    }
    virtual int numIfs() const
    {
        return 0;
    }
};

typedef vector<Instruction *> InstVector;

class BeginScope : public Instruction
{
public:
    string serialize() const
    {
        return "{";
    }
    int numScopes() const
    {
        return 1;
    }
};

class EndScope : public Instruction
{
public:
    string serialize() const
    {
        return "}";
    }
    int numScopes() const
    {
        return 1;
    }
};

class ElseStatement : public Instruction
{
public:
    string serialize() const
    {
        return "} else {";
    }
    int numScopes() const
    {
        return 1;
    }
};

class IfStringCond : public Instruction
{
public:
    IfStringCond(const string& condition_) : condition(condition_) {}
    string serialize() const
    {
        return " if ( "+condition+" ) { ";
    }
    int numIfs() const
    {
        return 1;
    }
private:
    const string condition;
};

class InlineCode : public Instruction
{
public:
    InlineCode(const string& code_) : code(code_) {}
    string serialize() const
    {
        return code;
    }
private:
    const string code;
};

enum MemRefType { LOAD_ALIGNED_VEC, LOAD_UNALIGNED_VEC, LOAD_MASKED_VEC, STORE_VEC, STREAM_VEC, STORE_MASKED_VEC, LOAD_NONVEC, L1_PREFETCH, NTA_PREFETCH, L2_PREFETCH, L1_EVICT, L2_EVICT, GATHER_VEC, SCATTER_VEC, GATHER_PREFETCH_L1, GATHER_PREFETCH_L2, GATHER_PREFETCH_NTA };

class MemRefInstruction : public Instruction
{
public:
    // Override virtual
    virtual bool hasAddress() const
    {
        return true;
    }
    virtual const Address* getAddress() const = 0;
    virtual MemRefType getType() const = 0;
    virtual bool hasGSAddress() const
    {
        return false;
    }
};

class DeclareFVec : public Instruction
{
public:
    DeclareFVec( const FVec &v_ ) : v(v_) { }
    string serialize() const;
    int numDeclarations() const
    {
        return 1;
    }

private:
    const FVec v;
};

class InitFVec : public Instruction
{
public:
    InitFVec( const FVec &v_ ) : v(v_) { }
    string serialize() const;
    int numDeclarations() const
    {
        return 0;
    }

private:
    const FVec v;
};

class DeclareMask : public Instruction
{
public:
    DeclareMask(string name_, string value_="") : name(name_), value(value_) { }
    string serialize() const;
    int numDeclarations() const
    {
        return (1);
    }
private:
    const string name;
    const string value;
};

class IntToMask : public Instruction
{
public:
    IntToMask(const string maskname, const string valname) : mask(maskname), value(valname) {}
    string serialize() const;
private:
    const string mask, value;
};

class DeclareOffsets : public Instruction
{
public:
    DeclareOffsets(string pname_, string vname_) : pname(pname_), vname(vname_) { }
    string serialize() const;
    int numDeclarations() const
    {
        return (1);
    }
private:
    const string vname;
    const string pname;
};

class IfAllOneCond : public Instruction
{
public:
    IfAllOneCond(const string& condition_) : condition(condition_) {}
    string serialize() const;
    int numIfs() const
    {
        return 1;
    }
private:
    const string condition;
};

class LoadFVec : public MemRefInstruction
{
public:
    LoadFVec( const FVec& v_, const Address* a_, const string mask_) : v(v_), a(a_), mask(mask_) {}
    string serialize() const;
    const Address* getAddress() const
    {
        return a;
    }
    MemRefType getType() const
    {
        return LOAD_ALIGNED_VEC;
    }
private:
    const FVec v;
    const Address* a;
    const string mask;

};

class StoreFVec : public MemRefInstruction
{
public:
    StoreFVec( const FVec& v_, const Address* a_, int isStreaming_) : v(v_), a(a_), isStreaming(isStreaming_) {}
    string serialize() const;
    const Address* getAddress() const
    {
        return a;
    }
    MemRefType getType() const
    {
        return STORE_VEC;
    }

private:
    const FVec v;
    const Address* a;
    int isStreaming;
};

class GatherFVec : public MemRefInstruction
{
public:
    GatherFVec( const FVec& v_, const GatherAddress* a_, const string mask_) : v(v_), a(a_), mask(mask_) {}
    string serialize() const;
    const Address* getAddress() const
    {
        return a;
    }
    MemRefType getType() const
    {
        return GATHER_VEC;
    }
    virtual bool hasGSAddress() const
    {
        return true;
    }
private:
    const FVec v;
    const GatherAddress* a;
    const string mask;
};

class ScatterFVec : public MemRefInstruction
{
public:
    ScatterFVec( const FVec& v_, const GatherAddress* a_) : v(v_), a(a_) {}
    string serialize() const;
    const Address* getAddress() const
    {
        return a;
    }
    MemRefType getType() const
    {
        return SCATTER_VEC;
    }
    virtual bool hasGSAddress() const
    {
        return true;
    }
private:
    const FVec v;
    const GatherAddress* a;
};

class LoadBroadcast : public MemRefInstruction
{
public:
    LoadBroadcast( const FVec& v_, const Address* a_) : v(v_), a(a_) {}
    string serialize() const;
    const Address* getAddress() const
    {
        return a;
    }
    MemRefType getType() const
    {
        return LOAD_NONVEC;
    }
private:
    const FVec v;
    const Address* a;
};


class PrefetchL1 : public MemRefInstruction
{
public:
    PrefetchL1( const Address* a_, int type = 0);

    string serialize() const
    {
        ostringstream stream;
        stream << " _mm_prefetch((const char *)( " << a->serialize() << " ), " << hint << ");" << endl;
        return stream.str();
    }

    MemRefType getType() const
    {
        return L1_PREFETCH;
    }
    const Address* getAddress() const
    {
        return a;
    }

private:
    const Address* a;
    string hint;
};

class PrefetchL2 : public MemRefInstruction
{
public:
    PrefetchL2( const Address* a_, int type = 0);
    string serialize() const
    {
        ostringstream stream;
        stream << " _mm_prefetch((const char *)( " << a->serialize() << " ), " << hint << ");" << endl;
        return stream.str();
    }

    MemRefType getType() const
    {
        return L2_PREFETCH;
    }
    const Address* getAddress() const
    {
        return a;
    }

private:
    const Address* a;
    string hint;
};

class GatherPrefetchL1 : public MemRefInstruction
{
public:
    GatherPrefetchL1( const GatherAddress* a_, int type = 0);
    string serialize() const;
    const Address* getAddress() const
    {
        return a;
    }
    MemRefType getType() const
    {
        return GATHER_PREFETCH_L1;
    }
    virtual bool hasGSAddress() const
    {
        return true;
    }
private:
    const GatherAddress* a;
    string hint;
};

class GatherPrefetchL2 : public MemRefInstruction
{
public:
    GatherPrefetchL2( const GatherAddress* a_, int type = 0);
    string serialize() const;
    const Address* getAddress() const
    {
        return a;
    }
    MemRefType getType() const
    {
        return GATHER_PREFETCH_L2;
    }
    virtual bool hasGSAddress() const
    {
        return true;
    }
private:
    const GatherAddress* a;
    string hint;
};


// Arithmetic Instructions
class SetZero : public Instruction
{
public:
    SetZero( const FVec& ret_) : ret(ret_) {}
    string serialize() const;
    int numArithmeticInst() const
    {
        return 1;
    }
private:
    const FVec ret;
};

class Mul : public Instruction
{
public:
    Mul( const FVec& ret_, const FVec& a_, const FVec& b_, const string& mask_) : ret(ret_), a(a_), b(b_), mask(mask_) {}
    string serialize() const;
    int numArithmeticInst() const
    {
        return 1;
    }
private:
    const FVec ret;
    const FVec a;
    const FVec b;
    const string mask;
};

class FnMAdd : public Instruction
{
public:
    FnMAdd( const FVec& ret_, const FVec& a_, const FVec& b_, const FVec& c_, const string& mask_) : ret(ret_), a(a_), b(b_), c(c_), mask(mask_) {}
    string serialize() const;
    int numArithmeticInst() const
    {
        return 1;
    }
private:
    const FVec ret;
    const FVec a;
    const FVec b;
    const FVec c;
    const string mask;
};

class FMAdd : public Instruction
{
public:
    FMAdd( const FVec& ret_, const FVec& a_, const FVec& b_, const FVec& c_, const string& mask_) : ret(ret_), a(a_), b(b_), c(c_), mask(mask_) {}
    string serialize() const;
    int numArithmeticInst() const
    {
        return 1;
    }
private:
    const FVec ret;
    const FVec a;
    const FVec b;
    const FVec c;
    const string mask;
};

class Add : public Instruction
{
public:
    Add( const FVec& ret_, const FVec& a_, const FVec& b_, const string& mask_) : ret(ret_), a(a_), b(b_), mask(mask_) {}
    string serialize() const;
    int numArithmeticInst() const
    {
        return 1;
    }
private:
    const FVec ret;
    const FVec a;
    const FVec b;
    const string mask;
};

class Sub : public Instruction
{
public:
    Sub( const FVec& ret_, const FVec& a_, const FVec& b_, const string& mask_) : ret(ret_), a(a_), b(b_), mask(mask_) {}
    string serialize() const;
    int numArithmeticInst() const
    {
        return 1;
    }
private:
    const FVec ret;
    const FVec a;
    const FVec b;
    const string mask;
};

class MovFVec : public Instruction
{
public:
    MovFVec( const FVec& ret_, const FVec& a_, const string mask_) : ret(ret_), a(a_), mask(mask_) {}
    string serialize() const;
    int numArithmeticInst() const
    {
        return 0;
    }
private:
    const FVec ret;
    const FVec a;
    const string mask;
};


void loadSOAFVec(InstVector& ivector, const FVec& ret, const Address *a, int soanum, int soalen, string mask);
void storeSOAFVec(InstVector& ivector, const FVec& ret, const Address *a, int soanum, int soalen);

void loadSplitSOAFVec(InstVector& ivector, const FVec& ret, const Address *a1, const Address *a2, int soanum, int soalen, int forward, string mask);
void unpackFVec(InstVector& ivector, const FVec& ret, Address *a, string mask, int possibleMask);
void packFVec(InstVector& ivector, const FVec& ret, Address *a, string mask, int possibleMask);

void gatherFVec(InstVector& ivector, const FVec& ret, GatherAddress *a, string mask);
void scatterFVec(InstVector& ivector, const FVec& ret, GatherAddress *a);
void gatherPrefetchL1(InstVector& ivector, GatherAddress *a, int type = 0);
void gatherPrefetchL2(InstVector& ivector, GatherAddress *a, int type = 0);

void transpose(InstVector& ivector, const FVec r[], const FVec f[], int soalen);


inline void movFVec(InstVector& ivector, const FVec& ret, const FVec& a, string mask)
{
    ivector.push_back(new MovFVec(ret, a, mask));
}

inline void beginScope(InstVector& ivector)
{
    ivector.push_back(new BeginScope());
}

inline void endScope(InstVector& ivector)
{
    ivector.push_back(new EndScope());
}

inline void ifStatement(InstVector& ivector, string condition)
{
    ivector.push_back( new IfStringCond(condition));
}

inline void elseStatement(InstVector& ivector)
{
    ivector.push_back( new ElseStatement());
}

inline void ifAllOneStatement(InstVector& ivector, string condition)
{
    ivector.push_back( new IfAllOneCond(condition));
}

inline void inlineCode(InstVector& ivector, string code)
{
    ivector.push_back( new InlineCode(code));
}

inline void declareFVecFromFVec(InstVector& ivector, const FVec& v)
{
    ivector.push_back(new DeclareFVec(v));
}

inline FVec declareFVec(InstVector& ivector, const std::string name)
{
    FVec tmp(name);
    ivector.push_back(new DeclareFVec( tmp ));
    return tmp;
}

inline void initFVec(InstVector& ivector, const FVec& ret)
{
    ivector.push_back(new InitFVec(ret));
}

inline void setZero(InstVector& ivector, const FVec& ret)
{
    ivector.push_back(new SetZero(ret));
}

inline void mulFVec(InstVector& ivector, const FVec& ret, const FVec& a, const FVec& b, string mask = "")
{
    ivector.push_back(new Mul(ret, a, b, mask));
}

inline void addFVec(InstVector& ivector, const FVec& ret, const FVec& a, const FVec& b, string mask = "")
{
    ivector.push_back(new Add(ret, a, b, mask));
}

inline void subFVec(InstVector& ivector, const FVec& ret, const FVec& a, const FVec& b, string mask = "")
{
    ivector.push_back(new Sub(ret, a, b, mask));
}

inline void fmaddFVec(InstVector& ivector, const FVec& ret, const FVec& a, const FVec& b, const FVec& c, string mask = "")
{
    ivector.push_back(new FMAdd(ret, a, b, c, mask));
}

inline void fnmaddFVec(InstVector& ivector, const FVec& ret, const FVec& a, const FVec& b, const FVec& c, string mask = "")
{
    ivector.push_back(new FnMAdd(ret, a, b, c, mask));
}

inline void loadFVec(InstVector& ivector, const FVec& ret, const Address *a, string mask)
{
    ivector.push_back( new LoadFVec(ret, a, mask));
}

inline void storeFVec(InstVector& ivector, const FVec& ret, const Address *a, int isStreaming)
{
    ivector.push_back( new StoreFVec(ret, a, isStreaming));
}

inline void loadBroadcastScalar(InstVector& ivector, const FVec& ret, string scalar_name, int type = 0)
{
    ivector.push_back( new LoadBroadcast(ret, new AddressOfScalar(scalar_name, type)));
}

inline void prefetchL1(InstVector& ivector, const Address *a, int type = 0)
{
    ivector.push_back( new PrefetchL1(a, type));
}

inline void prefetchL2(InstVector& ivector, const Address *a, int type = 0)
{
    ivector.push_back( new PrefetchL2(a, type));
}

inline void declareMask(InstVector& ivector, const string name, const string value="")
{
    ivector.push_back(new DeclareMask(name, value));
}

inline void intToMask(InstVector& ivector, const string maskname, const string intname)
{
    ivector.push_back(new IntToMask(maskname, intname));
}

#endif
