#include <sstream>
#include <fstream>
#include <iostream>
#include <vector>
#include <typeinfo>
#include <string>

using namespace std;

#include "dslash.h"

#if PRECISION == 1
#if VECLEN == 16
#ifdef AVX512
std::string ARCH_NAME="avx512";
#else
std::string ARCH_NAME="mic";
#endif
#elif VECLEN == 8
#ifdef AVX2
std::string ARCH_NAME="avx2";
#else
std::string ARCH_NAME="avx";
#endif
#elif VECLEN == 4
std::string ARCH_NAME="sse";
#elif VECLEN == 1
std::string ARCH_NAME="scalar";
#endif
#elif PRECISION == 2
#if VECLEN == 8
#ifdef AVX512
std::string ARCH_NAME="avx512";
#else
std::string ARCH_NAME="mic";
#endif
#elif VECLEN == 4
#ifdef AVX2
std::string ARCH_NAME="avx2";
#else
std::string ARCH_NAME="avx";
#endif
#elif VECLEN == 2
std::string ARCH_NAME="sse";
#elif VECLEN == 1
std::string ARCH_NAME="scalar";
#endif
#endif //PRECISION

void mergeIvectorWithL2Prefetches(InstVector& ivector, InstVector& l2prefs);
void dumpIVector(InstVector& ivector, string filename);

string dirname[2] = {"back", "forw"};
string dimchar[4] = {"X", "Y", "Z", "T"};



string basenames[8] = {"xyBase", "xyBase", "xyBase", "xyBase", "zbBase", "zfBase", "tbBase", "tfBase"};
string offsnames[8] = {"xbOffs", "xfOffs", "ybOffs", "yfOffs", "offs",   "offs",   "offs",   "offs"  };
string beta_names[8] = {"coeff_s", "coeff_s", "coeff_s", "coeff_s", "coeff_s", "coeff_s", "coeff_t_b", "coeff_t_f"};

extern FVec beta_vec;
string beta_name("beta");
string alpha_name("alpha");
string outBase("oBase");
string outOffs("offs");
string gBase("gBase");
string gOffs("gOffs");
string chiBase("chiBase");
string chiOffs("offs");
string clBase("clBase");
string clOffs("gOffs");

// Defines which dimensions are involved in SIMD blocking
// Currently just X and Y
bool requireAllOneCheck[4] = {true, true, false, false};
//bool requireAllOneCheck[4] = {false, false, false, false};

void generateFacePackL2Prefetches(InstVector& ivector, int dir)
{
    PrefetchL2HalfSpinorDir(ivector, "outbuf", "hsprefdist", dir, true, 2 /* Ex*/);
    PrefetchL2FullSpinorDirIn(ivector, "xyBase", "offs", "si_prefdist");
}

void generateFaceUnpackL2Prefetches(InstVector& ivector, int dir, bool compress12, bool clover)
{
    PrefetchL2HalfSpinorDir(ivector, "inbuf", "hsprefdist", dir, false, 0 /* None*/);
    PrefetchL2FullGaugeDirIn(ivector, "gBase", "gOffs", dir, "gprefdist", compress12);

    if(clover)	{
        PrefetchL2FullCloverIn(ivector, "clBase", "gOffs", "clprefdist");
    }

    PrefetchL2FullSpinorDirIn(ivector, outBase, "offs", "soprefdist");
}


// Generate all L2 prefetches
void generateL2Prefetches(InstVector& ivector, bool compress12, bool chi, bool clover)
{
    PrefetchL2FullSpinorDirIn(ivector, "xyBase", "pfyOffs", "siprefdist1");
    //PrefetchL2FullSpinorDirIn(ivector, "pfBase1", "offs", "siprefdist1");
    PrefetchL2FullSpinorDirIn(ivector, "pfBase2", "offs", "siprefdist2");
    PrefetchL2FullSpinorDirIn(ivector, "pfBase3", "offs", "siprefdist3");
    PrefetchL2FullSpinorDirIn(ivector, "pfBase4", "offs", "siprefdist4");

    if(clover)	{
        PrefetchL2FullCloverIn(ivector, "clBase", "gOffs", "clprefdist");
    }

    if(chi) {
        PrefetchL2FullSpinorDirIn(ivector, "pfBaseChi", "offs", "chiprefdist");
    }

    PrefetchL2FullGaugeIn(ivector, "gBase", "gOffs", "gprefdist", compress12);
    PrefetchL2FullSpinorOut(ivector, outBase, "offs", "siprefdist4");
}



#ifdef SERIAL_SPIN
void dslash_body(InstVector& ivector, bool compress12, proj_ops *ops, recons_ops *rec_ops_bw, recons_ops *rec_ops_fw, FVec outspinor[4][3][2])
{
    for(int dim = 0; dim < 4; dim++) {
        for(int dir = 0; dir < 2; dir++) {
            int d = dim * 2 + dir;
            stringstream d_str;
            d_str << d;
            string mask;
            bool adjMul;
            recons_ops rec_op;

            if(dir == 0) {
                adjMul = true;
                rec_op = rec_ops_bw[dim];
            }
            else {
                adjMul = false;
                rec_op = rec_ops_fw[dim];
            }

            ifStatement(ivector,"accumulate[" + d_str.str() + "]");
            {
                declareFVecFromFVec(ivector, beta_vec);
                loadBroadcastScalar(ivector, beta_vec, beta_names[d], SpinorType);

#ifdef NO_HW_MASKING

                if(requireAllOneCheck[dim]) {
                    ifAllOneStatement(ivector,"accumulate[" + d_str.str() + "]");
                    {
                        for(int s = 0; s < 2; s++) {
                            project(ivector, basenames[d], offsnames[d], ops[d], false, mask, d, s);

                            if(s==0) {
                                loadGaugeDir(ivector, d, compress12);
                            }

                            matMultVec(ivector, adjMul, s);
                            recons_add(ivector, rec_op, outspinor, mask, s);
                        }
                    }
                    elseStatement(ivector);
                }

#endif

                if(requireAllOneCheck[dim]) {
                    mask = "accMask";
                    declareMask(ivector, mask);
                    intToMask(ivector, mask, "accumulate[" + d_str.str() + "]");
                }

                for(int s = 0; s < 2; s++) {
                    project(ivector, basenames[d], offsnames[d], ops[d], false, mask, d, s);

                    if(s==0) {
                        loadGaugeDir(ivector, d, compress12);
                    }

                    matMultVec(ivector, adjMul, s);
                    recons_add(ivector, rec_op, outspinor, mask, s);
                }

#ifdef NO_HW_MASKING

                if(requireAllOneCheck[dim]) {
                    endScope(ivector);
                }

#endif
            }
            endScope(ivector);
        }
    }
}
#else // NO SERIAL_SPIN
void dslash_body(InstVector& ivector, bool compress12, proj_ops *ops, recons_ops *rec_ops_bw, recons_ops *rec_ops_fw, FVec outspinor[4][3][2])
{
    for(int dim = 0; dim < 4; dim++) {
        for(int dir = 0; dir < 2; dir++) {
            int d = dim * 2 + dir;
            stringstream d_str;
            d_str << d;
            string mask;
            bool adjMul;
            recons_ops rec_op;

            adjMul = (dir == 0 ? true : false);
            rec_op = (dir == 0 ? rec_ops_bw[dim] : rec_ops_fw[dim]);

            ifStatement(ivector,"accumulate[" + d_str.str() + "]");
            {
                declareFVecFromFVec(ivector, beta_vec);
                loadBroadcastScalar(ivector, beta_vec, beta_names[d], SpinorType);

#ifdef NO_HW_MASKING

                if(requireAllOneCheck[dim]) {
                    ifAllOneStatement(ivector,"accumulate[" + d_str.str() + "]");
                    {
                        project(ivector, basenames[d], offsnames[d], ops[d], false, mask, d);
                        loadGaugeDir(ivector, d, compress12);
                        matMultVec(ivector, adjMul);
                        recons_add(ivector, rec_op, outspinor, mask);
                    }
                    elseStatement(ivector);
                }

#endif

                if(requireAllOneCheck[dim]) {
                    mask = "accMask";
                    declareMask(ivector, mask);
                    intToMask(ivector, mask, "accumulate[" + d_str.str() + "]");
                }

                project(ivector, basenames[d], offsnames[d], ops[d], false, mask, d);
                loadGaugeDir(ivector, d, compress12);
                matMultVec(ivector, adjMul);
                recons_add(ivector, rec_op, outspinor, mask);
#ifdef NO_HW_MASKING

                if(requireAllOneCheck[dim]) {
                    endScope(ivector);
                }

#endif
            }
            endScope(ivector);
        }
    }
}
#endif // SERIAL_SPIN



// need xyBase, and offs to specify input spinor
// need outbuf for output half spinor
void pack_face_vec(InstVector& ivector, FVec spinor[2][3][2], proj_ops proj[], int dir)
{
    std::string intMask, mask;

    // Check if this dir has mask argument
    if(requireAllOneCheck[dir/2]) {
        intMask = "mask";
        mask = "accMask";
        declareMask(ivector, mask);
        intToMask(ivector, mask, "mask");
    }

    std::string out("outbuf");
    PrefetchL1HalfSpinorDir(ivector, out, dir, true, 2 /*Exclusive*/);

    // We need to reverse direction of projection for our neighbor
    int fb = (dir % 2 == 0 ? 1 : -1);
    project(ivector, "xyBase","offs", proj[dir+fb], true, mask, dir);

    // This will write it to outbuf
    PackHalfSpinor(ivector, spinor, out, dir, intMask);
}

// need inbuf pointer to half spinor
// need gBase and goffs to point to gauge
// need obase and offs to point to spinor to scatter.
void recons_add_face_vec(InstVector& ivector, bool compress12, bool adjMul, recons_ops rops[], int dir, int dim, bool clover)
{

    std::string in("inbuf");
    std::string mask, intMask;

    extern FVec out_spinor[4][3][2];
    extern FVec dout_spinor[4][3][2];
    extern FVec b_spinor[2][3][2];

    int gauge_index = dim * 2 + dir;

    // Check if this dir has mask argument
    if(requireAllOneCheck[dim]) {
        intMask = "mask";
        mask = "accMask";
        declareMask(ivector, mask);
        intToMask(ivector, mask, "mask");
    }

    declareFVecFromFVec(ivector, beta_vec);
    loadBroadcastScalar(ivector, beta_vec, beta_name, SpinorType);

    FVec (*outspinor)[4][3][2];

    if(clover) {
        outspinor = &dout_spinor;
        zeroResult(ivector, (*outspinor)[0][0]);
    }
    else {
        outspinor = &out_spinor;
    }

    PrefetchL1HalfSpinorDir(ivector, in, dir, false, 0 /*None*/);
    // Gather in the partial result
    PrefetchL1FullSpinorDirIn(ivector, outBase, outOffs, -1);
    LoadFullSpinor(ivector, out_spinor, outBase, outOffs, "");

    // load b-from inbuf
    UnpackHalfSpinor(ivector, b_spinor, in, gauge_index, intMask);

    loadGaugeDir(ivector, gauge_index, compress12);
    matMultVec(ivector, adjMul);
    recons_add(ivector, rops[dim], *outspinor, mask);

    if(clover) {
        clover_term(ivector, *outspinor, true);
    }

    // scatter it out
    StoreFullSpinor(ivector, out_spinor, outBase, outOffs);
}


string getTypeName(size_t s)
{
    if(s == 2) {
        return "half";
    }
    else if(s == 4) {
        return "float";
    }
    else if(s == 8) {
        return "double";
    }
    else {
        return "Unknown";
    }
}

void generate_code(void)
{
    InstVector ivector;
    InstVector l2prefs;
    bool compress12;

    const string SpinorTypeName = getTypeName(sizeof(SpinorBaseType));
    const string GaugeTypeName = getTypeName(sizeof(GaugeBaseType));
    const string CloverTypeName = getTypeName(sizeof(CloverBaseType));

    if(SOALEN == VECLEN) {
        requireAllOneCheck[1] = false;
    }

#ifdef NO_MASKS

    for(int i = 0; i < 4; i++) {
        requireAllOneCheck[i] = false;
    }

#endif

    for(int isign=-1; isign<=1; isign+=2) {
        bool isPlus = (isign == 1 ? true : false);
        string plusminus = (isPlus) ? "plus" : "minus";

        for(int clov = 0; clov < 2; clov++) {
            bool clover = (clov == 1 ? true : false);
            string clov_prefix = (clover ? "clov_"+CloverTypeName+"_" : "");

            for(int num_components=12; num_components <=18; num_components+=6) {
                compress12 = ( num_components==12 );

                std::ostringstream filename;
                filename << "./"<<ARCH_NAME<<"/" << clov_prefix << "dslash_"<<plusminus<<"_body_" << SpinorTypeName << "_" << GaugeTypeName << "_v"<< VECLEN <<"_s"<<SOALEN<<"_"<< num_components;
                l2prefs.resize(0);
                generateL2Prefetches(l2prefs,compress12, false, clover);
                //dumpIVector(l2prefs, "prefetches.out");

                // Dslash Plus
                cout << "GENERATING dslash_"<<plusminus<<"_vec body" << endl;
                // Flush instruction list
                ivector.resize(0);

                // Generate instructions
                dslash_plain_body(ivector,compress12,clover,isPlus);
                mergeIvectorWithL2Prefetches(ivector, l2prefs);
                dumpIVector(ivector,filename.str());

                filename.str("");
                filename.clear();
                filename << "./"<<ARCH_NAME<<"/" << clov_prefix << "dslash_achimbdpsi_"<<plusminus<<"_body_" << SpinorTypeName << "_" << GaugeTypeName << "_v"<< VECLEN <<"_s"<<SOALEN<<"_"<< num_components;

                l2prefs.resize(0);
                generateL2Prefetches(l2prefs,compress12, true, clover);
                //dumpIVector(l2prefs, "prefetches.out");
                cout << "GENERATING dslash_achimbdpsi_"<<plusminus<<"_vec body" << endl;
                // Flush instruction list
                ivector.resize(0);

                // Generate instructions
                dslash_achimbdpsi_body(ivector,compress12,clover,isPlus);
                mergeIvectorWithL2Prefetches(ivector, l2prefs);
                dumpIVector(ivector,filename.str());

                for(int dir = 0; dir < 2; dir++) {
                    for(int dim = 0; dim < 4; dim++) {
                        std::ostringstream filename;
                        filename << "./"<<ARCH_NAME<<"/" << clov_prefix << "dslash_face_unpack_from_"<<dirname[dir]<<"_"<<dimchar[dim]<<"_" << plusminus <<"_" << SpinorTypeName << "_" << GaugeTypeName << "_v"<< VECLEN <<"_s"<<SOALEN<<"_"<<num_components;
                        cout << "GENERATING face unpack file " << filename.str() << endl;
                        ivector.resize(0);
                        recons_add_face_from_dir_dim_vec(ivector, compress12,isPlus, dir, dim, clover);
                        l2prefs.resize(0);
                        generateFaceUnpackL2Prefetches(l2prefs, 2*dim+dir, compress12, clover);
                        mergeIvectorWithL2Prefetches(ivector, l2prefs);
                        dumpIVector(ivector,filename.str());
                    }
                }
            }
        }

        for(int dir = 0; dir < 2; dir++) {
            for(int dim = 0; dim < 4; dim++) {
                std::ostringstream filename;
                filename << "./"<<ARCH_NAME<<"/dslash_face_pack_to_"<<dirname[dir]<<"_"<<dimchar[dim]<<"_"<<plusminus<<"_" << SpinorTypeName << "_" << GaugeTypeName << "_v"<< VECLEN <<"_s"<<SOALEN;
                cout << "GENERATING face pack file " << filename.str() << endl;
                l2prefs.resize(0);
                generateFacePackL2Prefetches(l2prefs, 2*dim+dir);
                ivector.resize(0);
                pack_face_to_dir_dim_vec(ivector,isPlus,dir,dim);
                mergeIvectorWithL2Prefetches(ivector, l2prefs);
                dumpIVector(ivector,filename.str());
            }
        }
    }

    data_types<float,VECLEN,SOALEN,true>::Gauge cmped;
    data_types<float,VECLEN,SOALEN,false>::Gauge uncmped;

    cout << "Compressed Gauge size is " << sizeof(cmped) << endl;
    cout << "Uncompressed Gauge size is " << sizeof(uncmped) << endl;

}
