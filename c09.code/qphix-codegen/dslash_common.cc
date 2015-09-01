#include <cstdio>
#include <cstdlib>
using namespace std;

#include "dslash.h"

extern string beta_names[8];
extern string alpha_name;
extern string outBase;
extern string outOffs;
extern string gBase;
extern string gOffs;
extern string chiBase;
extern string chiOffs;
extern string clBase;
extern string clOffs;


FVec b_S0_C0_RE("b_S0_C0_RE");
FVec b_S0_C0_IM("b_S0_C0_IM");
FVec b_S0_C1_RE("b_S0_C1_RE");
FVec b_S0_C1_IM("b_S0_C1_IM");
FVec b_S0_C2_RE("b_S0_C2_RE");
FVec b_S0_C2_IM("b_S0_C2_IM");
FVec b_S1_C0_RE("b_S1_C0_RE");
FVec b_S1_C0_IM("b_S1_C0_IM");
FVec b_S1_C1_RE("b_S1_C1_RE");
FVec b_S1_C1_IM("b_S1_C1_IM");
FVec b_S1_C2_RE("b_S1_C2_RE");
FVec b_S1_C2_IM("b_S1_C2_IM");

FVec b_spinor[2][3][2] = {
    { {b_S0_C0_RE, b_S0_C0_IM}, {b_S0_C1_RE, b_S0_C1_IM}, {b_S0_C2_RE, b_S0_C2_IM} },
    { {b_S1_C0_RE, b_S1_C0_IM}, {b_S1_C1_RE, b_S1_C1_IM}, {b_S1_C2_RE, b_S1_C2_IM} }
};

FVec ub_S0_C0_RE("ub_S0_C0_RE");
FVec ub_S0_C0_IM("ub_S0_C0_IM");
FVec ub_S0_C1_RE("ub_S0_C1_RE");
FVec ub_S0_C1_IM("ub_S0_C1_IM");
FVec ub_S0_C2_RE("ub_S0_C2_RE");
FVec ub_S0_C2_IM("ub_S0_C2_IM");
FVec ub_S1_C0_RE("ub_S1_C0_RE");
FVec ub_S1_C0_IM("ub_S1_C0_IM");
FVec ub_S1_C1_RE("ub_S1_C1_RE");
FVec ub_S1_C1_IM("ub_S1_C1_IM");
FVec ub_S1_C2_RE("ub_S1_C2_RE");
FVec ub_S1_C2_IM("ub_S1_C2_IM");

FVec ub_spinor[2][3][2] = {
    { {ub_S0_C0_RE, ub_S0_C0_IM}, {ub_S0_C1_RE, ub_S0_C1_IM}, { ub_S0_C2_RE, ub_S0_C2_IM } },
    { {ub_S1_C0_RE, ub_S1_C0_IM}, {ub_S1_C1_RE, ub_S1_C1_IM}, { ub_S1_C2_RE, ub_S1_C2_IM } }
};

FVec out_S0_C0_RE("out_S0_C0_RE");
FVec out_S0_C0_IM("out_S0_C0_IM");
FVec out_S0_C1_RE("out_S0_C1_RE");
FVec out_S0_C1_IM("out_S0_C1_IM");
FVec out_S0_C2_RE("out_S0_C2_RE");
FVec out_S0_C2_IM("out_S0_C2_IM");
FVec out_S1_C0_RE("out_S1_C0_RE");
FVec out_S1_C0_IM("out_S1_C0_IM");
FVec out_S1_C1_RE("out_S1_C1_RE");
FVec out_S1_C1_IM("out_S1_C1_IM");
FVec out_S1_C2_RE("out_S1_C2_RE");
FVec out_S1_C2_IM("out_S1_C2_IM");
FVec out_S2_C0_RE("out_S2_C0_RE");
FVec out_S2_C0_IM("out_S2_C0_IM");
FVec out_S2_C1_RE("out_S2_C1_RE");
FVec out_S2_C1_IM("out_S2_C1_IM");
FVec out_S2_C2_RE("out_S2_C2_RE");
FVec out_S2_C2_IM("out_S2_C2_IM");
FVec out_S3_C0_RE("out_S3_C0_RE");
FVec out_S3_C0_IM("out_S3_C0_IM");
FVec out_S3_C1_RE("out_S3_C1_RE");
FVec out_S3_C1_IM("out_S3_C1_IM");
FVec out_S3_C2_RE("out_S3_C2_RE");
FVec out_S3_C2_IM("out_S3_C2_IM");

FVec out_spinor[4][3][2] = {
    { {out_S0_C0_RE, out_S0_C0_IM}, {out_S0_C1_RE, out_S0_C1_IM}, { out_S0_C2_RE, out_S0_C2_IM } },
    { {out_S1_C0_RE, out_S1_C0_IM}, {out_S1_C1_RE, out_S1_C1_IM}, { out_S1_C2_RE, out_S1_C2_IM } },
    { {out_S2_C0_RE, out_S2_C0_IM}, {out_S2_C1_RE, out_S2_C1_IM}, { out_S2_C2_RE, out_S2_C2_IM } },
    { {out_S3_C0_RE, out_S3_C0_IM}, {out_S3_C1_RE, out_S3_C1_IM}, { out_S3_C2_RE, out_S3_C2_IM } }
};

FVec clout_spinor[2][6][2] = {
    {   {out_S0_C0_RE, out_S0_C0_IM}, {out_S0_C1_RE, out_S0_C1_IM}, { out_S0_C2_RE, out_S0_C2_IM },
        {out_S1_C0_RE, out_S1_C0_IM}, {out_S1_C1_RE, out_S1_C1_IM}, { out_S1_C2_RE, out_S1_C2_IM }
    },
    {   {out_S2_C0_RE, out_S2_C0_IM}, {out_S2_C1_RE, out_S2_C1_IM}, { out_S2_C2_RE, out_S2_C2_IM },
        {out_S3_C0_RE, out_S3_C0_IM}, {out_S3_C1_RE, out_S3_C1_IM}, { out_S3_C2_RE, out_S3_C2_IM }
    }
};


FVec chi_S0_C0_RE("chi_S0_C0_RE");
FVec chi_S0_C0_IM("chi_S0_C0_IM");
FVec chi_S0_C1_RE("chi_S0_C1_RE");
FVec chi_S0_C1_IM("chi_S0_C1_IM");
FVec chi_S0_C2_RE("chi_S0_C2_RE");
FVec chi_S0_C2_IM("chi_S0_C2_IM");
FVec chi_S1_C0_RE("chi_S1_C0_RE");
FVec chi_S1_C0_IM("chi_S1_C0_IM");
FVec chi_S1_C1_RE("chi_S1_C1_RE");
FVec chi_S1_C1_IM("chi_S1_C1_IM");
FVec chi_S1_C2_RE("chi_S1_C2_RE");
FVec chi_S1_C2_IM("chi_S1_C2_IM");
FVec chi_S2_C0_RE("chi_S2_C0_RE");
FVec chi_S2_C0_IM("chi_S2_C0_IM");
FVec chi_S2_C1_RE("chi_S2_C1_RE");
FVec chi_S2_C1_IM("chi_S2_C1_IM");
FVec chi_S2_C2_RE("chi_S2_C2_RE");
FVec chi_S2_C2_IM("chi_S2_C2_IM");
FVec chi_S3_C0_RE("chi_S3_C0_RE");
FVec chi_S3_C0_IM("chi_S3_C0_IM");
FVec chi_S3_C1_RE("chi_S3_C1_RE");
FVec chi_S3_C1_IM("chi_S3_C1_IM");
FVec chi_S3_C2_RE("chi_S3_C2_RE");
FVec chi_S3_C2_IM("chi_S3_C2_IM");

FVec chi_spinor[4][3][2] = {
    { {chi_S0_C0_RE, chi_S0_C0_IM}, {chi_S0_C1_RE, chi_S0_C1_IM}, { chi_S0_C2_RE, chi_S0_C2_IM } },
    { {chi_S1_C0_RE, chi_S1_C0_IM}, {chi_S1_C1_RE, chi_S1_C1_IM}, { chi_S1_C2_RE, chi_S1_C2_IM } },
    { {chi_S2_C0_RE, chi_S2_C0_IM}, {chi_S2_C1_RE, chi_S2_C1_IM}, { chi_S2_C2_RE, chi_S2_C2_IM } },
    { {chi_S3_C0_RE, chi_S3_C0_IM}, {chi_S3_C1_RE, chi_S3_C1_IM}, { chi_S3_C2_RE, chi_S3_C2_IM } }
};

FVec dout_S0_C0_RE("dout_S0_C0_RE");
FVec dout_S0_C0_IM("dout_S0_C0_IM");
FVec dout_S0_C1_RE("dout_S0_C1_RE");
FVec dout_S0_C1_IM("dout_S0_C1_IM");
FVec dout_S0_C2_RE("dout_S0_C2_RE");
FVec dout_S0_C2_IM("dout_S0_C2_IM");
FVec dout_S1_C0_RE("dout_S1_C0_RE");
FVec dout_S1_C0_IM("dout_S1_C0_IM");
FVec dout_S1_C1_RE("dout_S1_C1_RE");
FVec dout_S1_C1_IM("dout_S1_C1_IM");
FVec dout_S1_C2_RE("dout_S1_C2_RE");
FVec dout_S1_C2_IM("dout_S1_C2_IM");
FVec dout_S2_C0_RE("dout_S2_C0_RE");
FVec dout_S2_C0_IM("dout_S2_C0_IM");
FVec dout_S2_C1_RE("dout_S2_C1_RE");
FVec dout_S2_C1_IM("dout_S2_C1_IM");
FVec dout_S2_C2_RE("dout_S2_C2_RE");
FVec dout_S2_C2_IM("dout_S2_C2_IM");
FVec dout_S3_C0_RE("dout_S3_C0_RE");
FVec dout_S3_C0_IM("dout_S3_C0_IM");
FVec dout_S3_C1_RE("dout_S3_C1_RE");
FVec dout_S3_C1_IM("dout_S3_C1_IM");
FVec dout_S3_C2_RE("dout_S3_C2_RE");
FVec dout_S3_C2_IM("dout_S3_C2_IM");

FVec dout_spinor[4][3][2] = {
    { {dout_S0_C0_RE, dout_S0_C0_IM}, {dout_S0_C1_RE, dout_S0_C1_IM}, { dout_S0_C2_RE, dout_S0_C2_IM } },
    { {dout_S1_C0_RE, dout_S1_C0_IM}, {dout_S1_C1_RE, dout_S1_C1_IM}, { dout_S1_C2_RE, dout_S1_C2_IM } },
    { {dout_S2_C0_RE, dout_S2_C0_IM}, {dout_S2_C1_RE, dout_S2_C1_IM}, { dout_S2_C2_RE, dout_S2_C2_IM } },
    { {dout_S3_C0_RE, dout_S3_C0_IM}, {dout_S3_C1_RE, dout_S3_C1_IM}, { dout_S3_C2_RE, dout_S3_C2_IM } }
};


FVec cl_diag_0("cl_diag_0");
FVec cl_diag_1("cl_diag_1");
FVec cl_diag_2("cl_diag_2");
FVec cl_diag_3("cl_diag_3");
FVec cl_diag_4("cl_diag_4");
FVec cl_diag_5("cl_diag_5");

FVec clov_diag[6] = { cl_diag_0, cl_diag_1, cl_diag_2, cl_diag_3, cl_diag_4, cl_diag_5 };

FVec cl_offdiag_0_RE("cl_offdiag_0_RE");
FVec cl_offdiag_0_IM("cl_offdiag_0_IM");
FVec cl_offdiag_1_RE("cl_offdiag_1_RE");
FVec cl_offdiag_1_IM("cl_offdiag_1_IM");
FVec cl_offdiag_2_RE("cl_offdiag_2_RE");
FVec cl_offdiag_2_IM("cl_offdiag_2_IM");
FVec cl_offdiag_3_RE("cl_offdiag_3_RE");
FVec cl_offdiag_3_IM("cl_offdiag_3_IM");
FVec cl_offdiag_4_RE("cl_offdiag_4_RE");
FVec cl_offdiag_4_IM("cl_offdiag_4_IM");
FVec cl_offdiag_5_RE("cl_offdiag_5_RE");
FVec cl_offdiag_5_IM("cl_offdiag_5_IM");
FVec cl_offdiag_6_RE("cl_offdiag_6_RE");
FVec cl_offdiag_6_IM("cl_offdiag_6_IM");
FVec cl_offdiag_7_RE("cl_offdiag_7_RE");
FVec cl_offdiag_7_IM("cl_offdiag_7_IM");
FVec cl_offdiag_8_RE("cl_offdiag_8_RE");
FVec cl_offdiag_8_IM("cl_offdiag_8_IM");
FVec cl_offdiag_9_RE("cl_offdiag_9_RE");
FVec cl_offdiag_9_IM("cl_offdiag_9_IM");
FVec cl_offdiag_10_RE("cl_offdiag_10_RE");
FVec cl_offdiag_10_IM("cl_offdiag_10_IM");
FVec cl_offdiag_11_RE("cl_offdiag_11_RE");
FVec cl_offdiag_11_IM("cl_offdiag_11_IM");
FVec cl_offdiag_12_RE("cl_offdiag_12_RE");
FVec cl_offdiag_12_IM("cl_offdiag_12_IM");
FVec cl_offdiag_13_RE("cl_offdiag_13_RE");
FVec cl_offdiag_13_IM("cl_offdiag_13_IM");
FVec cl_offdiag_14_RE("cl_offdiag_14_RE");
FVec cl_offdiag_14_IM("cl_offdiag_14_IM");

FVec clov_offdiag[15][2]= {
    { cl_offdiag_0_RE, cl_offdiag_0_IM },
    { cl_offdiag_1_RE, cl_offdiag_1_IM },
    { cl_offdiag_2_RE, cl_offdiag_2_IM },
    { cl_offdiag_3_RE, cl_offdiag_3_IM },
    { cl_offdiag_4_RE, cl_offdiag_4_IM },
    { cl_offdiag_5_RE, cl_offdiag_5_IM },
    { cl_offdiag_6_RE, cl_offdiag_6_IM },
    { cl_offdiag_7_RE, cl_offdiag_7_IM },
    { cl_offdiag_8_RE, cl_offdiag_8_IM },
    { cl_offdiag_9_RE, cl_offdiag_9_IM },
    { cl_offdiag_10_RE, cl_offdiag_10_IM },
    { cl_offdiag_11_RE, cl_offdiag_11_IM },
    { cl_offdiag_12_RE, cl_offdiag_12_IM },
    { cl_offdiag_13_RE, cl_offdiag_13_IM },
    { cl_offdiag_14_RE, cl_offdiag_14_IM }
};

FVec zero("zero");
FVec alpha_vec("alpha_vec");
FVec beta_vec("beta_vec");

FVec psi_S0_RE("psi_S0_RE");
FVec psi_S0_IM("psi_S0_IM");
FVec psi_S1_RE("psi_S1_RE");
FVec psi_S1_IM("psi_S1_IM");

FVec psi[2][2] = { {psi_S0_RE, psi_S0_IM}, {psi_S1_RE, psi_S1_IM} };

FVec tmp_1_re("tmp_1_re");
FVec tmp_1_im("tmp_1_im");

FVec tmp_2_re("tmp_2_re");
FVec tmp_2_im("tmp_2_im");

FVec tmp_3_re("tmp_3_re");
FVec tmp_3_im("tmp_3_im");

FVec tmp_4_re("tmp_4_re");
FVec tmp_4_im("tmp_4_im");

FVec tmp[4] = {tmp_1_re, tmp_2_re, tmp_3_re, tmp_4_re};

FVec u_00_re("u_00_re");
FVec u_00_im("u_00_im");
FVec u_10_re("u_10_re");
FVec u_10_im("u_10_im");
FVec u_20_re("u_20_re");
FVec u_20_im("u_20_im");

FVec u_01_re("u_01_re");
FVec u_01_im("u_01_im");
FVec u_11_re("u_11_re");
FVec u_11_im("u_11_im");
FVec u_21_re("u_21_re");
FVec u_21_im("u_21_im");

FVec u_02_re("u_02_re");
FVec u_02_im("u_02_im");
FVec u_12_re("u_12_re");
FVec u_12_im("u_12_im");
FVec u_22_re("u_22_re");
FVec u_22_im("u_22_im");

FVec u_gauge[3][3][2] = {
    { {u_00_re, u_00_im}, {u_01_re, u_01_im}, {u_02_re, u_02_im} },
    { {u_10_re, u_10_im}, {u_11_re, u_11_im}, {u_12_re, u_12_im} },
    { {u_20_re, u_20_im}, {u_21_re, u_21_im}, {u_22_re, u_22_im} }
};

void declare_b_Spins(InstVector& ivector)
{
    for(int s=0; s < 2; s++) {
        for(int c = 0; c < 3; c++) {
            declareFVecFromFVec(ivector, b_spinor[s][c][RE]);
            declareFVecFromFVec(ivector, b_spinor[s][c][IM]);
        }
    }
}

void declare_ub_Spins(InstVector& ivector)
{
    for(int s=0; s < 2; s++) {
        for(int c = 0; c < 3; c++) {
            declareFVecFromFVec(ivector, ub_spinor[s][c][RE]);
            declareFVecFromFVec(ivector, ub_spinor[s][c][IM]);
        }
    }
}

void declare_u_gaus(InstVector& ivector)
{
    for(int c1=0; c1 < 3; c1++) {
        for(int c2 = 0; c2 < 3; c2++) {
            declareFVecFromFVec(ivector, u_gauge[c1][c2][RE]);
            declareFVecFromFVec(ivector, u_gauge[c1][c2][IM]);
        }
    }
}

void declare_outs(InstVector& ivector)
{
    for(int s=0; s < 4; s++) {
        for(int c = 0; c < 3; c++) {
            declareFVecFromFVec(ivector, out_spinor[s][c][RE]);
            declareFVecFromFVec(ivector, out_spinor[s][c][IM]);
        }
    }
}

void declare_douts(InstVector& ivector)
{
    for(int s=0; s < 4; s++) {
        for(int c = 0; c < 3; c++) {
            declareFVecFromFVec(ivector, dout_spinor[s][c][RE]);
            declareFVecFromFVec(ivector, dout_spinor[s][c][IM]);
        }
    }
}

void declare_chi(InstVector& ivector)
{
    for(int s=0; s < 4; s++) {
        for(int c = 0; c < 3; c++) {
            declareFVecFromFVec(ivector, chi_spinor[s][c][RE]);
            declareFVecFromFVec(ivector, chi_spinor[s][c][IM]);
        }
    }
}

void declare_clover(InstVector& ivector)
{
    for(int s=0; s < 6; s++) {
        declareFVecFromFVec(ivector, clov_diag[s]);
    }

    for(int s=0; s < 15; s++) {
        declareFVecFromFVec(ivector, clov_offdiag[s][RE]);
        declareFVecFromFVec(ivector, clov_offdiag[s][IM]);
    }
}

void declare_misc(InstVector& ivector)
{
    declareFVecFromFVec(ivector, psi_S0_RE);
    declareFVecFromFVec(ivector, psi_S0_IM);
    declareFVecFromFVec(ivector, psi_S1_RE);
    declareFVecFromFVec(ivector, psi_S1_IM);

    declareFVecFromFVec(ivector, tmp_1_re);
    declareFVecFromFVec(ivector, tmp_1_im);
    declareFVecFromFVec(ivector, tmp_2_re);
    declareFVecFromFVec(ivector, tmp_2_im);
    declareFVecFromFVec(ivector, tmp_3_re);
    declareFVecFromFVec(ivector, tmp_3_im);
    declareFVecFromFVec(ivector, tmp_4_re);
    declareFVecFromFVec(ivector, tmp_4_im);

    declareFVecFromFVec(ivector, zero);
    setZero(ivector, zero);
}

void movCVec(InstVector& ivector, FVec *r, FVec *s1, string &mask)
{
    movFVec(ivector, r[RE], s1[RE], mask);
    movFVec(ivector, r[IM], s1[IM], mask);
}

void addCVec(InstVector& ivector, FVec *r, FVec *s1, FVec *s2, string &mask)
{
    addFVec(ivector, r[RE], s1[RE], s2[RE], mask);
    addFVec(ivector, r[IM], s1[IM], s2[IM], mask);
}

void subCVec(InstVector& ivector, FVec *r, FVec *s1, FVec *s2, string &mask)
{
    subFVec(ivector, r[RE], s1[RE], s2[RE], mask);
    subFVec(ivector, r[IM], s1[IM], s2[IM], mask);
}

void addiCVec(InstVector& ivector, FVec *r, FVec *s1, FVec *s2, string &mask)
{
    subFVec(ivector, r[RE], s1[RE], s2[IM], mask);
    addFVec(ivector, r[IM], s1[IM], s2[RE], mask);
}

void subiCVec(InstVector& ivector, FVec *r, FVec *s1, FVec *s2, string &mask)
{
    addFVec(ivector, r[RE], s1[RE], s2[IM], mask);
    subFVec(ivector, r[IM], s1[IM], s2[RE], mask);
}

// r[RE] = s1[RE]-beta_vec*s2[RE] = fnmadd(beta_vec,s2[RE],s1[RE])
// r[IM] = s1[IM]-beta_vec*s2[IM] = fnamdd(beta_vec,s2[IM],s1[IM])
void addCVec_mbeta(InstVector& ivector, FVec *r, FVec *s1, FVec *s2, string &mask)
{
    fnmaddFVec(ivector, r[RE], beta_vec, s2[RE], s1[RE], mask);
    fnmaddFVec(ivector, r[IM], beta_vec, s2[IM], s1[IM], mask);
}
// r[RE] = s1[RE] + beta_vec*s2[RE] = fmadd(beta_vec, s2[RE], s1[RE]);
// r[IM] = s1[IM] + beta_vec*s2[IM] = fmadd(beta_vec, s2[IM], s1[IM]);
void subCVec_mbeta(InstVector& ivector, FVec *r, FVec *s1, FVec *s2, string &mask)
{
    fmaddFVec(ivector, r[RE], beta_vec, s2[RE], s1[RE], mask);
    fmaddFVec(ivector, r[IM], beta_vec, s2[IM], s1[IM], mask);
}

// r[RE] = s1[RE] + beta_vec * s2[IM] = fmadd(beta_vec,s2[IM], s1[RE])
// r[IM] = s1[IM] - beta_vec * s2[RE] = fnmadd(beta_vec, s2[RE], s1[IM])

void addiCVec_mbeta(InstVector& ivector, FVec *r, FVec *s1, FVec *s2, string &mask)
{
    fmaddFVec(ivector, r[RE], beta_vec, s2[IM], s1[RE], mask);
    fnmaddFVec(ivector, r[IM], beta_vec, s2[RE], s1[IM], mask);
}

// r[RE] = s1[RE] - beta_vec*s2[IM] = fnmadd( beta_vec, s2[IM], s1[RE]);
// r[IM] = s1[IM] + beta_vec*s2[RE] = fmadd ( beta_vec, s2[RE], s1[IM]);
void subiCVec_mbeta(InstVector& ivector, FVec *r, FVec *s1, FVec *s2, string &mask)
{
    fnmaddFVec(ivector, r[RE], beta_vec, s2[IM], s1[RE], mask);
    fmaddFVec(ivector, r[IM], beta_vec, s2[RE], s1[IM], mask);
}

// r[RE] = s1[RE]+beta_vec*s2[RE] = fmadd(beta_vec,s2[RE],s1[RE])
// r[IM] = s1[IM]+beta_vec*s2[IM] = fmadd(beta_vec,s2[IM],s1[IM])
void addCVec_pbeta(InstVector& ivector, FVec *r, FVec *s1, FVec *s2, string &mask)
{
    fmaddFVec(ivector, r[RE], beta_vec, s2[RE], s1[RE], mask);
    fmaddFVec(ivector, r[IM], beta_vec, s2[IM], s1[IM], mask);
}
// r[RE] = s1[RE] - beta_vec*s2[RE] = fnmadd(beta_vec, s2[RE], s1[RE]);
// r[IM] = s1[IM] - beta_vec*s2[IM] = fnmadd(beta_vec, s2[IM], s1[IM]);
void subCVec_pbeta(InstVector& ivector, FVec *r, FVec *s1, FVec *s2, string &mask)
{
    fnmaddFVec(ivector, r[RE], beta_vec, s2[RE], s1[RE], mask);
    fnmaddFVec(ivector, r[IM], beta_vec, s2[IM], s1[IM], mask);
}

// r[RE] = s1[RE] - beta_vec * s2[IM] = fnmadd(beta_vec,s2[IM], s1[RE])
// r[IM] = s1[IM] + beta_vec * s2[RE] = fmadd(beta_vec, s2[RE], s1[IM])
void addiCVec_pbeta(InstVector& ivector, FVec *r, FVec *s1, FVec *s2, string &mask)
{
    fnmaddFVec(ivector, r[RE], beta_vec, s2[IM], s1[RE], mask);
    fmaddFVec(ivector, r[IM], beta_vec, s2[RE], s1[IM], mask);
}

// r[RE] = s1[RE] + beta_vec*s2[IM] = fmadd( beta_vec, s2[IM], s1[RE]);
// r[IM] = s1[IM] - beta_vec*s2[RE] = fnmadd ( beta_vec, s2[RE], s1[IM]);
void subiCVec_pbeta(InstVector& ivector, FVec *r, FVec *s1, FVec *s2, string &mask)
{
    fmaddFVec(ivector, r[RE], beta_vec, s2[IM], s1[RE], mask);
    fnmaddFVec(ivector, r[IM], beta_vec, s2[RE], s1[IM], mask);
}

// r = s1*s2
void mulCVec(InstVector& ivector, FVec *r, FVec *s1, FVec *s2, string &mask)
{
    mulFVec(ivector, r[RE], s1[RE], s2[RE], mask);
    fnmaddFVec(ivector, r[RE], s1[IM], s2[IM], r[RE], mask);
    mulFVec(ivector, r[IM], s1[RE], s2[IM], mask);
    fmaddFVec(ivector, r[IM], s1[IM], s2[RE], r[IM], mask);
}

// r = s1*s2+s3
void fmaddCVec(InstVector& ivector, FVec *r, FVec *s1, FVec *s2,  FVec *s3, string &mask)
{
    fmaddFVec(ivector, r[RE], s1[RE], s2[RE], s3[RE], mask);
    fnmaddFVec(ivector, r[RE], s1[IM], s2[IM], r[RE], mask);
    fmaddFVec(ivector, r[IM], s1[RE], s2[IM], s3[IM], mask);
    fmaddFVec(ivector, r[IM], s1[IM], s2[RE], r[IM], mask);
}

// r = s3-s1*s2
//r[RE] = (s3[RE]-s1[RE]*s2[RE])+(s1[IM]*s2[IM])
//r[IM] = (s3[IM]-s1[RE]*s2[IM])-(s1[IM]*s2[RE])
void fnmaddCVec(InstVector& ivector, FVec *r, FVec *s1, FVec *s2,  FVec *s3, string &mask)
{
    fnmaddFVec(ivector, r[RE], s1[RE], s2[RE], s3[RE], mask);
    fmaddFVec(ivector, r[RE], s1[IM], s2[IM], r[RE], mask);
    fnmaddFVec(ivector, r[IM], s1[RE], s2[IM], s3[IM], mask);
    fnmaddFVec(ivector, r[IM], s1[IM], s2[RE], r[IM], mask);
}

// r = (s1*s2-s3*s4)'
//r[RE] = (s1[RE]*s2[RE])-(s1[IM]*s2[IM])-(s3[RE]*s4[RE])+(s3[IM]*s4[IM])
//r[IM] = (s3[RE]*s4[IM])+(s3[IM]*s4[RE])-(s1[RE]*s2[IM])-(s1[IM]*s2[RE])
void Conj_CrossProd(InstVector& ivector, FVec *r, FVec *s1, FVec *s2, FVec *s3, FVec *s4, string &mask)
{
    mulFVec(ivector, r[RE], s1[RE], s2[RE], mask);
    fnmaddFVec(ivector, r[RE], s1[IM], s2[IM], r[RE], mask);
    fnmaddFVec(ivector, r[RE], s3[RE], s4[RE], r[RE], mask);
    fmaddFVec(ivector, r[RE], s3[IM], s4[IM], r[RE], mask);

    mulFVec(ivector, r[IM], s3[RE], s4[IM], mask);
    fmaddFVec(ivector, r[IM], s3[IM], s4[RE], r[IM], mask);
    fnmaddFVec(ivector, r[IM], s1[RE], s2[IM], r[IM], mask);
    fnmaddFVec(ivector, r[IM], s1[IM], s2[RE], r[IM], mask);
}

// r = s1'*s2
//r[RE] = (s1[RE]*s2[RE])+(s1[IM]*s2[IM])
//r[IM] = (s1[RE]*s2[IM])-(s1[IM]*s2[RE])
void mulConjCVec(InstVector& ivector, FVec *r, FVec *s1, FVec *s2, string &mask)
{
    mulFVec(ivector, r[RE], s1[RE], s2[RE], mask);
    fmaddFVec(ivector, r[RE], s1[IM], s2[IM], r[RE], mask);
    mulFVec(ivector, r[IM], s1[RE], s2[IM], mask);
    fnmaddFVec(ivector, r[IM], s1[IM], s2[RE], r[IM], mask);
}

// r = s1'*s2+s3
//r[RE] = (s3[RE]+s1[RE]*s2[RE])+(s1[IM]*s2[IM])
//r[IM] = (s3[IM]+s1[RE]*s2[IM])-(s1[IM]*s2[RE])
void fmaddConjCVec(InstVector& ivector, FVec *r, FVec *s1, FVec *s2,  FVec *s3, string &mask)
{
    fmaddFVec(ivector, r[RE], s1[RE], s2[RE], s3[RE], mask);
    fmaddFVec(ivector, r[RE], s1[IM], s2[IM], r[RE], mask);
    fmaddFVec(ivector, r[IM], s1[RE], s2[IM], s3[IM], mask);
    fnmaddFVec(ivector, r[IM], s1[IM], s2[RE], r[IM], mask);
}

proj_ops proj_ops_plus[] = {
    {"plus_X_back", {{0,3},{1,2}},{addiCVec,addiCVec}},
    {"minus_X_forw",{{0,3},{1,2}},{subiCVec,subiCVec}},
    {"plus_Y",      {{0,3},{1,2}},{subCVec, addCVec}},
    {"minus_Y",     {{0,3},{1,2}},{addCVec, subCVec}},
    {"plus_Z",      {{0,2},{1,3}},{addiCVec,subiCVec}},
    {"minus_Z",     {{0,2},{1,3}},{subiCVec,addiCVec}},
    {"plus_T",      {{0,2},{1,3}},{addCVec, addCVec}},
    {"minus_T",     {{0,2},{1,3}},{subCVec, subCVec}}
};

proj_ops proj_ops_minus[] = {
    {"minus_X_back",{{0,3},{1,2}},{subiCVec,subiCVec}},
    {"plus_X_forw", {{0,3},{1,2}},{addiCVec,addiCVec}},
    {"minus_Y",     {{0,3},{1,2}},{addCVec, subCVec}},
    {"plus_Y",      {{0,3},{1,2}},{subCVec, addCVec}},
    {"minus_Z",     {{0,2},{1,3}},{subiCVec,addiCVec}},
    {"plus_Z",      {{0,2},{1,3}},{addiCVec,subiCVec}},
    {"minus_T",     {{0,2},{1,3}},{subCVec, subCVec}},
    {"plus_T",      {{0,2},{1,3}},{addCVec, addCVec}}
};

/*
recons_ops rec_plus_ops[] = {
	{"recons_plus_X", 1,0, addCVec, subiCVec, subiCVec},
	{"recons_plus_Y", 1,0, addCVec, addCVec,  subCVec},
	{"recons_plus_Z", 0,1, addCVec, subiCVec, addiCVec},
	{"recons_plus_T", 0,1, addCVec, addCVec,  addCVec},
};

recons_ops rec_minus_ops[] = {
	{"recons_minus_X", 1,0, addCVec, addiCVec, addiCVec},
	{"recons_minus_Y", 1,0, addCVec, subCVec,  addCVec},
	{"recons_minus_Z", 0,1, addCVec, addiCVec, subiCVec},
	{"recons_minus_T", 0,1, addCVec, subCVec,  subCVec},
};
*/

recons_ops rec_plus_pbeta_ops[] = {
    {"recons_plus_X_pbeta", 1,0, addCVec_pbeta, subiCVec_pbeta, subiCVec_pbeta},
    {"recons_plus_Y_pbeta", 1,0, addCVec_pbeta, addCVec_pbeta,  subCVec_pbeta},
    {"recons_plus_Z_pbeta", 0,1, addCVec_pbeta, subiCVec_pbeta, addiCVec_pbeta},
    {"recons_plus_T_pbeta", 0,1, addCVec_pbeta, addCVec_pbeta, addCVec_pbeta},
};

recons_ops rec_minus_pbeta_ops[] = {
    {"recons_minus_X", 1,0, addCVec_pbeta, addiCVec_pbeta, addiCVec_pbeta},
    {"recons_minus_Y", 1,0, addCVec_pbeta, subCVec_pbeta,  addCVec_pbeta},
    {"recons_minus_Z", 0,1, addCVec_pbeta, addiCVec_pbeta, subiCVec_pbeta},
    {"recons_minus_T", 0,1, addCVec_pbeta, subCVec_pbeta,  subCVec_pbeta},
};

recons_ops rec_plus_mbeta_ops[] = {
    {"recons_plus_X_mbeta", 1,0, addCVec_mbeta, subiCVec_mbeta, subiCVec_mbeta},
    {"recons_plus_Y_mbeta", 1,0, addCVec_mbeta, addCVec_mbeta,  subCVec_mbeta},
    {"recons_plus_Z_mbeta", 0,1, addCVec_mbeta, subiCVec_mbeta, addiCVec_mbeta},
    {"recons_plus_T_mbeta", 0,1, addCVec_mbeta, addCVec_mbeta, addCVec_mbeta},
};

recons_ops rec_minus_mbeta_ops[] = {
    {"recons_minus_X", 1,0, addCVec_mbeta, addiCVec_mbeta, addiCVec_mbeta},
    {"recons_minus_Y", 1,0, addCVec_mbeta, subCVec_mbeta,  addCVec_mbeta},
    {"recons_minus_Z", 0,1, addCVec_mbeta, addiCVec_mbeta, subiCVec_mbeta},
    {"recons_minus_T", 0,1, addCVec_mbeta, subCVec_mbeta,  subCVec_mbeta},
};


// Uses ub_spinor as implicit input
// Uses psi[][] as temp and b_spinor[][] as implicit output
void project(InstVector& ivector, string base, string offset, proj_ops& ops, bool isFace, string mask, int dir)
{
    string tmask("");
    PrefetchL1FullSpinorDirIn(ivector, base, offset, dir);

    for(int s = 0; s < 2; s++) {
        for(int c = 0; c < 3; c++) {
            LoadSpinorElement(ivector, psi[0][RE], base, offset, ops.s[s][0], c, RE, isFace,  mask, dir);
            LoadSpinorElement(ivector, psi[0][IM], base, offset, ops.s[s][0], c, IM, isFace,  mask, dir);
            LoadSpinorElement(ivector, psi[1][RE], base, offset, ops.s[s][1], c, RE, isFace,  mask, dir);
            LoadSpinorElement(ivector, psi[1][IM], base, offset, ops.s[s][1], c, IM, isFace,  mask, dir);

            ops.CVecFunc[s](ivector, b_spinor[s][c], psi[0], psi[1], /*mask*/ tmask); // Not using mask here
        }
    }
}

// Serial Spin version
void project(InstVector& ivector, string base, string offset, proj_ops& ops, bool isFace, string mask, int dir, int s)
{
    string tmask("");

    if(s==0) {
        PrefetchL1FullSpinorDirIn(ivector, base, offset, dir);
    }

    for(int c = 0; c < 3; c++) {
        LoadSpinorElement(ivector, psi[0][RE], base, offset, ops.s[s][0], c, RE, isFace,  mask, dir);
        LoadSpinorElement(ivector, psi[0][IM], base, offset, ops.s[s][0], c, IM, isFace,  mask, dir);
        LoadSpinorElement(ivector, psi[1][RE], base, offset, ops.s[s][1], c, RE, isFace,  mask, dir);
        LoadSpinorElement(ivector, psi[1][IM], base, offset, ops.s[s][1], c, IM, isFace,  mask, dir);

        ops.CVecFunc[s](ivector, b_spinor[s][c], psi[0], psi[1], /*mask*/ tmask); // Not using mask here
    }
}

void recons_add(InstVector& ivector, recons_ops& ops, FVec outspinor[4][3][2], string &mask)
{
    for(int s=0; s < 2; s++) {
        for(int c = 0; c < 3; c++) {
            ops.CVecFuncTop2(ivector, outspinor[s][c], outspinor[s][c], ub_spinor[s][c], mask);
        }

        if(ops.s2 == s) {
            for(int c=0; c < 3; c++) {
                ops.CVecFunc1(ivector, outspinor[2][c], outspinor[2][c], ub_spinor[s][c], mask);
            }
        }
        else {
            for(int c=0; c < 3; c++) {
                ops.CVecFunc2(ivector, outspinor[3][c], outspinor[3][c], ub_spinor[s][c], mask);
            }
        }
    }
}

// Serial Spin version
void recons_add(InstVector& ivector, recons_ops& ops, FVec outspinor[4][3][2], string &mask, int s)
{
    for(int c = 0; c < 3; c++) {
        ops.CVecFuncTop2(ivector, outspinor[s][c], outspinor[s][c], ub_spinor[s][c], mask);
    }

    if(ops.s2 == s) {
        for(int c=0; c < 3; c++) {
            ops.CVecFunc1(ivector, outspinor[2][c], outspinor[2][c], ub_spinor[s][c], mask);
        }
    }
    else {
        for(int c=0; c < 3; c++) {
            ops.CVecFunc2(ivector, outspinor[3][c], outspinor[3][c], ub_spinor[s][c], mask);
        }
    }
}

void zeroResult(InstVector& ivector, FVec *outspinor)
{
    for(int i=0; i < 24; i++) {
        setZero(ivector,outspinor[i]);
    }
}

void clover_term(InstVector& ivector, FVec in_spinor[4][3][2], bool face, string _mask)
{
    FVec clout_tmp[2] = {tmp_1_re, tmp_1_im};

    for(int block=0; block < 2; block++) {
        PrefetchL1FullCloverBlockIn(ivector, clBase, clOffs, block);
        LoadFullCloverBlock(ivector, clov_diag, clov_offdiag, clBase, clOffs, block);

        for(int c1=0; c1 < 6; c1++) {
            int spin = 2*block+c1/3;
            int col = c1 % 3;
            bool acc = face;
            string mask = _mask;
            FVec *clout = out_spinor[spin][col];
            FVec *clin  = in_spinor[spin][col];
#ifdef NO_HW_MASKING

            if(_mask != "") {
                acc = false;
                clout = clout_tmp;
                mask = "";
            }

#endif

            if( acc ) {
                fmaddFVec( ivector, clout[RE], clov_diag[c1], clin[RE], clout[RE], mask);
                fmaddFVec( ivector, clout[IM], clov_diag[c1], clin[IM], clout[IM], mask);
            }
            else {
                mulFVec( ivector,  clout[RE], clov_diag[c1],  clin[RE], mask);
                mulFVec( ivector,  clout[IM], clov_diag[c1],  clin[IM], mask);
            }

            for(int c2=0; c2 < 6; c2++) {
                if(c1 == c2) {
                    continue;    // diagonal case
                }

                if(c1 < c2) {
                    int od = c2*(c2-1)/2+c1;
                    fmaddConjCVec(ivector, clout, clov_offdiag[od], in_spinor[2*block+c2/3][c2%3], clout, mask);
                }
                else {
                    int od = c1*(c1-1)/2+c2;
                    fmaddCVec(ivector, clout, clov_offdiag[od], in_spinor[2*block+c2/3][c2%3], clout, mask);
                }
            }

#ifdef NO_HW_MASKING

            if(_mask != "") {
                if(face) {
                    addCVec(ivector, out_spinor[spin][col], clout, clout_spinor[block][c1], _mask);
                }
                else {
                    movCVec(ivector, out_spinor[spin][col], clout, _mask);
                }
            }

#endif
        }
    }
}

void achiResult(InstVector& ivector, bool clover)
{
    PrefetchL1FullSpinorDirIn(ivector, chiBase, chiOffs, -1, 1 /*NTA*/);

    if(!clover) {
        for(int col=0; col < 3; col++) {
            for(int spin=0; spin < 4; spin++) {
                LoadSpinorElement(ivector, tmp_1_re, chiBase, chiOffs, spin, col, RE, false, "");
                LoadSpinorElement(ivector, tmp_1_im, chiBase, chiOffs, spin, col, IM, false, "");
                mulFVec(ivector, out_spinor[spin][col][RE], alpha_vec, tmp_1_re);
                mulFVec(ivector, out_spinor[spin][col][IM], alpha_vec, tmp_1_im);
            }
        }
    }
    else {
        for(int col=0; col < 3; col++) {
            for(int spin=0; spin < 4; spin++) {
                LoadSpinorElement(ivector, chi_spinor[spin][col][RE], chiBase, chiOffs, spin, col, RE, false, "");
                LoadSpinorElement(ivector, chi_spinor[spin][col][IM], chiBase, chiOffs, spin, col, IM, false, "");
            }
        }

        // Apply clover term, and store result in out spinor.
        // This is only on the AChi - bDPsi op (achimbdpsi = true)
        // This is only in body kernel (face = false)
        clover_term(ivector, chi_spinor, false);
    }
}

void loadGaugeDir(InstVector& ivector, int dir, bool compress12)
{
    string mask;

    PrefetchL1FullGaugeDirIn(ivector, gBase, gOffs, dir, compress12);
    LoadFullGaugeDir(ivector, u_gauge, gBase, gOffs, dir, compress12);

    if( compress12 ) {
        //printf("Using Compressed Gauges\n");
        for(int c = 0; c < 3; c++) {
            Conj_CrossProd(ivector, u_gauge[2][c], u_gauge[0][(c+1)%3], u_gauge[1][(c+2)%3], u_gauge[0][(c+2)%3], u_gauge[1][(c+1)%3], mask);
        }
    }
}

void matMultVec(InstVector& ivector, bool adjMul, int s)
{
    string mask;

    for(int c1 = 0; c1 < 3; c1++) {
        if(!adjMul) {
            mulCVec(ivector, ub_spinor[s][c1], u_gauge[0][c1], b_spinor[s][0], mask);
            fmaddCVec(ivector, ub_spinor[s][c1], u_gauge[1][c1], b_spinor[s][1], ub_spinor[s][c1], mask);
            fmaddCVec(ivector, ub_spinor[s][c1], u_gauge[2][c1], b_spinor[s][2], ub_spinor[s][c1], mask);
        }
        else {
            mulConjCVec(ivector, ub_spinor[s][c1], u_gauge[c1][0], b_spinor[s][0], mask);
            fmaddConjCVec(ivector, ub_spinor[s][c1], u_gauge[c1][1], b_spinor[s][1], ub_spinor[s][c1], mask);
            fmaddConjCVec(ivector, ub_spinor[s][c1], u_gauge[c1][2], b_spinor[s][2], ub_spinor[s][c1], mask);
        }
    }
}

void matMultVec(InstVector& ivector, bool adjMul)
{
    matMultVec(ivector, adjMul, 0);
    matMultVec(ivector, adjMul, 1);
}


void dslash_plain_body(InstVector& ivector, bool compress12, bool clover, bool isPlus)
{
    declare_b_Spins(ivector);
    declare_ub_Spins(ivector);
    declare_u_gaus(ivector);
    declare_misc(ivector);

    declare_outs(ivector);

    if(clover) {
        declare_douts(ivector);
        declare_clover(ivector);
    }

    FVec (*outspinor)[4][3][2];

    if(clover) {
        outspinor = &dout_spinor;
    }
    else {
        outspinor = &out_spinor;
    }

    zeroResult(ivector, (*outspinor)[0][0]);

    proj_ops *p_ops;
    recons_ops *rec_ops_bw;
    recons_ops *rec_ops_fw;

    if(isPlus) {
        p_ops 		= proj_ops_plus;
        rec_ops_bw	= rec_plus_pbeta_ops;
        rec_ops_fw 	= rec_minus_pbeta_ops;
    }
    else {
        p_ops 		= proj_ops_minus;
        rec_ops_bw	= rec_minus_pbeta_ops;
        rec_ops_fw 	= rec_plus_pbeta_ops;
    }

    dslash_body(ivector, compress12, p_ops, rec_ops_bw, rec_ops_fw, *outspinor);

    if(clover) {
        clover_term(ivector, *outspinor, false);
    }

    // Store
    StreamFullSpinor(ivector, out_spinor, outBase, outOffs);
}

// ***** ------- a chi - b D psi versions
void dslash_achimbdpsi_body(InstVector& ivector, bool compress12, bool clover, bool isPlus)
{
    declare_b_Spins(ivector);
    declare_ub_Spins(ivector);
    declare_u_gaus(ivector);
    declare_misc(ivector);

    declare_outs(ivector);
    declare_chi(ivector);

    if(clover) {
        declare_douts(ivector);
        declare_clover(ivector);
    }
    else {
        declareFVecFromFVec(ivector, alpha_vec);
        loadBroadcastScalar(ivector, alpha_vec, alpha_name, SpinorType);
    }

    // Fill result with a*chi
    achiResult(ivector, clover);

    proj_ops *p_ops;
    recons_ops *rec_ops_bw;
    recons_ops *rec_ops_fw;

    if(isPlus) {
        p_ops 		= proj_ops_plus;
        rec_ops_bw	= rec_plus_mbeta_ops;
        rec_ops_fw 	= rec_minus_mbeta_ops;
    }
    else {
        p_ops 		= proj_ops_minus;
        rec_ops_bw	= rec_minus_mbeta_ops;
        rec_ops_fw 	= rec_plus_mbeta_ops;
    }

    dslash_body(ivector, compress12, p_ops, rec_ops_bw, rec_ops_fw, out_spinor);

    // Store
    StreamFullSpinor(ivector, out_spinor, outBase, outOffs);
}

void pack_face_to_dir_dim_vec(InstVector& ivector, bool isPlus, int dir, int dim)
{
    declare_b_Spins(ivector);
    declare_misc(ivector);

    proj_ops *p_ops = (isPlus == true ? proj_ops_plus : proj_ops_minus);
    pack_face_vec(ivector, b_spinor, p_ops, 2*dim+dir);
}

void recons_add_face_from_dir_dim_vec(InstVector& ivector, bool compress12, bool isPlus, int dir, int dim, bool clover)
{
    declare_b_Spins(ivector);
    declare_ub_Spins(ivector);
    declare_outs(ivector);
    declare_u_gaus(ivector);
    declare_misc(ivector);

    if(clover) {
        declare_douts(ivector);
        declare_clover(ivector);
    }

    bool isBack = (dir == 0 ? true : false);
    recons_ops *rec_ops;

    if(clover) {
        rec_ops = (isPlus == isBack ? rec_plus_pbeta_ops : rec_minus_pbeta_ops);
    }
    else {
        rec_ops = (isPlus == isBack ? rec_plus_mbeta_ops : rec_minus_mbeta_ops);
    }

    recons_add_face_vec(ivector, compress12, isBack, rec_ops, dir, dim, clover);
}
