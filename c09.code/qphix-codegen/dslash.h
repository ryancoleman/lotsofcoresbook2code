#include <cstdio>
#include <cstdlib>

using namespace std;

#include "instructions.h"
#include "data_types.h"


typedef struct {
    const char *name;
    int s[2][2];
    void (*CVecFunc[2])(InstVector& ivector, FVec *r, FVec *s1, FVec *s2, string &mask);
//	void (*CVecFunc2)(InstVector& ivector, FVec *r, FVec *s1, FVec *s2, string &mask);
} proj_ops;

typedef struct {
    const char *name;
    int s2;
    int s3;
    void (*CVecFuncTop2)(InstVector& ivector, FVec *r, FVec *s1, FVec *s2, string &mask);
    void (*CVecFunc1)(InstVector& ivector, FVec *r, FVec *s1, FVec *s2, string &mask);
    void (*CVecFunc2)(InstVector& ivector, FVec *r, FVec *s1, FVec *s2, string &mask);
} recons_ops;


// Uses psi[][] as temp and b_spinor[][] as implicit output
void project(InstVector& ivector, string base, string offset, proj_ops& ops, bool isFace, string mask, int dir);
// Serial Spin version
void project(InstVector& ivector, string base, string offset, proj_ops& ops, bool isFace, string mask, int dir, int s);
void recons_add(InstVector& ivector, recons_ops& ops, FVec outspinor[4][3][2], string &mask);
// Serial Spin version
void recons_add(InstVector& ivector, recons_ops& ops, FVec outspinor[4][3][2], string &mask, int s);
void zeroResult(InstVector& ivector, FVec *outspinor);
void clover_term(InstVector& ivector, FVec in_spinor[4][3][2], bool face, string _mask="");
void achiResult(InstVector& ivector, bool clover);
void loadGaugeDir(InstVector& ivector, int dir, bool compress12);
void matMultVec(InstVector& ivector, bool adjMul, int s);
void matMultVec(InstVector& ivector, bool adjMul);
void dslash_plain_body(InstVector& ivector, bool compress12, bool clover, bool isPlus);
// ***** ------- a chi - b D psi versions
void dslash_achimbdpsi_body(InstVector& ivector, bool compress12, bool clover, bool isPlus);
void pack_face_to_dir_dim_vec(InstVector& ivector, bool isPlus, int dir, int dim);
void recons_add_face_from_dir_dim_vec(InstVector& ivector, bool compress12, bool isPlus, int dir, int dim, bool clover);

void dslash_body(InstVector& ivector, bool compress12, proj_ops *ops, recons_ops *rec_ops_bw, recons_ops *rec_ops_fw, FVec outspinor[4][3][2]);
void pack_face_vec(InstVector& ivector, FVec spinor[2][3][2], proj_ops proj[], int dir);
void recons_add_face_vec(InstVector& ivector, bool compress12, bool adjMul, recons_ops rops[], int dir, int dim, bool clover);
