#ifndef QDP_SSE_LINALG_MM_SU3_DOUBLE_H
#define QDP_SSE_LINALG_MM_SU3_DOUBLE_H

#include "qdp_precision.h"
namespace QDP { 

  /* M = a*M  a is scalar */
  void ssed_m_eq_scal_m(REAL64* m2, REAL64* a, REAL64 *m1, int n_mat);

  /* M *= a,  a is a scalar */
  void ssed_m_muleq_scal(REAL64* m, REAL64* a, int n_mat);

  /* M2 += M1 */
  void ssed_m_peq_m(REAL64* m2, REAL64* m1, int n_mat);

  /* M2 -= M1 */
  void ssed_m_meq_m(REAL64* m2, REAL64* m1, int nmat);

  /* M2 += adj(M1) */
  void ssed_m_peq_h(REAL64* m2, REAL64* m1, int n_mat);

  /* M2 -= adj(M1) */
  void ssed_m_meq_h(REAL64* m2, REAL64* m1, int n_mat);

  /* M3 = M1*M2 */
  void ssed_m_eq_mm(REAL64* m3, REAL64* m1, REAL64* m2, int n_mat);
  void ssed_m_eq_mm_u(REAL64* m3, REAL64* m1, REAL64* m2, int n_mat);
 
 /* M3 += a M1*M2 */
  void ssed_m_peq_amm(REAL64* m3, REAL64* a, REAL64* m1, REAL64* m2, int n_mat); 
  void ssed_m_peq_amm_u(REAL64* m3, REAL64* a, REAL64* m1, REAL64* m2, int n_mat); 

  /* M3 = M1*adj(M2) */
  void ssed_m_eq_mh(REAL64* m3, REAL64* m1, REAL64* m2, int n_mat);
  void ssed_m_eq_mh_u(REAL64* m3, REAL64* m1, REAL64* m2, int n_mat);
  
  /* M3 += a M1*adj(M2) */
  void ssed_m_peq_amh(REAL64* m3, REAL64* a, REAL64* m1, REAL64* m2, int n_mat);
  void ssed_m_peq_amh_u(REAL64* m3, REAL64* a, REAL64* m1, REAL64* m2, int n_mat);

  
  /* M3 = adj(M1)*M2 */
  void ssed_m_eq_hm(REAL64* m3, REAL64* m1, REAL64* m2, int n_mat);
  void ssed_m_eq_hm_u(REAL64* m3, REAL64* m1, REAL64* m2, int n_mat);

  /* M3 += a adj(M1)*M2 */
  void ssed_m_peq_ahm(REAL64* m3, REAL64*a, REAL64* m1, REAL64* m2, int n_mat);
  void ssed_m_peq_ahm_u(REAL64* m3, REAL64*a, REAL64* m1, REAL64* m2, int n_mat);
  
  /* M3 = adj(M1)*adj(M2) */
  void ssed_m_eq_hh(REAL64* m3, REAL64* m1, REAL64* m2, int n_mat);
  void ssed_m_eq_hh_u(REAL64* m3, REAL64* m1, REAL64* m2, int n_mat);

  /* M3 += a adj(M1)*adj(M2) */
  void ssed_m_peq_ahh(REAL64* m3, REAL64*a, REAL64* m1, REAL64* m2, int n_mat);
  void ssed_m_peq_ahh_u(REAL64* m3, REAL64*a, REAL64* m1, REAL64* m2, int n_mat);


  


}
#endif
