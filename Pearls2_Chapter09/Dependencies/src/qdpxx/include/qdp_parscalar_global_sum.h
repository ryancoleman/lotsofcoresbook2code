#ifndef QDP_GLOBAL_SUM_H
#define QDP_GLOBAL_SUM_H
 
#include <qmp.h>
 
namespace QDPGlobalSums {
  QMP_status_t QDP_sum_int(int *i);
  QMP_status_t QDP_sum_float_array(float *x, int length);
  QMP_status_t QDP_sum_double_array(double *x, int length);
};
 
#endif
