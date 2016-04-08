/* Generic MV Switchbox */
#ifndef GENERIC_MV_SWITCHBOX_H
#define GENERIC_MV_SWITCHBOX_H

#include "scalarsite_generic/generic_mat_vec.h"

#define _inline_mult_su3_mat_vec(aa,bb,cc) \
{\
  _inline_generic_mult_su3_mat_vec(aa,bb,cc) \
}

#define _inline_mult_adj_su3_mat_vec(aa,bb,cc) \
{\
 _inline_generic_mult_adj_su3_mat_vec(aa,bb,cc) \
}

#endif
