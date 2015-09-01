#ifndef GENERIC_SPIN_AGGREGATE_H
#define GENERIC_SPIN_AGGREGATE_H

// Include this first to define the mat mult ops (possibly assembler?)
#include "scalarsite_generic/generic_mv_switchbox.h"

// The inline project, recon and add recon ops.
#include "scalarsite_generic/generic_spin_proj_inlines.h"
#include "scalarsite_generic/generic_spin_recon_inlines.h"

// These are now just hooks that call the above
// #include "scalarsite_generic/generic_spin_proj.h"
// #include "scalarsite_generic/generic_spin_recon.h"
// #include "scalarsite_generic/generic_fused_spin_proj.h"
// #include "scalarsite_generic/generic_fused_spin_recon.h"

// These are the 'evaluates'
#include "scalarsite_generic/qdp_generic_spin_project_evaluates.h"
#include "scalarsite_generic/qdp_generic_fused_spin_proj_evaluates.h"
#include "scalarsite_generic/qdp_generic_fused_spin_recon_evaluates.h"

#endif
