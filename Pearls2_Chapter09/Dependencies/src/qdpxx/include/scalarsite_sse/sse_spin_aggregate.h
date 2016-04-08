#ifndef SSE_SPIN_AGGREGATE
#define SSE_SPIN_AGGREGATE

// The inline project, recon and add recon ops.
#include "scalarsite_sse/sse_spin_proj_inlines.h"
#include "scalarsite_sse/sse_spin_recon_inlines.h"

// These are now just hooks that call the above
//#include "scalarsite_sse/sse_spin_proj.h"
//#include "scalarsite_sse/sse_spin_recon.h"
//#include "scalarsite_sse/sse_fused_spin_proj.h"
//#include "scalarsite_sse/sse_fused_spin_recon.h"

// These are the 'evaluates'
#include "scalarsite_sse/qdp_sse_spin_evaluates.h"
#include "scalarsite_sse/qdp_sse_fused_spin_proj_evaluates.h"
#include "scalarsite_sse/qdp_sse_fused_spin_recon_evaluates.h"

#endif
