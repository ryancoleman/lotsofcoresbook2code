// -*- C++ -*-

/*! \file
 * \brief Primitive classes
 *
 * Primitives are the various types on the fibers at the lattice sites
 */


#ifndef QDP_PRIMITIVE_H
#define QDP_PRIMITIVE_H


/*! \defgroup fiber Fiber only types and operations
 * \ingroup fiberbundle
 *
 * Primitives are the various types on the fibers at the lattice sites.
 *
 * The primitive indices, including Reality (also known as complex or real),
 * is represented as a tensor product over various vector spaces. Different
 * kinds of object can transform in those vector spaces, like Scalar, Vector, and
 * Matrix.
 */

#include "qdp_primscalar.h"
#include "qdp_primmatrix.h"
#include "qdp_primvector.h"
#include "qdp_primseed.h"
#include "qdp_primcolormat.h"
#include "qdp_primcolorvec.h"
#include "qdp_primgamma.h"
#include "qdp_primspinmat.h"
#include "qdp_primspinvec.h"

#endif
