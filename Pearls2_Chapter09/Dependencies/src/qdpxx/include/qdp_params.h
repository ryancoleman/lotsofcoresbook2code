// -*- C++ -*-

/*! @file
 * @brief Fundamental parameters
 */

#ifndef QDP_PARAMS_H
#define QDP_PARAMS_H

namespace QDP {


/*! @defgroup params Fundamental parameters for QDP
 *
 * The user can change the compile time number of dimensions,
 * colors and spin components -- note these are now determined
 * by configure
 *
 * @{
 */

#include <qdp_config.h>

const int Nd = QDP_ND;
const int Nc = QDP_NC;
const int Ns = QDP_NS;

/*! @} */  // end of group params

} // namespace QDP

#endif


