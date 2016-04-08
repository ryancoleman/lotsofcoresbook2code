// Copyright (c) 2015 Matthias Noack (ma.noack.pr@gmail.com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef kernel_hpp
#define kernel_hpp

#include "common.hpp"

void commutator_reference(complex_t* sigma_in, complex_t* sigma_out, complex_t* hamiltonian, size_t dim, size_t num_sigma, real_t hbar, real_t dt);

#define SCALAR_PARAMETERS real_t const* restrict sigma_in,    \
                           real_t* restrict sigma_out,         \
                           real_t const* restrict hamiltonian, \
                           const int num, const int dim,       \
                           const real_t hbar, const real_t dt
// auto: 

void commutator_omp_aosoa( SCALAR_PARAMETERS );

void commutator_omp_aosoa_constants( SCALAR_PARAMETERS );

void commutator_omp_aosoa_direct( SCALAR_PARAMETERS );

void commutator_omp_aosoa_constants_direct( SCALAR_PARAMETERS );

void commutator_omp_aosoa_constants_direct_perm( SCALAR_PARAMETERS );

// no OpenCL equivalent:
void commutator_omp_aosoa_constants_direct_perm2to3( SCALAR_PARAMETERS );

void commutator_omp_aosoa_constants_direct_perm2to5( SCALAR_PARAMETERS );

# define VECTOR_PARAMETERS real_vec_t const* restrict sigma_in, \
                           real_vec_t* restrict sigma_out,      \
                           real_t const* restrict hamiltonian,  \
                           const int num, const int dim,        \
                           const real_t hbar, const real_t dt
// manual:
void commutator_omp_manual_aosoa( VECTOR_PARAMETERS );

void commutator_omp_manual_aosoa_constants( VECTOR_PARAMETERS );

void commutator_omp_manual_aosoa_constants_perm( VECTOR_PARAMETERS );

void commutator_omp_manual_aosoa_direct( VECTOR_PARAMETERS );

void commutator_omp_manual_aosoa_constants_direct( VECTOR_PARAMETERS );

void commutator_omp_manual_aosoa_constants_direct_perm( VECTOR_PARAMETERS );

// no OpenCL equivalent:
void commutator_omp_manual_aosoa_constants_direct_perm4to5( VECTOR_PARAMETERS );

void commutator_omp_manual_aosoa_constants_direct_unrollhints( VECTOR_PARAMETERS );

void commutator_omp_manual_aosoa_constants_direct_perm_unrollhints( VECTOR_PARAMETERS );

#undef SCALAR_PARAMETERS
#undef VECTOR_PARAMETERS
#endif // kernel_hpp

