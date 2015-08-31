// Copyright (c) 2015 Matthias Noack (ma.noack.pr@gmail.com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include "common.hpp"

// unoptimised/readable reference implementation for correctness validation
void commutator_reference(complex_t* sigma_in, complex_t* sigma_out, complex_t* hamiltonian, size_t dim, size_t num_sigma, real_t hbar, real_t dt)
{
	const size_t size_sigma = dim * dim;
	real_t hdt = dt / hbar;

	// iterate over all sigma matrices
	#pragma omp parallel for
	for (size_t n = 0; n < num_sigma; ++n)
	{
		size_t sigma_id = n * size_sigma;
		// compute commutator term: i * dt / hbar * (hamiltonian * sigma - sigma * hamiltonian)
		#pragma novector
		for (size_t i = 0; i < dim; ++i)
		{
			#pragma novector
			for (size_t j = 0; j < dim; ++j)
			{
				complex_t tmp = 0.0; //(0.0, 0.0);
				#pragma novector
				for (size_t k = 0; k < dim; ++k)
				{
					tmp += hamiltonian[i * dim + k] * sigma_in[sigma_id + k * dim + j]
					     - sigma_in[sigma_id + i * dim + k] * hamiltonian[k * dim + j];
				}
				sigma_out[sigma_id + i * dim + j] -= complex_t(0.0,1.0) * hdt * tmp;
			}
		}
	}
}

