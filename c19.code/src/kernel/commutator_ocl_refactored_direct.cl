// Copyright (c) 2015 Matthias Noack (ma.noack.pr@gmail.com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include "kernel/common.cl"

// Same as commutator_ocl_refactored.cl, but with direct stores and pre-scaled
// Hamiltonian ( * dt / hbar).
__kernel
void commutator_ocl_refactored_direct(__global real_t const* restrict sigma_in, 
                                      __global real_t* restrict sigma_out, 
                                      __global real_t const* restrict hamiltonian,
                                      const int num, const int dim,
                                      const real_t hbar, const real_t dt)
{
	// number of package to process == get_global_id(0)
	#define sigma_real(i, j) (sigma_id + 2 * ((i) * dim + (j)))
	#define sigma_imag(i, j) (sigma_id + 2 * ((i) * dim + (j)) + 1)
	
	#define ham_real(i, j) (2 * ((i) * dim + (j)))
	#define ham_imag(i, j) (2 * ((i) * dim + (k)) + 1)

	int sigma_id = get_global_id(0) * dim * dim * 2;

	// compute commutator: (hamiltonian * sigma_in[sigma_id] - sigma_in[sigma_id] * hamiltonian)
	int i, j, k;
	for (i = 0; i < dim; ++i)
	{
		for (j = 0; j < dim; ++j)
		{
			for (k = 0; k < dim; ++k)
			{
				sigma_out[sigma_imag(i, j)] -= hamiltonian[ham_real(i, k)] * sigma_in[sigma_real(k, j)];
				sigma_out[sigma_imag(i, j)] += sigma_in[sigma_real(i, k)] * hamiltonian[ham_real(k, j)];
				sigma_out[sigma_imag(i, j)] += hamiltonian[ham_imag(i, k)] * sigma_in[sigma_imag(k, j)];
				sigma_out[sigma_imag(i, j)] -= sigma_in[sigma_imag(i, k)] * hamiltonian[ham_imag(k, j)];
				sigma_out[sigma_real(i, j)] += hamiltonian[ham_real(i, k)] * sigma_in[sigma_imag(k, j)];
				sigma_out[sigma_real(i, j)] -= sigma_in[sigma_real(i, k)] * hamiltonian[ham_imag(k, j)];
				sigma_out[sigma_real(i, j)] += hamiltonian[ham_imag(i, k)] * sigma_in[sigma_real(k, j)];
				sigma_out[sigma_real(i, j)] -= sigma_in[sigma_imag(i, k)] * hamiltonian[ham_real(k, j)];
			}
		}
	}
}
