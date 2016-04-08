// Copyright (c) 2015 Matthias Noack (ma.noack.pr@gmail.com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include "kernel/common.cl"

__kernel __attribute__((vec_type_hint(real_vec_t)))
void commutator_ocl_manual_aosoa_direct(__global real_vec_t const* restrict sigma_in, 
                                        __global real_vec_t* restrict sigma_out, 
                                        __global real_t const* restrict hamiltonian, 
                                        const int num, const int dim,
                                        const real_t hbar, const real_t dt)
{
	// number of package to process == get_global_id(0)
	#define package_id (get_global_id(0) * dim * dim * 2)
	
	#define sigma_real(i, j) (package_id + 2 * (dim * (i) + (j)))
	#define sigma_imag(i, j) (package_id + 2 * (dim * (i) + (j)) + 1)
	
	#define ham_real(i, j) ((i) * dim + (j))
	#define ham_imag(i, j) (dim * dim + (i) * dim + (j))
	
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

