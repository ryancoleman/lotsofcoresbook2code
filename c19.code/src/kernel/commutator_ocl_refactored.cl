// Copyright (c) 2015 Matthias Noack (ma.noack.pr@gmail.com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include "kernel/common.cl"

// This kernel completely avoids the use of real_2_t as complex supplement.
// Instead, it only makes use of real_t, and a factor 2 in the indices.
// In comparison to the initial kernel, the index computations are put into
// macros to keep the arithmetic readable.
// The latter is reformulated into 8 lines that make the expected outcome of
// of 8 FMA instructions for the inner loop obvious. This allows to easily 
// analyse the assembly makes it potentially easier for the compiler to do the
// right thing.
__kernel
void commutator_ocl_refactored(__global real_t const* restrict sigma_in, 
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
	real_t hdt = dt / hbar;

	// compute commutator: -i * dt/hbar * (hamiltonian * sigma_in[sigma_id] - sigma_in[sigma_id] * hamiltonian)
	int i, j, k;
	for (i = 0; i < dim; ++i)
	{
		for (j = 0; j < dim; ++j)
		{
			real_t tmp_real = 0.0;
			real_t tmp_imag = 0.0;
			for (k = 0; k < dim; ++k)
			{
				tmp_real += hamiltonian[ham_real(i, k)] * sigma_in[sigma_real(k, j)];
				tmp_real -= sigma_in[sigma_real(i, k)] * hamiltonian[ham_real(k, j)];
				tmp_real -= hamiltonian[ham_imag(i, k)] * sigma_in[sigma_imag(k, j)];
				tmp_real += sigma_in[sigma_imag(i, k)] * hamiltonian[ham_imag(k, j)];
				tmp_imag += hamiltonian[ham_real(i, k)] * sigma_in[sigma_imag(k, j)];
				tmp_imag -= sigma_in[sigma_real(i, k)] * hamiltonian[ham_imag(k, j)];
				tmp_imag += hamiltonian[ham_imag(i, k)] * sigma_in[sigma_real(k, j)]; 
				tmp_imag -= sigma_in[sigma_imag(i, k)] * hamiltonian[ham_real(k, j)];
			}
			// multiply with -i dt/hbar
			sigma_out[sigma_real(i, j)] += hdt * tmp_imag;
			sigma_out[sigma_imag(i, j)] -= hdt * tmp_real;
		}
	}
}
