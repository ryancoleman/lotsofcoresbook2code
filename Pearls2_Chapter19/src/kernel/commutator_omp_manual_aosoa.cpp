// Copyright (c) 2015 Matthias Noack (ma.noack.pr@gmail.com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include "common.hpp"

void commutator_omp_manual_aosoa(real_vec_t const* restrict sigma_in, 
                                 real_vec_t* restrict sigma_out, 
                                 real_t const* restrict hamiltonian, 
                                 const int num, const int dim,
                                 const real_t hbar, const real_t dt)
{
	// OpenCL work-groups are mapped to threads
	#pragma omp parallel for
	#pragma novector // NOTE: we do not want any implicit vectorisation in this kernel
	for (int global_id = 0; global_id < (num / VEC_LENGTH); ++global_id)
	{
		// original OpenCL kernel begins here
		#define package_id (global_id * dim * dim * 2)

		#define sigma_real(i, j) (package_id + 2*(dim * i + j))
		#define sigma_imag(i, j) (package_id + 2*(dim * i + j) + 1)

		#define ham_real(i, j) (i * dim + j)
		#define ham_imag(i, j) (dim * dim + i * dim + j)
	
		// compute commutator: (hamiltonian * sigma_in[sigma_id] - sigma_in[sigma_id] * hamiltonian)
		int i, j, k;
		for (i = 0; i < dim; ++i)
		{
			for (j = 0; j < dim; ++j)
			{
#ifdef USE_INITZERO
				real_vec_t tmp_real(0.0);
				real_vec_t tmp_imag(0.0);
# else
				real_vec_t tmp_real = sigma_out[sigma_real(i, j)];
				real_vec_t tmp_imag = sigma_out[sigma_imag(i, j)];
#endif
				for (k = 0; k < dim; ++k)
				{
					tmp_imag -= sigma_in[sigma_real(k, j)] * hamiltonian[ham_real(i, k)];
					tmp_imag += sigma_in[sigma_real(i, k)] * hamiltonian[ham_real(k, j)];
					tmp_imag += sigma_in[sigma_imag(k, j)] * hamiltonian[ham_imag(i, k)];
					tmp_imag -= sigma_in[sigma_imag(i, k)] * hamiltonian[ham_imag(k, j)];
					tmp_real += sigma_in[sigma_imag(k, j)] * hamiltonian[ham_real(i, k)];
					tmp_real -= sigma_in[sigma_real(i, k)] * hamiltonian[ham_imag(k, j)];
					tmp_real += sigma_in[sigma_real(k, j)] * hamiltonian[ham_imag(i, k)];
					tmp_real -= sigma_in[sigma_imag(i, k)] * hamiltonian[ham_real(k, j)];
				}
#ifdef USE_INITZERO
				sigma_out[sigma_real(i, j)] += tmp_real;
				sigma_out[sigma_imag(i, j)] += tmp_imag;
#else
				sigma_out[sigma_real(i, j)] = tmp_real;
				sigma_out[sigma_imag(i, j)] = tmp_imag;
#endif
			}
		}
	}
}

