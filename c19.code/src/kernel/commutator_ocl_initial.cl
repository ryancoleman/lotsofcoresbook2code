// Copyright (c) 2015 Matthias Noack (ma.noack.pr@gmail.com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include "kernel/common.cl"

__kernel
void commutator_ocl_initial(__global real_2_t const* restrict sigma_in, 
                            __global real_2_t* restrict sigma_out, 
                            __global real_2_t const* restrict hamiltonian, 
                            const int num, const int dim,
                            const real_t hbar, const real_t dt)
{
	int sigma_id = get_global_id(0) * dim * dim;
	int i, j, k;
	real_t hdt = dt / hbar;

	// compute commutator: -i * dt/hbar * (hamiltonian * sigma_in[sigma_id] - sigma_in[sigma_id] * hamiltonian)
	for (i = 0; i < dim; ++i)
	{
		for (j = 0; j < dim; ++j)
		{
			real_2_t tmp;
			tmp.x = 0.0;
			tmp.y = 0.0;
			for (k = 0; k < dim; ++k)
			{
				tmp.x += (hamiltonian[i * dim + k].x * sigma_in[sigma_id + k * dim + j].x - sigma_in[sigma_id + i * dim + k].x * hamiltonian[k * dim + j].x);
				tmp.x -= (hamiltonian[i * dim + k].y * sigma_in[sigma_id + k * dim + j].y - sigma_in[sigma_id + i * dim + k].y * hamiltonian[k * dim + j].y);
				tmp.y += (hamiltonian[i * dim + k].x * sigma_in[sigma_id + k * dim + j].y - sigma_in[sigma_id + i * dim + k].x * hamiltonian[k * dim + j].y);
				tmp.y += (hamiltonian[i * dim + k].y * sigma_in[sigma_id + k * dim + j].x - sigma_in[sigma_id + i * dim + k].y * hamiltonian[k * dim + j].x);
			}
			// multiply with -i * dt / hbar
			sigma_out[sigma_id + i * dim + j].x += hdt * tmp.y;
			sigma_out[sigma_id + i * dim + j].y -= hdt * tmp.x;
		}
	}
}
