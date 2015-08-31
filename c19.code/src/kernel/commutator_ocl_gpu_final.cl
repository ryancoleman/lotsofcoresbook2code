// Copyright (c) 2015 Florian Wende (wende@zib.de)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include "kernel/common.cl"

#define id_2d_to_1d(i,j) ((i) * DIM + (j))
#define sigma_id(i,j,m) ((m) * DIM * DIM + ((i) * DIM + (j)))

#define MIN(X,Y) ((X) < (Y) ? (X) : (Y))

// This kernel is specific to using the real_2_t data type for internal computations
__kernel
void commutator_ocl_gpu_final(__global real_2_t const* restrict sigma_in,
                              __global real_2_t* restrict sigma_out,
                              __global real_2_t const* restrict hamiltonian,
                              const int num, const int dim,
                              const real_t hbar, const real_t dt)
{
	// Determine matrix index (i,j) this work item is responsible for
	int ij = get_local_id(0);
	int i = ij / DIM; // Matrix index 'i' to be processed by this work item in any of 'start -> stop' matrices
	int j = ij % DIM; // Matrix index 'j' to be processed by this work item in any of 'start -> stop' matrices

	// Determine working set : Each work item participates in processing CHUNK_SIZE matrices : 'start -> stop'
	int sub_group_id = get_local_id(1); // Local matrix ID within work group
	int start = get_group_id(0) * NUM_SUB_GROUPS * CHUNK_SIZE + sub_group_id * CHUNK_SIZE; // Global matrix ID : start
	int stop = MIN(num, start + CHUNK_SIZE); // Global matrix ID : stop

	// Local variables
	real_2_t snew1_ij, snew2_ij;
	real_2_t s1, s2;

	// Local memory: shared between all work items in the same work group
	// Remark 1: 2-way shared memory bank conflicts will occur for real_t = double
	// Remark 2: real parts and imaginary parts are stored separately to avoid 4-way bank conflicts in case of real_2_t = double2
	__local real_t ham_local_real[DIM * DIM]; // Hamiltonian: real part
	__local real_t ham_local_imag[DIM * DIM]; // Hamiltonian: imaginary part
	__local real_t sigma_local_real[2][NUM_SUB_GROUPS][DIM * DIM]; // Input sigma matrix: real part (2 matrices are processed at once)
	__local real_t sigma_local_imag[2][NUM_SUB_GROUPS][DIM * DIM]; // Input sigma matrix: imaginary part (2 matrices are processed at once)

	// Load Hamiltonian into local memory: only the first sub-group participates
	if (ij < (DIM * DIM) && sub_group_id == 0)
	{
		const real_2_t h = hamiltonian[ij];
		ham_local_real[ij] = h.x;
		ham_local_imag[ij] = h.y;
	}

	// Process all CHUNK_SIZE matrices: two matrices are processed at once (therefore increment 2)
	for (int m = start; m < stop; m += 2)
	{
		barrier(CLK_LOCAL_MEM_FENCE);
		if (ij < (DIM * DIM)) 
		{ // Load input sigma matrix into local memory: only threads with valid IDs participate
			s1 = sigma_in[sigma_id(i, j, m)]; // Real and imaginary part of matrix 'm', element (i,j)
			sigma_local_real[0][sub_group_id][ij] = s1.x;
			sigma_local_imag[0][sub_group_id][ij] = s1.y;

			s2 = sigma_in[sigma_id(i, j, m + 1)]; // Real and imaginary part of matrix 'm+1', element (i,j)
			sigma_local_real[1][sub_group_id][ij] = s2.x;
			sigma_local_imag[1][sub_group_id][ij] = s2.y;

			s1 = sigma_out[sigma_id(i, j, m)]; // Prefetch real and imaginary part of output sigma matrix 'm', element (i,j)
			snew1_ij.x = s1.x;
			snew2_ij.x = s1.y;

			s2 = sigma_out[sigma_id(i, j, m + 1)]; // Prefetch real and imaginary part of output sigma matrix 'm+1', element (i,j)
			snew1_ij.y = s2.x;
			snew2_ij.y = s2.y;
		}
		barrier(CLK_LOCAL_MEM_FENCE);

		if (ij < (DIM * DIM))
		{
			// Compute commutator: [H,sigma] = H * sigma - sigma * H <=> [H,sigma]_ij = \sum_k ( H_ik * sigma_kj - sigma_ik * H_kj )
			for (int k = 0; k < DIM; ++k)
			{
				const int ik = id_2d_to_1d(i, k);
				const int kj = id_2d_to_1d(k, j);

				// Reassemble real_2_t elements from local memory: 'vector processing' gives better performance here
				s1 = (real_2_t)(sigma_local_real[0][sub_group_id][kj], sigma_local_real[1][sub_group_id][kj]);
				s2 = (real_2_t)(sigma_local_imag[0][sub_group_id][kj], sigma_local_imag[1][sub_group_id][kj]);
				snew1_ij += ham_local_real[ik] * s2;
				snew1_ij += ham_local_imag[ik] * s1;
				snew2_ij -= ham_local_real[ik] * s1;
				snew2_ij += ham_local_imag[ik] * s2;

				// Reassemble real_2_t elements from local memory: 'vector processing' gives better performance here
				s1 = (real_2_t)(sigma_local_real[0][sub_group_id][ik], sigma_local_real[1][sub_group_id][ik]);
				s2 = (real_2_t)(sigma_local_imag[0][sub_group_id][ik], sigma_local_imag[1][sub_group_id][ik]);
				snew1_ij -= ham_local_real[kj] * s2;
				snew1_ij += ham_local_imag[kj] * s1;
				snew2_ij += ham_local_real[kj] * s1;
				snew2_ij -= ham_local_imag[kj] * s2;
			}

			// Write output sigma matrices 'm' and 'm+1', element (i,j)
			sigma_out[sigma_id(i, j, m)] = (double2)(snew1_ij.x, snew2_ij.x);
			sigma_out[sigma_id(i, j, m + 1)] = (double2)(snew1_ij.y, snew2_ij.y);
		}
	}
}

