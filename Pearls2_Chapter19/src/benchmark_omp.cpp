// Copyright (c) 2015 Matthias Noack (ma.noack.pr@gmail.com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <iostream>

#include <cstring> // memcpy
#include <cmath>

#include "ham/util/time.hpp" // ham::util::time

#include "common.hpp"
#include "kernel/kernel.hpp"

using namespace ham::util;

int main(void)
{
	print_compile_config(std::cerr);

	// constants
	const size_t dim = DIM;
	const size_t num = NUM;
	const real_t hbar = 1.0 / std::acos(-1.0); // == 1 / Pi
	const real_t dt = 1.0e-3; 

	real_t deviation = 0.0;

	// allocate memory
	size_t size_hamiltonian = dim * dim;
	size_t size_sigma = size_hamiltonian * num;
	size_t size_sigma_byte = sizeof(complex_t) * size_sigma;

	complex_t* hamiltonian = allocate_aligned<complex_t>(size_hamiltonian);
	complex_t* sigma_in = allocate_aligned<complex_t>(size_sigma);
	complex_t* sigma_out = allocate_aligned<complex_t>(size_sigma);
	complex_t* sigma_reference = allocate_aligned<complex_t>(size_sigma);
	complex_t* sigma_reference_transformed = allocate_aligned<complex_t>(size_sigma);

	// initialise memory
	initialise_hamiltonian(hamiltonian, dim);
	initialise_sigma(sigma_in, sigma_out, dim, num);

	// print output header
	std::cout << "name\t" << time::statistics::header_string() << std::endl;
	
	// perform reference computation for correctness analysis
	benchmark_kernel(
		[&]() // lambda expression
		{
			commutator_reference(sigma_in, sigma_out, hamiltonian, dim, num, hbar, dt);
		},
		"commutator_reference",
		NUM_ITERATIONS,
		NUM_WARMUP);

	// copy reference results
	std::memcpy(sigma_reference, sigma_out, size_sigma_byte);

	// Lambda to: transform memory, benchmark, compare results
	auto benchmark = [&](std::function<void()> kernel,
	                     std::string name,
	                     decltype(&transform_matrices_aos_to_aosoa) transformation_sigma,
	                     bool scale_hamiltonian,
	                     decltype(&transform_matrix_aos_to_soa) transformation_hamiltonian)
	{
		initialise_hamiltonian(hamiltonian, dim);
		if (scale_hamiltonian) 
			transform_matrix_scale_aos(hamiltonian, dim, dt / hbar); // pre-scale hamiltonian
		if (transformation_hamiltonian)
			transformation_hamiltonian(hamiltonian, dim);	
	
		initialise_sigma(sigma_in, sigma_out, dim, num);
		std::memcpy(sigma_reference_transformed, sigma_reference, size_sigma_byte);
		// transform memory layout if a transformation is specified
		if (transformation_sigma)
		{
			// transform reference for comparison
			transformation_sigma(sigma_reference_transformed, dim, num, VEC_LENGTH);
			// tranform sigma
			transformation_sigma(sigma_in, dim, num, VEC_LENGTH);
		}
		
		benchmark_kernel(kernel, name, NUM_ITERATIONS, NUM_WARMUP);
		
		// compute deviation from reference	(small deviations are expected)
		deviation = compare_matrices(sigma_out, sigma_reference_transformed, dim, num);
		std::cerr << "Deviation:\t" << deviation << std::endl;
	};
	
	
#define SCALAR_ARGUMENTS reinterpret_cast<real_t*>(sigma_in),    \
			 reinterpret_cast<real_t*>(sigma_out),   \
			 reinterpret_cast<real_t*>(hamiltonian), \
			 num, dim, 0.0, 0.0
	// BENCHMARK
	benchmark(
		[&]() // lambda expression
		{
			commutator_omp_aosoa( SCALAR_ARGUMENTS );
		},
		"commutator_omp_aosoa",
		&transform_matrices_aos_to_aosoa, SCALE_HAMILT, &transform_matrix_aos_to_soa);


	// BENCHMARK
	benchmark(
		[&]() // lambda expression
		{
			commutator_omp_aosoa_constants( SCALAR_ARGUMENTS );
		},
		"commutator_omp_aosoa_constants",
		&transform_matrices_aos_to_aosoa, SCALE_HAMILT, &transform_matrix_aos_to_soa);

	// BENCHMARK
	benchmark(
		[&]() // lambda expression
		{
			commutator_omp_aosoa_direct( SCALAR_ARGUMENTS );
		},
		"commutator_omp_aosoa_direct",
		&transform_matrices_aos_to_aosoa, SCALE_HAMILT, &transform_matrix_aos_to_soa);

	// BENCHMARK
	benchmark(
		[&]() // lambda expression
		{
			commutator_omp_aosoa_constants_direct( SCALAR_ARGUMENTS );
		},
		"commutator_omp_aosoa_constants_direct",
		&transform_matrices_aos_to_aosoa, SCALE_HAMILT, &transform_matrix_aos_to_soa);

	// BENCHMARK
	benchmark(
		[&]() // lambda expression
		{
			commutator_omp_aosoa_constants_direct_perm( SCALAR_ARGUMENTS );
		},
		"commutator_omp_aosoa_constants_direct_perm",
		&transform_matrices_aos_to_aosoa, SCALE_HAMILT, &transform_matrix_aos_to_soa);

	// BENCHMARK
	benchmark(
		[&]() // lambda expression
		{
			commutator_omp_aosoa_constants_direct_perm2to3( SCALAR_ARGUMENTS );
		},
		"commutator_omp_aosoa_constants_direct_perm2to3",
		&transform_matrices_aos_to_aosoa, SCALE_HAMILT, &transform_matrix_aos_to_soa);

	// BENCHMARK
	benchmark(
		[&]() // lambda expression
		{
			commutator_omp_aosoa_constants_direct_perm2to5( SCALAR_ARGUMENTS );
		},
		"commutator_omp_aosoa_constants_direct_perm2to5",
		&transform_matrices_aos_to_aosoa, SCALE_HAMILT, &transform_matrix_aos_to_soa);		
		

	// manually vectorised kernels

#define VECTOR_ARGUMENTS reinterpret_cast<real_vec_t*>(sigma_in),  \
			 reinterpret_cast<real_vec_t*>(sigma_out), \
			 reinterpret_cast<real_t*>(hamiltonian),   \
			 num, dim, 0.0, 0.0
	// BENCHMARK: 
	benchmark(
		[&]() // lambda expression
		{
			commutator_omp_manual_aosoa( VECTOR_ARGUMENTS );
		},
		"commutator_omp_manual_aosoa",
		&transform_matrices_aos_to_aosoa, SCALE_HAMILT, &transform_matrix_aos_to_soa);

	// BENCHMARK: 
	benchmark(
		[&]() // lambda expression
		{
			commutator_omp_manual_aosoa_constants( VECTOR_ARGUMENTS );
		},
		"commutator_omp_manual_aosoa_constants",
		&transform_matrices_aos_to_aosoa, SCALE_HAMILT, &transform_matrix_aos_to_soa);

	// BENCHMARK: 
	benchmark(
		[&]() // lambda expression
		{
			commutator_omp_manual_aosoa_constants_perm( VECTOR_ARGUMENTS );
		},
		"commutator_omp_manual_aosoa_constants_perm",
		&transform_matrices_aos_to_aosoa, SCALE_HAMILT, &transform_matrix_aos_to_soa);

	// BENCHMARK: 
	benchmark(
		[&]() // lambda expression
		{
			commutator_omp_manual_aosoa_direct( VECTOR_ARGUMENTS );
		},
		"commutator_omp_manual_aosoa_direct",
		&transform_matrices_aos_to_aosoa, SCALE_HAMILT, &transform_matrix_aos_to_soa);

	// BENCHMARK: 
	benchmark(
		[&]() // lambda expression
		{
			commutator_omp_manual_aosoa_constants_direct( VECTOR_ARGUMENTS );
		},
		"commutator_omp_manual_aosoa_constants_direct",
		&transform_matrices_aos_to_aosoa, SCALE_HAMILT, &transform_matrix_aos_to_soa);

	// BENCHMARK: 
	benchmark(
		[&]() // lambda expression
		{
			commutator_omp_manual_aosoa_constants_direct_perm( VECTOR_ARGUMENTS );
		},
		"commutator_omp_manual_aosoa_constants_direct_perm",
		&transform_matrices_aos_to_aosoa, SCALE_HAMILT, &transform_matrix_aos_to_soa);

	// BENCHMARK: 
	benchmark(
		[&]() // lambda expression
		{
			commutator_omp_manual_aosoa_constants_direct_unrollhints( VECTOR_ARGUMENTS );
		},
		"commutator_omp_manual_aosoa_constants_direct_unrollhints",
		&transform_matrices_aos_to_aosoa, SCALE_HAMILT, &transform_matrix_aos_to_soa);		
		
	// BENCHMARK: 
	benchmark(
		[&]() // lambda expression
		{
			commutator_omp_manual_aosoa_constants_direct_perm_unrollhints( VECTOR_ARGUMENTS );
		},
		"commutator_omp_manual_aosoa_constants_direct_perm_unrollhints",
		&transform_matrices_aos_to_aosoa, SCALE_HAMILT, &transform_matrix_aos_to_soa);		
		


	delete hamiltonian;
	delete sigma_in;
	delete sigma_out;
	delete sigma_reference;

	return 0;
}

