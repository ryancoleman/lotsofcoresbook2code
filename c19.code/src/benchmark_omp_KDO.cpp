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
#define USE_FASTREFERENCE
#ifdef USE_FASTREFERENCE
		1,
		0);
#else
		NUM_ITERATIONS,
		NUM_WARMUP);
#endif

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
	
	
	// BENCHMARK: run optimised automatically vectorised kernel
	benchmark(
		[&]() // lambda expression
		{
			commutator_omp_auto_direct(
				reinterpret_cast<real_t*>(sigma_in),
				reinterpret_cast<real_t*>(sigma_out),
				reinterpret_cast<real_t*>(hamiltonian),
				num, 0, 0.0, 0.0);
		},
		"commutator_omp_auto_direct",
		&transform_matrices_aos_to_aosoa, SCALE_HAMILT, &transform_matrix_aos_to_soa);
	
	// BENCHMARK: run optimised automatically vectorised kernel
	benchmark(
		[&]() // lambda expression
		{
			commutator_omp_auto_direct_perm2to3(
				reinterpret_cast<real_t*>(sigma_in),
				reinterpret_cast<real_t*>(sigma_out),
				reinterpret_cast<real_t*>(hamiltonian),
				num, 0, 0.0, 0.0);
		},
		"commutator_omp_auto_direct_perm2to3",
		&transform_matrices_aos_to_aosoa, SCALE_HAMILT, &transform_matrix_aos_to_soa);
	
	// BENCHMARK: run optimised automatically vectorised kernel
	benchmark(
		[&]() // lambda expression
		{
			commutator_omp_auto_direct_perm2to5(
				reinterpret_cast<real_t*>(sigma_in),
				reinterpret_cast<real_t*>(sigma_out),
				reinterpret_cast<real_t*>(hamiltonian),
				num, 0, 0.0, 0.0);
		},
		"commutator_omp_auto_direct_perm2to5",
		&transform_matrices_aos_to_aosoa, SCALE_HAMILT, &transform_matrix_aos_to_soa);

	// BENCHMARK: run optimised manually vectorised kernel
	benchmark(
		[&]() // lambda expression
		{
			commutator_omp_manual_aosoa_constants_direct(
				reinterpret_cast<real_vec_t*>(sigma_in),
				reinterpret_cast<real_vec_t*>(sigma_out),
				reinterpret_cast<real_t*>(hamiltonian),
				num, 0, 0.0, 0.0);
		},
		"commutator_omp_manual_aosoa_constants_direct",
		&transform_matrices_aos_to_aosoa, SCALE_HAMILT, &transform_matrix_aos_to_soa);

	// BENCHMARK: run optimised manually vectorised kernel
	benchmark(
		[&]() // lambda expression
		{
			commutator_omp_manual_direct(
				reinterpret_cast<real_vec_t*>(sigma_in),
				reinterpret_cast<real_vec_t*>(sigma_out),
				reinterpret_cast<real_t*>(hamiltonian),
				num, 0, 0.0, 0.0);
		},
		"commutator_omp_manual_direct",
		&transform_matrices_aos_to_aosoa, SCALE_HAMILT, &transform_matrix_aos_to_soa);

	// BENCHMARK: run optimised manually vectorised kernel
	benchmark(
		[&]() // lambda expression
		{
			commutator_omp_manual_direct_perm4to5(
				reinterpret_cast<real_vec_t*>(sigma_in),
				reinterpret_cast<real_vec_t*>(sigma_out),
				reinterpret_cast<real_t*>(hamiltonian),
				num, 0, 0.0, 0.0);
		},
		"commutator_omp_manual_direct_perm4to5",
		&transform_matrices_aos_to_aosoa, SCALE_HAMILT, &transform_matrix_aos_to_soa);


	delete hamiltonian;
	delete sigma_in;
	delete sigma_out;
	delete sigma_reference;

	return 0;
}

