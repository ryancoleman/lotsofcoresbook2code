// Copyright (c) 2015 Matthias Noack (ma.noack.pr@gmail.com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include "common.hpp"

#include <cstring> // memcpy
#include <cmath> // abs
#include <iostream>

#include "ham/util/time.hpp"

using namespace ham::util;

void print_compile_config(std::ostream& out)
{
	out << "DIM: " << DIM << std::endl;
	out << "NUM: " << NUM << std::endl;
	out << "NUM_ITERATIONS: " << NUM_ITERATIONS << std::endl;
	out << "NUM_WARMUP: " << NUM_WARMUP << std::endl;
	out << "PRECISION: ";
	#ifdef SINGLE_PRECISION
		out << "SINGLE";
	#else
		out << "DOUBLE";
	#endif
	out << std::endl;
	out << "DEFAULT_ALIGNMENT: " << DEFAULT_ALIGNMENT << std::endl;
	out << "VEC_LIB: ";
	#if defined(VEC_INTEL)
		out << "VEC_INTEL";
	#elif defined(VEC_VC)
		out << "VEC_VC";
	#elif defined(VEC_VCL)
		out << "VEC_VCL";
	#else
		out << "NO_VEC_LIB";
	#endif
	out << std::endl;
        #if defined(USE_INITZERO)
             	out << "USE_INITZERO" << std::endl;
        #endif
	out << "VEC_LENGTH: " << VEC_LENGTH << std::endl;
}

void initialise_sigma(complex_t* sigma_in, complex_t* sigma_out, size_t dim, size_t num)
{
	size_t size_sigma = dim * dim;
	#pragma omp parallel for
	for (size_t sigma_id = 0; sigma_id < num; ++sigma_id)
		for (size_t i = 0; i < size_sigma; ++i)
		{
			real_t x = static_cast<real_t>(sigma_id) / num;
			real_t y = static_cast<real_t>(i) / size_sigma;
			sigma_in[sigma_id * size_sigma + i] = complex_t(x - y, y - x);
			sigma_out[sigma_id * size_sigma + i] = complex_t(0.0, 0.0);
		}
}

void initialise_hamiltonian(complex_t* hamiltonian, size_t dim)
{
	size_t size = dim * dim;
	for (size_t i = 0; i < size; ++i)
	{
		hamiltonian[i] = 1.0 - static_cast<real_t>(i) / size;
	}
}

void transform_matrix_scale_aos(complex_t* matrix, size_t dim, real_t factor)
{
	for (size_t i = 0; i < dim * dim; ++i)
		matrix[i] *= factor;
}

void transform_matrix_aos_to_soa(complex_t* matrix, size_t dim)
{
	size_t size = dim * dim;

	// create a temporary copy of matrix
	complex_t* matrix_tmp = new complex_t[size];
	std::memcpy(matrix_tmp, matrix, sizeof(complex_t) * size);
	
	// copy back with new layout
	real_t* matrix_r = reinterpret_cast<real_t*>(matrix);
	#pragma omp parallel for
	for (size_t i = 0; i < dim; ++i)
		for (size_t j = 0; j < dim; ++j)
		{
			// NOTE: matrix_r is real, matrix_tmp is complex (i.e. a tuple of 2 reals)
			matrix_r[i * dim + j] = matrix_tmp[i * dim + j].real();
			matrix_r[size + i * dim + j] = matrix_tmp[i * dim + j].imag();
		}

	delete [] matrix_tmp;
}

void transform_matrices_aos_to_aosoa(complex_t* matrices, size_t dim, size_t num, size_t vec_length)
{
	// indexing macros for the transformed data layout
	//#define package_id ((n / vec_length) * vec_length * 2 * dim * dim)
	auto package_id = [&](size_t m) { return (m / vec_length) * vec_length * 2 * dim * dim; };
	//#define sigma_id (n % vec_length)
	auto sigma_id = [&](size_t m) { return m % vec_length; };
	//#define sigma_real(i, j) (package_id + 2 * vec_length * (dim * (i) + (j)) + sigma_id)
	auto sigma_real = [&](size_t i, size_t j, size_t m) { return package_id(m) + 2 * vec_length * (dim * i + j) + sigma_id(m); };
	//#define sigma_imag(i, j) (package_id + 2 * vec_length * (dim * (i) + (j)) + vec_length + sigma_id)
	auto sigma_imag = [&](size_t i, size_t j, size_t m) { return package_id(m) + 2 * vec_length * (dim * i + j) + vec_length + sigma_id(m); };

	size_t size = num * dim * dim;

	// create a temporary copy of matrix
	complex_t* matrices_tmp = new complex_t[size];
	std::memcpy(matrices_tmp, matrices, sizeof(complex_t) * size);

	// copy back with new layout
	real_t* matrices_r = reinterpret_cast<real_t*>(matrices);
	#pragma omp parallel for
	for (size_t m = 0; m < num; ++m)
	{
		for (size_t i = 0; i < dim; ++i)
		{
			for (size_t j = 0; j < dim; ++j)
			{
				matrices_r[sigma_real(i, j, m)] = matrices_tmp[m * dim * dim + i * dim + j].real();
				matrices_r[sigma_imag(i, j, m)] = matrices_tmp[m * dim * dim + i * dim + j].imag();				
			}
		}
	}

	delete [] matrices_tmp;
}

void transform_matrices_aos_to_aosoa_gpu(complex_t* matrices, size_t dim, size_t num, size_t vec_length)
{
	// lambdas for indexing
	auto package_id = [&](size_t m) { return (m / vec_length) * vec_length * 2 * dim * dim; };
	auto sigma_id = [&](size_t m) { return m % vec_length; };
	auto sigma_real = [&](size_t i, size_t j, size_t m) { return package_id(m) + vec_length * (dim * i + j) + sigma_id(m); };
	auto sigma_imag = [&](size_t i, size_t j, size_t m) { return package_id(m) + dim * dim * vec_length + vec_length * (dim * i + j) + sigma_id(m); };

	size_t size = num * dim * dim;

	// create a temporary copy of matrix
	complex_t* matrices_tmp = new complex_t[size];
	std::memcpy(matrices_tmp, matrices, sizeof(complex_t) * size);

	// copy back with new layout
	real_t* matrices_r = reinterpret_cast<real_t*>(matrices);
	#pragma omp parallel for
	for (size_t m = 0; m < num; ++m)
	{
		for (size_t i = 0; i < dim; ++i)
		{
			for (size_t j = 0; j < dim; ++j)
			{
				matrices_r[sigma_real(i, j, m)] = matrices_tmp[m * dim * dim + i * dim + j].real();
				matrices_r[sigma_imag(i, j, m)] = matrices_tmp[m * dim * dim + i * dim + j].imag();				
			}
		}
	}

	delete [] matrices_tmp;
}

real_t compare_matrices(complex_t* a, complex_t* b, size_t dim, size_t num)
{
	real_t deviation = 0.0;
	size_t size = num * dim * dim;
	
	// add up absolute values of the element-wise differences
	for (size_t i = 0; i < size; ++i)
	{
		deviation += abs(a[i] - b[i]);
	}

	return deviation;
}

void benchmark_kernel(std::function<void()> kernel, std::string name, size_t overall_runs, size_t warmup_runs)
{
	time::statistics stats(overall_runs, warmup_runs);
	for (size_t i = 0; i < overall_runs; ++i)
	{
		time::timer t;
		kernel();
		stats.add(t);
	}
	std::cout << name << "\t" << stats.string() << std::endl;
}

