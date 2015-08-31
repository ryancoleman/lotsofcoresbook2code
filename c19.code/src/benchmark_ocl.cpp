// Copyright (c) 2015 Matthias Noack (ma.noack.pr@gmail.com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <iostream>

#include <cstring> // memcpy
#include <cmath>
#include <cstddef>
#include <string>

#include "clu_runtime/clu.h"
#include "ham/util/time.hpp" // ham::util::time

#include "common.hpp"
#include "kernel/kernel.hpp"

using namespace ham::util;

#ifndef DEVICE_TYPE
	// CL_DEVICE_TYPE_CPU  
	// CL_DEVICE_TYPE_ACCELERATOR
	// CL_DEVICE_TYPE_GPU
	#define DEVICE_TYPE CL_DEVICE_TYPE_ACCELERATOR
#endif
#ifndef INTEL_PREFETCH_LEVEL
	#define INTEL_PREFETCH_LEVEL 1
#endif
#ifndef PACKAGES_PER_WG
	#define PACKAGES_PER_WG 4
#endif
#ifndef NUM_SUB_GROUPS
	#define NUM_SUB_GROUPS 2
#endif
#ifndef CHUNK_SIZE
	#define CHUNK_SIZE 16
#endif
#ifndef WARP_SIZE
	#define WARP_SIZE 32
#endif

bool ocl_error_handler(cl_int err, const std::string& function_name, bool terminate = true)
{
	bool error = err != CL_SUCCESS;
	if (error)
	{
		std::cerr << "Error: "<< function_name << " returned " << cluPrintError(err) << std::endl;
		if (terminate)
			exit(-1);
	}
	return error;
}

void benchmark_ocl_kernel(cl_kernel kernel, std::string name, clu_nd_range range, size_t num, size_t overall_runs, size_t warmup_runs)
{
	cl_int err = 0;
	cl_event event;
	cl_ulong t_start = 0;
	cl_ulong t_end = 0;
	clu_enqueue_params params = { range, CLU_DEFAULT_Q, 0, nullptr, &event };

	// benchmark loop
	time::statistics stats(overall_runs, warmup_runs);
	for (size_t i = 0; i < overall_runs; ++i)
	{
		// execute kernel
		cluEnqueue(kernel, &params);
		ocl_error_handler(err, "cluEnqueue()");
		err = clWaitForEvents(1, &event);
		ocl_error_handler(err, "clWaitForEvents()");
		
		// get times for statistics
		err = clGetEventProfilingInfo(event, CL_PROFILING_COMMAND_START, sizeof(cl_ulong), &t_start, nullptr);
		ocl_error_handler(err, "clGetEventProfilingInfo()");
		err = clGetEventProfilingInfo(event, CL_PROFILING_COMMAND_END, sizeof(cl_ulong), &t_end, nullptr);
		ocl_error_handler(err, "clGetEventProfilingInfo()");
		stats.add(static_cast<time::rep>(t_end - t_start));
	}

	std::cout << name << "\t" << stats.string() << std::endl;
}

int main(void)
{
	print_compile_config(std::cerr);
	std::cerr << "VEC_LENGTH_AUTO: " << VEC_LENGTH_AUTO << std::endl;
	std::cerr << "DEVICE_TYPE: " << DEVICE_TYPE << std::endl;

	// constants
	const size_t dim = DIM;
	const size_t num = NUM;
	const real_t hbar = 1.0 / std::acos(-1.0); // == 1 / Pi
	const real_t dt = 1.0e-3;
	const real_t hdt = dt / hbar;

	real_t deviation = 0.0;

	// allocate memory
	size_t size_hamiltonian = dim * dim;
	size_t size_hamiltonian_byte = sizeof(complex_t) * size_hamiltonian;
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

	// setup OpenCL using CLU
#if (DEVICE_TYPE == CL_DEVICE_TYPE_CPU)
	const std::string compile_options_impl = " -cl-mad-enable -auto-prefetch-level=" STR(INTEL_PREFETCH_LEVEL) " "; 
#elif (DEVICE_TYPE == CL_DEVICE_TYPE_ACCELERATOR)
	const std::string compile_options_impl = " -cl-mad-enable -auto-prefetch-level=" STR(INTEL_PREFETCH_LEVEL) " "; // -cl-finite-math-only -cl-no-signed-zeros "; 
#elif (DEVICE_TYPE == CL_DEVICE_TYPE_GPU)
	const std::string compile_options_impl = ""; // -cl-nv-verbose -cl-nv-opt-level=3 -cl-mad-enable -cl-strict-aliasing -cl-nv-arch sm_35 -cl-nv-maxrregcount=64 "; 
#endif 
	const std::string compile_options_common = "-Iinclude -DNUM=" STR(NUM) " -DDIM=" STR(DIM) + compile_options_impl;
	const std::string compile_options_auto = compile_options_common + " -DVEC_LENGTH=" STR(VEC_LENGTH_AUTO) " -DPACKAGES_PER_WG=" STR(PACKAGES_PER_WG);
	const std::string compile_options_manual = compile_options_common + " -DVEC_LENGTH=" STR(VEC_LENGTH) " -DPACKAGES_PER_WG=" STR(PACKAGES_PER_WG);
	const std::string compile_options_gpu = compile_options_common + " -DVEC_LENGTH=2 -DCHUNK_SIZE=" STR(CHUNK_SIZE) " -DNUM_SUB_GROUPS=" STR(NUM_SUB_GROUPS);

	cl_int err = 0;
	clu_initialize_params init_params = {}; // NOTE: null initialise!
	init_params.default_queue_props = CL_QUEUE_PROFILING_ENABLE;
	init_params.preferred_device_type = DEVICE_TYPE; // set device type
	init_params.compile_options = compile_options_common.c_str(); // set default compile options
	cl_int status = cluInitialize(&init_params);

	// output the used device
	cl_device_id dev_id;
	err = clGetCommandQueueInfo(CLU_DEFAULT_Q, CL_QUEUE_DEVICE, sizeof(cl_device_id), &dev_id, nullptr);
	ocl_error_handler(err, "clGetCommandQueueInfo()");
	clu_device_info dev_info = cluGetDeviceInfo(dev_id, &err);
	ocl_error_handler(err, "cluGetDeviceInfo()");
	std::cerr << "Using device: " << dev_info.device_name << std::endl;

	// allocate OpenCL device memory
	cl_mem hamiltonian_ocl = clCreateBuffer(CLU_CONTEXT, CL_MEM_READ_ONLY, size_hamiltonian_byte, 0, &err);
	ocl_error_handler(err, "clCreateBuffer(hamiltonian_ocl)");
	cl_mem sigma_in_ocl = clCreateBuffer(CLU_CONTEXT, CL_MEM_READ_WRITE, size_sigma_byte, 0, &err);
	ocl_error_handler(err, "clCreateBuffer(sigma_in_ocl)");
	cl_mem sigma_out_ocl = clCreateBuffer(CLU_CONTEXT, CL_MEM_READ_WRITE, size_sigma_byte, 0, &err);
	ocl_error_handler(err, "clCreateBuffer(sigma_out_ocl)");

	// function to build and set-up a kernel
	auto prepare_kernel = [&](const std::string& file_name, const std::string& kernel_name, const std::string& compile_options)
	{
		// build kernel
		cl_program prog = cluBuildSourceFromFile(file_name.c_str(), compile_options.c_str(), &err);
		if (ocl_error_handler(err, "cluBuildSourceFromFile()", false))
		{
			std::cerr << "OpenCL build log for: " << file_name << std::endl
			          << cluGetBuildErrors(prog) << std::endl
			          << "--- end of build log ---" << std::endl;
			exit(-1);
		}
		cl_kernel kernel = clCreateKernel(prog, kernel_name.c_str(), &err);
		ocl_error_handler(err, "clCreateKernel()");

		// set kernel arguments
		err = clSetKernelArg(kernel, 0, sizeof(cl_mem), static_cast<const void*>(&sigma_in_ocl));
		ocl_error_handler(err, "clSetKernelArg(0)");
		err = clSetKernelArg(kernel, 1, sizeof(cl_mem), static_cast<const void*>(&sigma_out_ocl));
		ocl_error_handler(err, "clSetKernelArg(1)");
		err = clSetKernelArg(kernel, 2, sizeof(cl_mem), static_cast<const void*>(&hamiltonian_ocl));
		ocl_error_handler(err, "clSetKernelArg(2)");
		// GCC bug work-around (#ifdef __GNUC__) & type conversion
		// http://stackoverflow.com/questions/19616610/c11-lambda-doesnt-take-const-variable-by-reference-why
		int32_t num_tmp = static_cast<int32_t>(num);
		err = clSetKernelArg(kernel, 3, sizeof(int32_t), static_cast<const void*>(&num_tmp));
		ocl_error_handler(err, "clSetKernelArg(3)");
		int32_t dim_tmp = static_cast<int32_t>(dim);
		err = clSetKernelArg(kernel, 4, sizeof(int32_t), static_cast<const void*>(&dim_tmp));
		ocl_error_handler(err, "clSetKernelArg(4)");
		err = clSetKernelArg(kernel, 5, sizeof(real_t), static_cast<const void*>(&hbar));
		ocl_error_handler(err, "clSetKernelArg(5)");
		err = clSetKernelArg(kernel, 6, sizeof(real_t), static_cast<const void*>(&dt));
		ocl_error_handler(err, "clSetKernelArg(6)");

		return kernel;
	}; // prepare_kernel

	auto write_hamiltonian = [&]()
	{
		err = clEnqueueWriteBuffer(CLU_DEFAULT_Q, hamiltonian_ocl, CL_TRUE, 0, size_hamiltonian_byte, hamiltonian, 0, nullptr, nullptr);
		ocl_error_handler(err, "clEnqueueWriteBuffer(hamiltonian)");
	}; // write_hamiltonian

	// lambda to write data to the device
	auto write_sigma = [&]()
	{
		// write data to device
		err = clEnqueueWriteBuffer(CLU_DEFAULT_Q, sigma_in_ocl, CL_TRUE, 0, size_sigma_byte, sigma_in, 0, nullptr, nullptr);
		ocl_error_handler(err, "clEnqueueWriteBuffer(sigma_in_ocl)");
		err = clEnqueueWriteBuffer(CLU_DEFAULT_Q, sigma_out_ocl, CL_TRUE, 0, size_sigma_byte, sigma_out, 0, nullptr, nullptr);
		ocl_error_handler(err, "clEnqueueWriteBuffer(sigma_out_ocl)");
	}; // write_sigma

	// lambda to get the result from the device and compare it with the reference
	auto read_and_compare_sigma = [&]()
	{
		// read data from device	
		err = clEnqueueReadBuffer(CLU_DEFAULT_Q, sigma_out_ocl, CL_TRUE, 0, size_sigma_byte, sigma_out, 0, nullptr, nullptr);
		ocl_error_handler(err, "clEnqueueReadBuffer(sigma_out_ocl)");
		// compute deviation from reference	(small deviations are expected)
		deviation = compare_matrices(sigma_out, sigma_reference_transformed, dim, num);
		std::cerr << "Deviation:\t" << deviation << std::endl;
	}; // read_and_compare_sigma

	// Lambda to: transform memory, benchmark, compare results
	// NOTE:
	// typedef struct
	// {
	//     cl_uint   dim; 
	//     size_t    global[3];
	//     size_t    local[3];
	//     size_t    offset[3];
	// } clu_nd_range;
	auto benchmark = [&](const std::string& file_name, const std::string& kernel_name,
	                     const std::string& compile_options, size_t vec_length, clu_nd_range range,
	                     decltype(&transform_matrices_aos_to_aosoa) transformation_sigma,
	                     bool scale_hamiltonian,
	                     decltype(&transform_matrix_aos_to_soa) transformation_hamiltonian)
	{
		initialise_hamiltonian(hamiltonian, dim);
		if (scale_hamiltonian) 
			transform_matrix_scale_aos(hamiltonian, dim, dt / hbar); // pre-scale hamiltonian
		if (transformation_hamiltonian)
			transformation_hamiltonian(hamiltonian, dim);	
		write_hamiltonian();
	
		initialise_sigma(sigma_in, sigma_out, dim, num);
		std::memcpy(sigma_reference_transformed, sigma_reference, size_sigma_byte);
		// transform memory layout if a transformation is specified
		if (transformation_sigma)
		{
			// transform reference for comparison
			transformation_sigma(sigma_reference_transformed, dim, num, vec_length);
			// tranform sigma
			transformation_sigma(sigma_in, dim, num, vec_length);
		}
		write_sigma();

		cl_kernel kernel = prepare_kernel(file_name, kernel_name, compile_options);
		benchmark_ocl_kernel(kernel, kernel_name, range, num, NUM_ITERATIONS, NUM_WARMUP);
		
		read_and_compare_sigma();
	}; // benchmark

	// BENCHMARK: initial kernel
	benchmark("src/kernel/commutator_ocl_initial.cl", "commutator_ocl_initial",
	          compile_options_common, VEC_LENGTH,
	          { 1, // NDRange dimension
	            { num }, // global size
	            { }, // local size
	            { } // offset
	          }, NO_TRANSFORM, NO_SCALE_HAMILT, NO_TRANSFORM);

	// BENCHMARK: refactored initial kernel
	benchmark("src/kernel/commutator_ocl_refactored.cl", "commutator_ocl_refactored",
	          compile_options_auto, VEC_LENGTH,
	          { 1, // NDRange dimension
	            { num }, // global size
	            { }, // local size
	            { } // offset
	          }, NO_TRANSFORM, NO_SCALE_HAMILT, NO_TRANSFORM);

	// BENCHMARK: refactored initial kernel with direct store
	benchmark("src/kernel/commutator_ocl_refactored_direct.cl", "commutator_ocl_refactored_direct",
	          compile_options_auto, VEC_LENGTH,
	          { 1, // NDRange dimension
	            { num }, // global size
	            { }, // local size
	            { } // offset
	          }, NO_TRANSFORM, SCALE_HAMILT, NO_TRANSFORM);

	// BENCHMARK: automatically vectorised kernel with naive NDRange and indexing
	benchmark("src/kernel/commutator_ocl_aosoa_naive.cl", "commutator_ocl_aosoa_naive",
	          compile_options_auto, VEC_LENGTH_AUTO,
	          { 1, // NDRange dimension
	            { num }, // global size
	            { }, // local size
	            { } // offset
	          }, &transform_matrices_aos_to_aosoa, SCALE_HAMILT, &transform_matrix_aos_to_soa);

	
	// BENCHMARK: automatically vectorised kernel with naive NDRange and indexing and compile time constants
	benchmark("src/kernel/commutator_ocl_aosoa_naive_constants.cl", "commutator_ocl_aosoa_naive_constants",
	          compile_options_auto, VEC_LENGTH_AUTO,
	          { 1, // NDRange dimension
	            { num }, // global size
	            { }, // local size
	            { } // offset
	          }, &transform_matrices_aos_to_aosoa, SCALE_HAMILT, &transform_matrix_aos_to_soa);

	// BENCHMARK: automatically vectorised kernel with naive NDRange and indexing and direct store
	benchmark("src/kernel/commutator_ocl_aosoa_naive_direct.cl", "commutator_ocl_aosoa_naive_direct",
	          compile_options_auto, VEC_LENGTH_AUTO,
	          { 1, // NDRange dimension
	            { num }, // global size
	            { }, // local size
	            { } // offset
	          }, &transform_matrices_aos_to_aosoa, SCALE_HAMILT, &transform_matrix_aos_to_soa);
	
	// BENCHMARK: automatically vectorised kernel with naive NDRange and indexing, compile time constants, and direct store
	benchmark("src/kernel/commutator_ocl_aosoa_naive_constants_direct.cl", "commutator_ocl_aosoa_naive_constants_direct",
	          compile_options_auto, VEC_LENGTH_AUTO,
	          { 1, // NDRange dimension
	            { num }, // global size
	            { }, // local size
	            { } // offset
	          }, &transform_matrices_aos_to_aosoa, SCALE_HAMILT, &transform_matrix_aos_to_soa);
	
	// BENCHMARK: automatically vectorised kernel with compiler-friendly NDRange and indexing 
	benchmark("src/kernel/commutator_ocl_aosoa.cl", "commutator_ocl_aosoa",
	          compile_options_auto, VEC_LENGTH_AUTO,
	          { 2, // NDRange dimension
	            { VEC_LENGTH_AUTO, num / (VEC_LENGTH_AUTO) }, // global size
	            { (VEC_LENGTH_AUTO), PACKAGES_PER_WG }, // local size
	            { } // offset
	          }, &transform_matrices_aos_to_aosoa, SCALE_HAMILT, &transform_matrix_aos_to_soa);
	
	// BENCHMARK: automatically vectorised kernel with compiler-friendly NDRange and indexing, and compile time constants
	benchmark("src/kernel/commutator_ocl_aosoa_constants.cl", "commutator_ocl_aosoa_constants",
	          compile_options_auto, VEC_LENGTH_AUTO,
	          { 2, // NDRange dimension
	            { VEC_LENGTH_AUTO, num / (VEC_LENGTH_AUTO) }, // global size
	            { (VEC_LENGTH_AUTO), PACKAGES_PER_WG }, // local size
	            { } // offset
	          }, &transform_matrices_aos_to_aosoa, SCALE_HAMILT, &transform_matrix_aos_to_soa);

	// BENCHMARK: automatically vectorised kernel with compiler-friendly NDRange and indexing, and direct store
	benchmark("src/kernel/commutator_ocl_aosoa_direct.cl", "commutator_ocl_aosoa_direct",
	          compile_options_auto, VEC_LENGTH_AUTO,
	          { 2, // NDRange dimension
	            { VEC_LENGTH_AUTO, num / (VEC_LENGTH_AUTO) }, // global size
	            { (VEC_LENGTH_AUTO), PACKAGES_PER_WG }, // local size
	            { } // offset
	          }, &transform_matrices_aos_to_aosoa, SCALE_HAMILT, &transform_matrix_aos_to_soa);
	
	// BENCHMARK: automatically vectorised kernel with compiler-friendly NDRange and indexing, compile time constants, and direct store
	benchmark("src/kernel/commutator_ocl_aosoa_constants_direct.cl", "commutator_ocl_aosoa_constants_direct",
	          compile_options_auto, VEC_LENGTH_AUTO,
	          { 2, // NDRange dimension
	            { VEC_LENGTH_AUTO, num / (VEC_LENGTH_AUTO) }, // global size
	            { (VEC_LENGTH_AUTO), PACKAGES_PER_WG }, // local size
	            { } // offset
	          }, &transform_matrices_aos_to_aosoa, SCALE_HAMILT, &transform_matrix_aos_to_soa);
	
	// BENCHMARK: automatically vectorised kernel with compiler-friendly NDRange and indexing, compile time constants, direct store, and permuted loops with temporaries
	benchmark("src/kernel/commutator_ocl_aosoa_constants_direct_perm.cl", "commutator_ocl_aosoa_constants_direct_perm",
	          compile_options_auto, VEC_LENGTH_AUTO,
	          { 2, // NDRange dimension
	            { VEC_LENGTH_AUTO, num / (VEC_LENGTH_AUTO) }, // global size
	            { (VEC_LENGTH_AUTO), PACKAGES_PER_WG }, // local size
	            { } // offset
	          }, &transform_matrices_aos_to_aosoa, SCALE_HAMILT, &transform_matrix_aos_to_soa);

	// BENCHMARK: manually vectorised kernel
	benchmark("src/kernel/commutator_ocl_manual_aosoa.cl", "commutator_ocl_manual_aosoa",
	          compile_options_manual, VEC_LENGTH,
	          { 1, // NDRange dimension
	            { num / VEC_LENGTH}, // global size
	            { }, // local size
	            { } // offset
	          }, &transform_matrices_aos_to_aosoa, SCALE_HAMILT, &transform_matrix_aos_to_soa);

	// BENCHMARK: manually vectorised kernel with compile time constants
	benchmark("src/kernel/commutator_ocl_manual_aosoa_constants.cl", "commutator_ocl_manual_aosoa_constants",
	          compile_options_manual, VEC_LENGTH,
	          { 1, // NDRange dimension
	            { num / VEC_LENGTH}, // global size
	            { }, // local size
	            { } // offset
	          }, &transform_matrices_aos_to_aosoa, SCALE_HAMILT, &transform_matrix_aos_to_soa);

	// BENCHMARK: manually vectorised kernel with compile time constants
	benchmark("src/kernel/commutator_ocl_manual_aosoa_constants_prefetch.cl", "commutator_ocl_manual_aosoa_constants_prefetch",
	          compile_options_manual, VEC_LENGTH,
	          { 1, // NDRange dimension
	            { num / VEC_LENGTH}, // global size
	            { }, // local size
	            { } // offset
	          }, &transform_matrices_aos_to_aosoa, SCALE_HAMILT, &transform_matrix_aos_to_soa);

	
	// BENCHMARK: manually vectorised kernel with direct store
	benchmark("src/kernel/commutator_ocl_manual_aosoa_direct.cl", "commutator_ocl_manual_aosoa_direct",
	          compile_options_manual, VEC_LENGTH,
	          { 1, // NDRange dimension
	            { num / VEC_LENGTH}, // global size
	            { }, // local size
	            { } // offset
	          }, &transform_matrices_aos_to_aosoa, SCALE_HAMILT, &transform_matrix_aos_to_soa);

	// BENCHMARK: manually vectorised kernel with compile time constants and direct store
	benchmark("src/kernel/commutator_ocl_manual_aosoa_constants_direct.cl", "commutator_ocl_manual_aosoa_constants_direct",
	          compile_options_manual, VEC_LENGTH,
	          { 1, // NDRange dimension
	            { num / VEC_LENGTH}, // global size
	            { }, // local size
	            { } // offset
	          }, &transform_matrices_aos_to_aosoa, SCALE_HAMILT, &transform_matrix_aos_to_soa);

	// BENCHMARK: manually vectorised kernel with compile time constants and direct store
	benchmark("src/kernel/commutator_ocl_manual_aosoa_constants_direct_prefetch.cl", "commutator_ocl_manual_aosoa_constants_direct_prefetch",
	          compile_options_manual, VEC_LENGTH,
	          { 1, // NDRange dimension
	            { num / VEC_LENGTH}, // global size
	            { }, // local size
	            { } // offset
	          }, &transform_matrices_aos_to_aosoa, SCALE_HAMILT, &transform_matrix_aos_to_soa);

	
	// BENCHMARK: manually vectorised kernel with compile time constants, direct store, and permuted loops with temporaries
	benchmark("src/kernel/commutator_ocl_manual_aosoa_constants_direct_perm.cl", "commutator_ocl_manual_aosoa_constants_direct_perm",
	          compile_options_manual, VEC_LENGTH,
	          { 1, // NDRange dimension
	            { num / VEC_LENGTH}, // global size
	            { }, // local size
	            { } // offset
	          }, &transform_matrices_aos_to_aosoa, SCALE_HAMILT, &transform_matrix_aos_to_soa);
	
	// BENCHMARK: final GPGPU kernel, optimised for Nvidia K40
	{ // keep things local
	auto ceil_n = [](size_t x, size_t n) { return ((x + n - 1) / n) * n; };
	size_t block_dim_x = ceil_n(dim * dim, WARP_SIZE);
	size_t block_dim_y = NUM_SUB_GROUPS;
	
	benchmark("src/kernel/commutator_ocl_gpu_final.cl", "commutator_ocl_gpu_final", compile_options_gpu,
	          2, // NOTE: vec_length has a fix value of 2 for this kernel
	          { 2, // NDRange dimension
	            { (NUM / (block_dim_y * CHUNK_SIZE)) * block_dim_x, block_dim_y }, // global size
	            { block_dim_x, block_dim_y }, // local size
	            { } // offset
	          }, NO_TRANSFORM, SCALE_HAMILT, NO_TRANSFORM);
	}


	// de-init CLU
	cluRelease(); 

	delete hamiltonian;
	delete sigma_in;
	delete sigma_out;
	delete sigma_reference;

	return 0;
}

