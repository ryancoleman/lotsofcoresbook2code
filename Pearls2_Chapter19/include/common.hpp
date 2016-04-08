// Copyright (c) 2015 Matthias Noack (ma.noack.pr@gmail.com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef common_hpp
#define common_hpp

#include <complex>
#include <cstdint> // fixed width integers
#include <functional>
#include <iostream>

// SIMD vector libraries
#if defined(VEC_INTEL) || defined(VEC_VC) || defined(VEC_VCL)
	#include <Vc/vector.h>
	#include <vectorclass.h>
	#ifdef __MIC__
	#include <micvec.h>
	#else
	#include <dvec.h>
	#endif

	// compile time constant with defaults
	// assumed SIMD width (4 for Xeon, 8 for Xeon Phi)
	#ifndef VEC_LENGTH
		#ifdef SINGLE_PRECISION
			#define VEC_LENGTH VC_FLOAT_V_SIZE // use Vc
		#else
			#define VEC_LENGTH VC_DOUBLE_V_SIZE // use Vc
		#endif
	#endif

	// use one library to define real_vec_t
	#ifdef VEC_INTEL
//		#warning "Using Intel vector classes"
		#ifdef __MIC__
		typedef F64vec8 double_v;
		typedef F32vec16 float_v;
		#else
		typedef F64vec4 double_v;;
		typedef F32vec8 float_v;
		#endif
	#elif defined(VEC_VC)
//		#warning "Using Vc"
		using Vc::double_v;
		using Vc::float_v;
	#elif defined(VEC_VCL)
#undef USE_VCL_ORIGINAL
#if !defined(USE_VCL_ORIGINAL)
#ifdef __MIC__
	#warning "Using VCL (Vector Class Library) with Vec8dMod"
class Vec8dMod : public Vec8d {
	public:
		Vec8dMod(double d) : Vec8d(d) {}
		friend Vec8d operator *(const Vec8dMod &a, const double &b) { return Vec8d(_mm512_mul_pd(a, _mm512_set1_pd(b))); } 
};
#define Vec8d Vec8dMod
#else
class Vec4dMod : public Vec4d {
	public:
		Vec4dMod(double d) : Vec4d(d) {}
		friend Vec4d operator *(const Vec4dMod &a, const double &b) { return Vec4d(_mm256_mul_pd(a, _mm256_set1_pd(b))); } 
};
#define Vec4d Vec4dMod
#endif
#endif
//		#warning "Using VCL (Vector Class Library)"
		#ifdef __MIC__
		typedef Vec8d double_v;
		typedef Vec16f float_v;
		#else
		typedef Vec4d double_v;
		typedef Vec8f float_v;
		#endif
	#endif
#else
	// we do not use any SIMD library
//	#warning "No vector library defined"
	#ifndef VEC_LENGTH
		#ifdef __MIC__
			#define VEC_LENGTH 8
		#else
			#define VEC_LENGTH 4
		#endif
	#endif
	struct float_v { float data[VEC_LENGTH]; };
	struct double_v { double data[VEC_LENGTH]; };
#endif

#ifdef SINGLE_PRECISION
//	#warning "Using single precision (float)"
	using real_t = float;
	using real_vec_t = float_v;
#else
	using real_t = double;
	using real_vec_t = double_v;
#endif

using complex_t = std::complex<real_t>;

// alignment for memory allocations
#ifndef DEFAULT_ALIGNMENT
	#define DEFAULT_ALIGNMENT 64
#endif
// number of kernel iterations (including warmup)
#ifndef NUM_ITERATIONS
	#define NUM_ITERATIONS 26
#endif
// number of warmup iterations not taken into statistics
#ifndef NUM_WARMUP
	#define NUM_WARMUP 1
#endif	
// matrix dimension (based on actual application value)
#ifndef DIM
	#define DIM 7
#endif
// number of matrices in the sigma vectors (based on actual application value)
#ifndef NUM
	#define NUM 512*1024
#endif

#define STR_HELPER(x) #x
#define STR(x) STR_HELPER(x)

constexpr std::nullptr_t NO_TRANSFORM = nullptr;
constexpr bool NO_SCALE_HAMILT = false;
constexpr bool SCALE_HAMILT = true;

void print_compile_config(std::ostream& out);

// aligned_alloc from C++11 is not available for the Phi
template<typename T>
T* allocate_aligned(size_t size, size_t alignment = DEFAULT_ALIGNMENT)
{
	T* ptr = nullptr;
	int err = 0;
	if ((err = posix_memalign((void**)&ptr, alignment, size * sizeof(T))))
		std::cerr << "Error: posix_memalign() returned: " << err << std::endl;
	return ptr;
}

// intitialise sigma matrices
void initialise_sigma(complex_t* sigma_in, complex_t* sigma_out, size_t dim, size_t num);

// initialise hamiltonian
void initialise_hamiltonian(complex_t* hamiltonian, size_t dim);

// scale matrix by factor
//void transform_matrix_scale_aos(complex_t* matrix, size_t dim, complex_t factor);
void transform_matrix_scale_aos(complex_t* matrix, size_t dim, real_t factor);

// transform matrix format of a complex matrix from array of structs (AoS)
// RIRIRI... to struct of array (SoA) RRR...III...
// AoS: 
//     struct complex_t { real x, y; };
//     complex_t matrix[size];
// SoA: 
//     struct { real x[size], y[size]; } matrix;
void transform_matrix_aos_to_soa(complex_t* matrix, size_t dim);

// transform a vector of complex AoS matrices into an interleaved hybrid SoA
// format (AoSoA) with an inner size of the SIMD-width specified by VEC_LENGTH
// RIRIRI...RIRIRI... => RRR...III...RRR...III..
// AoS:
//     struct complex_t { real x, y; }; 
//     complex_t matrices[size * num];
// AoSoA: 
//     struct complex_t { real x[VEC_LENGTH], y[VEC_LENGTH]; }; 
//     complex_t matrix[size * num / VEC_LENGTH];
void transform_matrices_aos_to_aosoa(complex_t* matrices, size_t dim, size_t num, size_t vec_length = VEC_LENGTH);

// similar to transform_matrices_aos_to_aosoa, but packs real and imaginary
// parts of complex numbers differently:
// stores packages of interleaved matrices, with all real parts of the package
// preceding all the imaginare parts
void transform_matrices_aos_to_aosoa_gpu(complex_t* matrices, size_t dim, size_t num, size_t vec_length = VEC_LENGTH);

// returns the sum of the absolute values of the element-wise differences as
// measure of deviation
real_t compare_matrices(complex_t* a, complex_t* b, size_t dim, size_t num);

// benchmarks a kernel function and prints statistics to stdout
void benchmark_kernel(std::function<void()> kernel, std::string name, size_t overall_runs, size_t warmup_runs);

#endif // common_hpp

