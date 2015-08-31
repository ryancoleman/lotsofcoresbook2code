// Copyright (c) 2015 Matthias Noack (ma.noack.pr@gmail.com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma OPENCL EXTENSION cl_khr_fp64 : enable

#ifndef DIM
	#define DIM 7
#endif
#ifndef VEC_LENGTH
	#define VEC_LENGTH 8
#endif

#ifdef SINGLE_PRECISION
	#define FLOATVEC_HELPER(n) float ## n
	#define FLOATVEC(n) FLOATVEC_HELPER(n)
	typedef FLOATVEC(VEC_LENGTH) real_vec_t;
	typedef FLOATVEC(2) real_2_t;
	typedef float real_t;
#else
	#define DOUBLEVEC_HELPER(n) double ## n
	#define DOUBLEVEC(n) DOUBLEVEC_HELPER(n)
	typedef DOUBLEVEC(VEC_LENGTH) real_vec_t;
	typedef DOUBLEVEC(2) real_2_t;
	typedef double real_t;
#endif

