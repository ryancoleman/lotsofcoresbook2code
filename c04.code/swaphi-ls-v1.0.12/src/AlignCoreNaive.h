/*
 * AlignCoreNaive.h
 *
 *  Created on: Feb 25, 2014
 *      Author: yongchao
 */

#ifndef ALIGNCORENAIVE_H_
#define ALIGNCORENAIVE_H_

#include "AlignCore.h"

class AlignCoreNaive: public AlignCore {
public:
	AlignCoreNaive(Align* align, int dev) :
			AlignCore(align, dev) {
	}
	~AlignCoreNaive() {
	}

	/*alignment*/
	void align();

private:
	/*scale and vectorize*/
	__attribute__((target(mic)))     int32_t _swKernel(
			const int8_t* __restrict__ seq1, const int32_t _s1Length,
			const int8_t* __restrict__ seq2, const int32_t _s2Length,
			const int32_t gapOE, const int32_t gapExtend, const int32_t match,
			const int32_t mismatch, const int32_t numThreads);
};

#endif /* ALIGNCORENAIVE_H_ */
