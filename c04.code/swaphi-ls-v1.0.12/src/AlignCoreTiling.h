/*
 * AlignCoreTiling.h
 *
 *  Created on: Feb 25, 2014
 *      Author: yongchao
 */

#ifndef ALIGNCORETILING_H_
#define ALIGNCORETILING_H_

#include "AlignCore.h"

class AlignCoreTiling: public AlignCore {
public:
	AlignCoreTiling(Align* align, int dev) :
			AlignCore(align, dev) {
	}
	~AlignCoreTiling() {
	}

	/*alignment*/
	void align();

private:
	/*distributed-memory version*/
	void _alignDistributed(const int32_t vMpiNumBlocks,
			const int32_t hMpiNumBlocks, const int32_t realNumProcs);
	void _shiftHFBuffers(const int32_t* mpiHF, const int32_t* mpiHFC,
			const int32_t hfMpiVecSize, const int32_t realNumProcs);
	__attribute__((target(mic)))       int32_t _swKernelDistributed(
			const int8_t* __restrict__ seq1, const int32_t s1Length,
			const int8_t* __restrict__ seq2Rev, const int32_t s2Length,
			int2* __restrict__ mpiHE, const int32_t* __restrict__ mpiH,
			int32_t* __restrict__ mpiHC, const int32_t* __restrict__ mpiF,
			int32_t* __restrict__ mpiFC, const int32_t mpiHD,
			int32_t* __restrict__ mpiHDC, const int32_t gapOE,
			const int32_t gapExtend, const int32_t match,
			const int32_t mismatch, const int32_t numThreads,
			const int32_t alpha, const int32_t beta);

	/*shared-memory version*/
	void _alignShared();
	__attribute__((target(mic)))       int32_t _swKernelShared(
			const int8_t* __restrict__ seq1, const int32_t s1Length,
			const int8_t* __restrict__ seq2Rev, const int32_t s2Length,
			const int32_t gapOE, const int32_t gapExtend, const int32_t match,
			const int32_t mismatch, const int32_t numThreads,
			const int32_t alpha, const int32_t beta);

};

#endif /* ALIGNCORETILING_H_ */
