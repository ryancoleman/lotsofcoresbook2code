/*
 * AlignCoreNaive.cpp
 *
 *  Created on: Feb 25, 2014
 *      Author: yongchao
 */

#include "AlignCoreNaive.h"

void AlignCoreNaive::align() {
	int32_t maxScore;
	double runtime = 0;
	int8_t* seq1 = _seq1;
	int8_t* seq2Rev = _seq2;/*NOTE: Here seq2 has been reversed*/

	Utils::log("Number of device threads: %d\n", _numMicThreads);

	/*get the start time*/
	runtime = omp_get_wtime();

	/*start to offload*/
#pragma offload target(mic:_micIndex) \
		in(seq1: length(_seq1Length) MIC_ALLOC_FREE) \
		in(seq2Rev: length(_seq2Length) MIC_ALLOC_FREE)
	{

		maxScore = _swKernel(seq1, _seq1Length, seq2Rev, _seq2Length, _gapOE,
				_gapExtend, _match, _mismatch, _numMicThreads);

	}

	/*calculate the CPU time on the device*/
	runtime = omp_get_wtime() - runtime;

	if (maxScore == -1) {
		Utils::exit("Errors have occurred in the offload region\n");
	}
	Utils::log("Computation time on the Xeon Phi: %f seconds\n", runtime);

	/*compute GCUPS*/
	double product = _seq1Length;
	product *= _seq2Length;
	Utils::log("Performance in GCUPS: %f GCUPS\n",
			product / runtime / 1000000 / 1000);

	Utils::log("Optimal local alignment score: %d\n", maxScore);
}
/*scale and vectorize*/
__attribute__((target(mic)))                    int32_t AlignCoreNaive::_swKernel(
		const int8_t* __restrict__ seq1, const int32_t s1Length,
		const int8_t* __restrict__ seq2Rev, const int32_t s2Length,
		const int32_t gapOE, const int32_t gapExtend, const int32_t match,
		const int32_t mismatch, const int32_t numThreads) {
	int32_t maxScore;
#ifdef __MIC__
	__m512i *vecMaxScores;
	const register int32_t vecSize = s1Length + 2;
	int32_t* __restrict__ addrHA;
	int32_t* __restrict__ addrHB;
	int32_t* __restrict__ addrHC;
	int32_t* __restrict__ addrF;
	int32_t* __restrict__ addrEO;
	int32_t* __restrict__ addrEC;

	/*allocate space*/
	addrHA = (int32_t*)_mm_malloc(vecSize * sizeof(int32_t), 64);
	addrHB = (int32_t*)_mm_malloc(vecSize * sizeof(int32_t), 64);
	addrHC = (int32_t*)_mm_malloc(vecSize * sizeof(int32_t), 64);
	addrF = (int32_t*)_mm_malloc(vecSize * sizeof(int32_t), 64);
	addrEO = (int32_t*)_mm_malloc(vecSize * sizeof(int32_t), 64);
	addrEC = (int32_t*)_mm_malloc(vecSize * sizeof(int32_t), 64);
	if(!addrHA || !addrHB || !addrHC || !addrF || !addrEO || !addrEC) {
		printf("Memory allocation failed\n");
		return -2;
	}

	/*allocate maximum score buffer and initialize them*/
	vecMaxScores = (__m512i*)_mm_malloc(sizeof(__m512i) * numThreads, 64);
	for(int32_t i = 0; i < numThreads; ++i) {
		vecMaxScores[i] = _mm512_setzero_epi32();
	}

	/*initialize invariant values*/
	const register __m512i vecGapOE = _mm512_set1_epi32(gapOE);
	const register __m512i vecGapExtend = _mm512_set1_epi32(gapExtend);
	const register __m512i vecMatch = _mm512_set1_epi32(match);
	const register __m512i vecMismatch = _mm512_set1_epi32(mismatch);
	const register __m512i vecIndices = _mm512_set_epi32(15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0); /*from high to low*/
	const register __m512i vecZero = _mm512_setzero_epi32();
	const register __m512i vecPermuteIndex = _mm512_set_epi32(14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0, 15); /*from high to low*/

	/*initialize the intermediate buffers. use separate read/write buffers to avoid false sharing*/
	for(int32_t i = 0; i < vecSize; ++i) {
		addrHA[i] = 0;
		addrHB[i] = 0;
		addrHC[i] = 0;
		addrF[i] = 0;
		addrEO[i] = 0;
		addrEC[i] = 0;
	}

	/*set the number of threads*/
	omp_set_num_threads(numThreads);

	/*sequence one (seq1) is placed vertically and sequence two (seq2) horizontally*/
#pragma omp parallel firstprivate(addrHA, addrHB, addrHC, addrF, addrEO, addrEC) default(shared)
	for(int32_t diagonal = 2; diagonal <= s1Length + s2Length; diagonal++) { /*for each minor diagonal*/
		/* startRow and endRow are indexed from 1;
		 * calculate the starting and end column indices by the function x + y = diagonal (0<=x<=s2Length, 0<=y<=s1Length);
		 * */
		const register int32_t startRow = max(1, diagonal - s2Length); /*intersection point with line y = 1 or x = s2Length*/
		const register int32_t endRow = min(diagonal - 1, s1Length); /*intersection point with line x = 1 or y = s1Length*/
		const register __m512i vecEndRow = _mm512_set1_epi32(endRow);
		//printf("diagonal: %d startRow: %d endRow: %d\n", diagonal, startRow, endRow);

		/*get the thread ID*/
		const register int32_t tid = omp_get_thread_num();
		register int32_t col, row2;
		register __mmask16 vecMask;
		register __m512i vecData1, vecData2;
		register __m512i vecV, vecH, vecD;
		register __m512i localVecMaxScore = vecMaxScores[tid];

		/*synchronize all threads*/
#pragma omp barrier

		/*start the core loop*/
		const register int32_t numVectors = (endRow - startRow + 16) >> 4;
#pragma omp for schedule(static, (numVectors + numThreads - 1) / numThreads) nowait
		for(int32_t row = startRow; row <= endRow; row += 16) {
			/*mask the lanes which are out of the range*/
			vecMask = _mm512_cmple_epi32_mask(_mm512_add_epi32(vecIndices, _mm512_set1_epi32(row)), vecEndRow);
			//printf("============vecMask 0x%x row: %d endRow: %d tid %d\n", _mm512_mask2int(vecMask), row, endRow, tid);

			/*from left cell*/
			vecData1 = _mm512_mask_loadunpacklo_epi32(vecZero, vecMask, addrF + row);
			vecData1 = _mm512_mask_loadunpackhi_epi32(vecData1, vecMask, addrF + row + 16);
			vecData2 = _mm512_mask_loadunpacklo_epi32(vecZero, vecMask, addrHB + row);
			vecData2 = _mm512_mask_loadunpackhi_epi32(vecData2, vecMask, addrHB + row + 16);

			vecH = _mm512_mask_max_epi32(vecZero, vecMask, _mm512_mask_sub_epi32(vecZero, vecMask, vecData1, vecGapExtend), _mm512_mask_sub_epi32(vecZero, vecMask, vecData2, vecGapOE));

			/*from upper cell*/
			row2 = row - 1;
#if 0
			vecData1 = _mm512_mask_loadunpacklo_epi32(vecZero, vecMask, addrHB + row2);
			vecData1 = _mm512_mask_loadunpackhi_epi32(vecData1, vecMask, addrHB + row2 + 16);
#else
			/*we can ensure that the (row-1)th element is effective*/
			vecData1 = _mm512_mask_mov_epi32(_mm512_permutevar_epi32(vecPermuteIndex, vecData2), _mm512_int2mask(1), _mm512_set1_epi32(addrHB[row2]));
#endif
			vecData2 = _mm512_mask_loadunpacklo_epi32(vecZero, vecMask, addrEO + row2);
			vecData2 = _mm512_mask_loadunpackhi_epi32(vecData2, vecMask, addrEO + row2 + 16);

			vecV = _mm512_mask_max_epi32(vecZero, vecMask, _mm512_mask_sub_epi32(vecZero, vecMask, vecData2, vecGapExtend), _mm512_mask_sub_epi32(vecZero, vecMask, vecData1, vecGapOE));

			/*get characters from sequence 1*/
			vecData1 = _mm512_mask_extloadunpacklo_epi32(vecZero, vecMask, seq1 + row2, _MM_UPCONV_EPI32_SINT8, _MM_HINT_NONE);
			vecData1 = _mm512_mask_extloadunpackhi_epi32(vecData1, vecMask, seq1 + row2 + 64, _MM_UPCONV_EPI32_SINT8, _MM_HINT_NONE);

			/*get characters from sequence 2*/
			col = s2Length - (diagonal - row);
			vecData2 = _mm512_mask_extloadunpacklo_epi32(vecZero, vecMask, seq2Rev + col, _MM_UPCONV_EPI32_SINT8, _MM_HINT_NONE);
			vecData2 = _mm512_mask_extloadunpackhi_epi32(vecData2, vecMask, seq2Rev + col + 64, _MM_UPCONV_EPI32_SINT8, _MM_HINT_NONE);
			/*printf("col: %d diagonal - row: %d row: %d\n", col,diagonal - row, row);
			 printVector(vecData1, "sequence 1: ");
			 printVector(vecData2, "sequence 2: ");*/

			/*get the substitution score*/
			vecData1= _mm512_mask_mov_epi32(vecMismatch, _mm512_mask_cmpeq_epi32_mask(vecMask, vecData1, vecData2), vecMatch);
			//printVector(vecData1, "substitution score: ");

			/*diagonal cell*/
			vecData2 = _mm512_mask_loadunpacklo_epi32(vecZero, vecMask, addrHA + row2);
			vecData2 = _mm512_mask_loadunpackhi_epi32(vecData2, vecMask, addrHA + row2 + 16);

			vecD = _mm512_mask_add_epi32(vecZero, vecMask, vecData1, vecData2);
			vecD = _mm512_mask_max_epi32(vecZero, vecMask, _mm512_mask_max_epi32(vecZero, vecMask, vecH, vecV), _mm512_mask_max_epi32(vecZero, vecMask, vecD, vecZero));
			//printVector(vecD, "new alignment score: ");

			/*calculate the maximum alignment score*/
			localVecMaxScore = _mm512_mask_max_epi32(localVecMaxScore, vecMask, localVecMaxScore, vecD);

			/*save alignment scores*/
			_mm512_mask_packstorelo_epi32(addrHC + row, vecMask, vecD);
			_mm512_mask_packstorehi_epi32(addrHC + row + 16, vecMask, vecD);
			_mm512_mask_packstorelo_epi32(addrF + row, vecMask, vecH);
			_mm512_mask_packstorehi_epi32(addrF + row + 16, vecMask, vecH);
			_mm512_mask_packstorelo_epi32(addrEC + row, vecMask, vecV);
			_mm512_mask_packstorehi_epi32(addrEC + row + 16, vecMask, vecV);
		} /*there is an implicit barrier at the end of the loop*/

		/*save the maximum alignment score*/
		vecMaxScores[tid] = localVecMaxScore;

		//swap the buffers A <- B; B <- C;
		int32_t* __restrict__ tmp = addrHA;
		addrHA = addrHB;
		addrHB = addrHC;
		addrHC = tmp;

		tmp = addrEC;
		addrEC = addrEO;
		addrEO = tmp;
	}
	/*get the maximum alignment score*/
	maxScore = 0;
	for(int32_t i = 0; i < numThreads; ++i) {
		maxScore = max(maxScore, _mm512_reduce_max_epi32(vecMaxScores[i]));
	}

	/*release memory*/
	_mm_free(vecMaxScores);
	_mm_free(addrHA);
	_mm_free(addrHB);
	_mm_free(addrHC);
	_mm_free(addrEO);
	_mm_free(addrEC);
	_mm_free(addrF);

#else
	maxScore = -1; /*executed on the host*/
#endif
	return maxScore;
}

