/*
 * AlignCoreTiling.cpp
 *
 *  Created on: Feb 25, 2014
 *      Author: yongchao
 */
#include "AlignCoreTiling.h"

/*alignment entry function*/
void AlignCoreTiling::align() {

	/*checking if all the conditions have been met*/
	if ((_alpha & 15) != 0) {
		if (_rank == 0)
			Utils::log("The vertical tile size must be a multiple of 16\n");
		return;
	}
	if (_beta && (_beta & 15) != 0) {
		if (_rank == 0)
			Utils::log("The horizontal tile size must be a multiple of 16\n");
		return;
	}

	if ((_seq1Length & 15) != 0 || (_seq2Length & 15) != 0) {
		if (_rank == 0)
			Utils::log(
					"ERROR: both of the sequence lengths must be multipls of 16\n");
		return;
	}

#ifdef MPI_PARALLEL
	/*adjust the vertical block size*/
	_mpiAlpha = min(_mpiAlpha, _seq1Length); /*must be a multiple of 16*/
	_mpiBeta = (_seq2Length + _numProcs - 1) / _numProcs;
	_mpiBeta = (_mpiBeta + 15) >> 4 << 4; /*tile size is aligned to 16*/

	/*the number of blocks in the vertical direction*/
	const int32_t vMpiNumBlocks = (_seq1Length + _mpiAlpha - 1) / _mpiAlpha; /*vertical*/
	/*the number of blocks in the horizontal direction*/
	const int32_t hMpiNumBlocks = (_seq2Length + _mpiBeta - 1) / _mpiBeta; /*horizontal*/
	/*compute the number of effective processes*/
	const int32_t realNumProcs = min(hMpiNumBlocks, _numProcs);

	/*calculate the vertical tile size*/
	if(_beta == 0) {
		_beta = (_mpiBeta + _numMicThreads - 1) / _numMicThreads;
		_beta = (_beta + 15) >> 4 << 4;
	}
	/*statistical information*/
	if (_rank == 0) {
		Utils::log("Horizontal block size: %d\nVertical block size:%d\n", _mpiBeta, _mpiAlpha);
		Utils::log("Number of horizontal blocks: %d\nNumber of vertical blocks: %d\n", hMpiNumBlocks,
				vMpiNumBlocks);
		Utils::log("Horizontal tile size: %d\nVertical tile size: %d\n", _beta, _alpha);
		Utils::log("Number of horizontal tiles per block: %d\nNumber of vertical tiles per block: %d\n", (_mpiBeta + _beta - 1) / _beta, (_mpiAlpha + _alpha - 1) / _alpha);
		Utils::log("Number of effective processes: %d\n", realNumProcs);
	}

	/*invoke the kernels*/
	if (realNumProcs == 1) {
		_alignShared();
	} else {
		_alignDistributed(vMpiNumBlocks, hMpiNumBlocks, realNumProcs);
	}
#else
	/*calculate the vertical tile size*/
	if (_beta == 0) {
		_beta = (_seq2Length + _numMicThreads - 1) / _numMicThreads;
		_beta = (_beta + 15) >> 4 << 4;
	}

	/*statistical information*/
	Utils::log("Horizontal tile size: %d\nVertical tile size: %d\n", _beta,
			_alpha);
	Utils::log("Number of horizontal tiles: %d\nNumber of vertical tiles: %d\n",
			(_seq2Length + _beta - 1) / _beta,
			(_seq1Length + _alpha - 1) / _alpha);
	Utils::log("User-specified number of threads: %d\n", _numMicThreads);
	_alignShared();
#endif
}

void AlignCoreTiling::_alignDistributed(const int32_t vMpiNumBlocks,
		const int32_t hMpiNumBlocks, const int32_t realNumProcs) {
#ifdef MPI_PARALLEL
	double runtime;
	int32_t* __restrict__ mpiHE;
	int32_t* __restrict__ mpiHF;
	int32_t* __restrict__ mpiHFC;

	/*get the runtime*/
	runtime = MPI_Wtime();

	/*sequences*/
	int8_t *seq1 = _seq1;
	int8_t *seq2 = _seq2;

	/*the horizontal H and F vector sizes for the current process*/
	const int32_t heMpiVecSize = _mpiBeta;
	/*the vertical H and E vector for the current process*/
	const int32_t hfMpiVecSize = _mpiAlpha;

	/*allocate global buffers on the host*/
	mpiHE = (int32_t*) _mm_malloc(heMpiVecSize * sizeof(int2), 64);
	mpiHF = (int32_t*) _mm_malloc((hfMpiVecSize * 2 + 1) * sizeof(int32_t), 64);
	mpiHFC = (int32_t*) _mm_malloc((hfMpiVecSize * 2 + 1) * sizeof(int32_t), 64);

	/*initialize the buffer*/
#pragma simd
	for(int32_t i = 0; i < heMpiVecSize * 2; ++i) {
		mpiHE[i] = 0;
	}
#pragma simd
	for(int32_t i = 0; i <= hfMpiVecSize * 2; ++i) {
		mpiHF[i] = 0;
		mpiHFC[i] = 0;
	}

	/*transfer the sequences*/
#pragma offload_transfer target(mic:_micIndex) \
    in(seq1: length(_seq1Length) MIC_ALLOC) \
    in(seq2: length(_seq2Length) MIC_ALLOC)

	/*allocate buffers on the device*/
#pragma offload_transfer target(mic:_micIndex) \
	in(mpiHE: length(heMpiVecSize * 2) MIC_ALLOC) \
	in(mpiHF: length(hfMpiVecSize * 2 + 1) MIC_ALLOC) \
	nocopy(mpiHFC: length(hfMpiVecSize * 2 + 1) MIC_ALLOC)

	/*the core loop*/
	int32_t maxScore = 0;
	int32_t localMaxScore = 0;
	for (int32_t diagonal = 0; diagonal < vMpiNumBlocks + hMpiNumBlocks - 1;
			diagonal++) { /*for each minor diagonal; indexed from zero*/
		/* startBlockRow and endBlockRow are indexed from 1;
		 * calculate the starting and end column indices by the function x + y = diagonal (0<=x<s2Length, 0<=y<s1Length);
		 * */
		const int32_t startMpiRow = max(0, diagonal - hMpiNumBlocks + 1); /*intersection point with line y = 0 or x = s1Length - 1*/
		const int32_t endMpiRow = min(diagonal, vMpiNumBlocks - 1); /*intersection point with line x = 0 or y = s2Length - 1*/
		const int32_t startMpiCol = diagonal - endMpiRow;
		const int32_t endMpiCol = diagonal - startMpiRow;
		//printf("startMpiCol: %d endMpiCol: %d\n", startMpiCol, endMpiCol);

		/*get the block coordinate for each process*/
		const int32_t currentMpiCol = _rank;
		const int32_t currentMpiRow = diagonal - currentMpiCol;

		/*synchronize all processes*/
		MPI_Barrier (MPI_COMM_WORLD);

		/*invoke the kernel*/
		if(currentMpiCol >= startMpiCol && currentMpiCol <= endMpiCol) {
			/*invoke the kernel*/
#pragma offload target(mic:_micIndex) \
		nocopy(seq1: MIC_REUSE) \
		nocopy(seq2: MIC_REUSE) \
		nocopy(mpiHE: MIC_REUSE) \
		in(mpiHF: length (hfMpiVecSize * 2 + 1) MIC_REUSE)\
		out(mpiHFC: length(hfMpiVecSize * 2 + 1) MIC_REUSE)
			{
				int32_t score, mpiHDC;
				/*compute the global offset*/
				const int32_t globalRowOff = currentMpiRow * _mpiAlpha;
				const int32_t globalColOff = currentMpiCol * _mpiBeta;
				//printf("rank: %d currentMpiRow: %d currentMpiCol: %d globalRowOff: %d globalColOff: %d\n", _rank, currentMpiRow, currentMpiCol, globalRowOff, globalColOff);

				/*compute the segment lengths for S1 and S2*/
				const int32_t seq1SegmentLength = min(_seq1Length - globalRowOff, _mpiAlpha);
				const int32_t seq2SegmentLength = min(_seq2Length - globalColOff, _mpiBeta);
				//printf("rank: %d seq1SegmentLength: %d seq2SegmentLength: %d\n", _rank, seq1SegmentLength, seq2SegmentLength);

				/*invoke the kernel*/
				score = _swKernelDistributed(seq1 + globalRowOff, seq1SegmentLength, seq2 + globalColOff, seq2SegmentLength,
						(int2*)mpiHE, mpiHF, mpiHFC, mpiHF + hfMpiVecSize, mpiHFC + hfMpiVecSize, mpiHF[hfMpiVecSize * 2], &mpiHDC, _gapOE, _gapExtend, _match, _mismatch,
						_numMicThreads, _alpha, _beta);
				//printf("rank: %d score: %d\n", _rank, score);

				/*calculate the maximum local alignment score*/
				localMaxScore = max(localMaxScore, score);

				/*write the returned mpiHD at the end of mpiHFC buffer*/
				mpiHFC[hfMpiVecSize * 2] = mpiHDC;
			}/*offload region*/
		}

		/*communicate the new horizontal H and F values, as well as the diagonal value*/
		_shiftHFBuffers(mpiHF, mpiHFC, hfMpiVecSize, realNumProcs);
	}

	/*release the sequence on the device*/
#pragma offload_transfer target(mic:_micIndex) \
	nocopy(seq1: MIC_FREE) \
	nocopy(seq2: MIC_FREE) \
	nocopy(mpiHE : MIC_FREE) \
	nocopy(mpiHF : MIC_FREE) \
	nocopy(mpiHFC: MIC_FREE)

	/*release host buffers*/
	_mm_free(mpiHE);
	_mm_free(mpiHF);
	_mm_free(mpiHFC);

	/*Reduce the get the optimal local alignment scores*/
	if (MPI_Reduce(&localMaxScore, &maxScore, 1, MPI_INT, MPI_MAX, 0,
					MPI_COMM_WORLD) != MPI_SUCCESS) {
		Utils::exit("MPI_Reduce failed for process %d at line %d in file %s\n",
				_rank, __LINE__, __FILE__);
	}

	/*get the runtime*/
	runtime = MPI_Wtime() - runtime;

	/*report the result*/
	if (_rank == 0) {
		Utils::log("Computation time on the Xeon Phi: %f seconds\n", runtime);

		/*compute GCUPS*/
		double product = _seq1Length;
		product *= _seq2Length;
		Utils::log("Performance in GCUPS: %f GCUPS\n",
				product / runtime / 1000000 / 1000);
		Utils::log("Optimal local alignment score: %d\n", maxScore);
	}
#else
	Utils::log("Invoked the wrong function\n");
#endif
}

/*scale and vectorize optimization one*/
__attribute__((target(mic)))           int32_t AlignCoreTiling::_swKernelDistributed(
		const int8_t* __restrict__ seq1, const int32_t s1Length,
		const int8_t* __restrict__ seq2, const int32_t s2Length,
		int2* __restrict__ mpiHE, const int32_t* __restrict__ mpiH,
		int32_t* __restrict__ mpiHC, const int32_t* __restrict__ mpiF,
		int32_t* __restrict__ mpiFC, const int32_t mpiHD,
		int32_t* __restrict__ mpiHDC, const int32_t gapOE,
		const int32_t gapExtend, const int32_t match, const int32_t mismatch,
		const int32_t numThreads, const int32_t alpha, const int32_t beta) {

	int32_t maxScore;
#ifdef __MIC__
	__m512i* vecMaxScores;
	int2* __restrict__ globalHE;
	int32_t* __restrict__ globalH;
	int32_t* __restrict__ globalF;
	int32_t* __restrict__ globalHC;
	int32_t* __restrict__ globalFC;
	int32_t* __restrict__ globalHD; /*each tile has an diagonal value for the left-upper cell*/
	int32_t* __restrict__ globalHDC; /*each tile has an diagonal value for the left-upper cell*/
	const register int32_t vNumBlocks = (s1Length + alpha - 1) / alpha; /*vertical*/
	const register int32_t hNumBlocks = (s2Length + beta - 1) / beta; /*horizontal*/
	const register int32_t hfVecSize = (hNumBlocks + 1) * alpha;

	globalHE = mpiHE; /*no need to allocate new memory*/
	globalH = (int32_t*)_mm_malloc(hfVecSize * sizeof(int32_t), 64);
	globalHC = (int32_t*)_mm_malloc(hfVecSize * sizeof(int32_t), 64);
	globalF = (int32_t*)_mm_malloc(hfVecSize * sizeof(int32_t), 64);
	globalFC = (int32_t*)_mm_malloc(hfVecSize * sizeof(int32_t), 64);
	globalHD = (int32_t*)_mm_malloc((hNumBlocks + 1) * sizeof(int32_t), 64);
	globalHDC = (int32_t*)_mm_malloc((hNumBlocks + 1) * sizeof(int32_t), 64);

	/*initialize invariant values*/
	const register __m512i vecGapOE = _mm512_set1_epi32(gapOE);
	const register __m512i vecGapExtend = _mm512_set1_epi32(gapExtend);
	const register __m512i vecMatch = _mm512_set1_epi32(match);
	const register __m512i vecMismatch = _mm512_set1_epi32(mismatch);
	const register __m512i vecZero = _mm512_setzero_epi32();
	const register __m512i vecIndices = _mm512_set_epi32(15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0); /*from high to low*/
	const register __m512i vecDummy = _mm512_set1_epi32(DUMMY_NUCLEOTIDE);

	/*set the number of threads*/
	const int32_t numThreadsUsed = min(hNumBlocks, numThreads);
	omp_set_num_threads(numThreadsUsed);

	/*allocate maximum score buffer and initialize them*/
	vecMaxScores = (__m512i*)_mm_malloc(sizeof(__m512i) * numThreadsUsed, 64);
	for(int32_t i = 0; i < numThreadsUsed; ++i) {
		vecMaxScores[i] = _mm512_setzero_epi32();
	}

	/*save the next HD value*/
	*mpiHDC = globalHE[s2Length - 1]._x;

	/*sequence one (seq1) is placed vertically and sequence two (seq2) horizontally; indexed from zero*/
#pragma omp parallel firstprivate(globalH, globalHC, globalF, globalFC, globalHD, globalHDC) default(shared)
	for(int32_t outdiag = 0; outdiag < vNumBlocks + hNumBlocks - 1; outdiag++) { /*for each minor diagonal; indexed from zero*/
		/* startBlockRow and endBlockRow are indexed from 1;
		 * calculate the starting and end column indices by the function x + y = diagonal (0<=x<s2Length, 0<=y<s1Length);
		 * */
		const register int32_t startBlockRow = max(0, outdiag - hNumBlocks + 1); /*intersection point with line y = 0 or x = s1Length - 1*/
		const register int32_t endBlockRow = min(outdiag, vNumBlocks - 1); /*intersection point with line x = 0 or y = s2Length - 1*/

		//printf("startBlockRow: %d endBlockRow: %d\n", startBlockRow, endBlockRow);
		/*get the thread ID*/
		const register int32_t tid = omp_get_thread_num();
		register int32_t row, vTileSize, hTileSize, tileRow, tileCol, startGlobalCol;
		int8_t* __restrict__ seq2Ptr;
		register __mmask16 vecMask;
		register __m512i vecSeq1Symbols, vecSeq2Symbols, vecSubScore;
		register __m512i vecF, vecH, vecE, vecEUp, vecHUp, vecHUp2;
		register __m512i localVecMaxScore = vecMaxScores[tid];
		int32_t buffer[16] __attribute__((aligned(64)));
		int32_t buffer2[16] __attribute__((aligned(64)));

		/*initialize the first columns of the related buffers*/
#pragma omp single
		if(outdiag - endBlockRow == 0) {
			row = endBlockRow * alpha;
			for(tileRow = 0; tileRow < alpha; tileRow += 16) {
				_mm512_store_epi32(globalH + tileRow, _mm512_load_epi32(mpiH + row + tileRow));
				_mm512_store_epi32(globalF + tileRow, _mm512_load_epi32(mpiF + row + tileRow));
			}
			globalHD[0] = endBlockRow ? mpiH[row - 1] : mpiHD;
		} /*an implicit barrier for all threads*/

		/*start the core loop*/
#pragma omp for schedule(static, 1) nowait
		for(int32_t blockRow = startBlockRow; blockRow <= endBlockRow; ++blockRow) {
			/*for the current thread, we calculate the row range for its corresponding tile*/
			const int32_t blockCol = outdiag - blockRow;
			const int32_t startGlobalRow = blockRow * alpha;
			/*for each vector length*/
			hTileSize = min(beta, s2Length - blockCol * beta);
			vTileSize = min(alpha, s1Length - blockRow * alpha);

			/*write the diagonal value to the buffer of the right neighboring tile*/
			globalHDC[blockCol + 1] = globalHE[blockCol * beta + hTileSize - 1]._x;

			/*compute the global column index*/
			startGlobalCol = blockCol * beta;
			for(tileRow = 0; tileRow < vTileSize; tileRow += 16) {

				/*get characters from sequence 1*/
				row = blockRow * alpha + tileRow;
				vecSeq1Symbols = _mm512_extloadunpacklo_epi32(vecZero, seq1 + row, _MM_UPCONV_EPI32_SINT8, _MM_HINT_NONE);
				vecSeq1Symbols = _mm512_extloadunpackhi_epi32(vecSeq1Symbols, seq1 + row + 64, _MM_UPCONV_EPI32_SINT8, _MM_HINT_NONE);
				vecSeq2Symbols = vecDummy;

				/*load the left apron values*/
				row = blockCol * alpha + tileRow; /*get the H and F row index*/
				vecH = _mm512_load_epi32(globalH + row);
				vecF = _mm512_load_epi32(globalF + row);
				//printf("startGlobalCol: %d startGlobalRow: %d tileRow: %d\n", startGlobalCol, startGlobalRow, tileRow);
				//printVector(vecH, "Horizontal H\n");
				//printVector(vecF, "Horizontal F\n");
				/*form the diagonal values*/
				vecHUp2 = tileRow ? _mm512_set1_epi32(globalH[tileRow - 1]) : _mm512_set1_epi32(globalHD[blockCol]);
				//printf("globalHDC[0]: %d\n", globalHDC[0]);
				//printf("blockRow: %d blockCol: %d tileRow: %d mpiH index %d mpiH %d\n", blockRow, blockCol, tileRow, startGlobalRow + tileRow - 1, mpiH[startGlobalRow + tileRow - 1]);
				//printVector(vecHUp2, "Diagonal value\n");
				vecE = vecZero;

				/*set sequence 2 pointer*/
				seq2Ptr = (int8_t*)&seq2[startGlobalCol];
				for(tileCol = 0; tileCol < 16; ++tileCol, ++startGlobalCol, ++seq2Ptr) {

					/*calculate a mask*/
					vecMask = _mm512_cmple_epi32_mask(vecIndices, _mm512_set1_epi32(tileCol));
					//printf("vecMask: 0x%x, tileCol: %d\n", vecMask, tileCol);

					/*load the upper H and E values*/
					vecHUp = _mm512_set1_epi32(globalHE[startGlobalCol]._x);
					vecEUp = _mm512_set1_epi32(globalHE[startGlobalCol]._y);
					//printf("Load startGlobalCol: %d E: %d Hup: %d\n", startGlobalCol, globalHE[startGlobalCol]._y, globalHE[startGlobalCol]._x);

					/*shift the vector: put the new E and H value at the lowest element of E and H*/
					vecEUp = _mm512_alignr_epi32(vecE, vecEUp, 15);
					vecHUp = _mm512_alignr_epi32(vecH, vecHUp, 15);

					/*update E value*/
					vecE = _mm512_max_epi32(_mm512_sub_epi32(vecEUp, vecGapExtend), _mm512_sub_epi32(vecHUp, vecGapOE));

					/*update F value*/
					vecF = _mm512_max_epi32(_mm512_sub_epi32(vecF, vecGapExtend), _mm512_sub_epi32(vecH, vecGapOE));

					/*update the H value*/
					/*get characters from sequence 2*/
					vecSeq2Symbols = _mm512_alignr_epi32(vecSeq2Symbols, _mm512_set1_epi32(*seq2Ptr), 15);

					//printVector(vecSeq1Symbols, "sequence 1\n");
					//printVector(vecSeq2Symbols, "sequence 2\n");

					/*get the substitution score*/
					vecSubScore = _mm512_mask_mov_epi32(vecMismatch, _mm512_mask_cmpeq_epi32_mask(vecMask, vecSeq1Symbols, vecSeq2Symbols), vecMatch);
					vecSubScore = _mm512_mask_mov_epi32(vecSubScore, _mm512_kor(_mm512_mask_cmpeq_epi32_mask(vecMask, vecSeq1Symbols, vecDummy), _mm512_mask_cmpeq_epi32_mask(vecMask, vecSeq2Symbols, vecDummy)), vecZero);
					//printVector(vecData1, "substitution score\n");

					//printVector(vecHUp2, "diagonal score\n");
					vecH = _mm512_mask_add_epi32(vecH, vecMask, vecHUp2, vecSubScore);
					vecH = _mm512_mask_max_epi32(vecZero, vecMask, _mm512_mask_max_epi32(vecZero, vecMask, vecE, vecF), _mm512_max_epi32(vecH, vecZero));
					//printVector(vecH, "Kernel one: alignment score\n");

					/*calculate the maximum alignment score*/
					localVecMaxScore = _mm512_mask_max_epi32(localVecMaxScore, vecMask, localVecMaxScore, vecH);

					/*save the H value for the next diagonal*/
					vecHUp2 = _mm512_mask_mov_epi32(_mm512_set1_epi32(globalH[row + tileCol]), vecMask, vecHUp);
				}
				/*save the H and E values for the first triangle*/
				_mm512_store_epi32(buffer, vecH);
				_mm512_store_epi32(buffer2, vecE);
				globalHE[startGlobalCol - 16] = int2(buffer[15], buffer2[15]);
				//printf("startGlobalCol: %d saved E %d H %d\n", startGlobalCol - 16, globalHE[startGlobalCol - 16]._y, globalHE[startGlobalCol - 16]._x);

				//printf("in the middle\n");
				for(; tileCol < hTileSize; ++tileCol, ++startGlobalCol, ++seq2Ptr) {

					/*load the upper H and E values*/
					vecHUp = _mm512_set1_epi32(globalHE[startGlobalCol]._x);
					vecEUp = _mm512_set1_epi32(globalHE[startGlobalCol]._y);
					//printf("Load startGlobalCol: %d E: %d Hup: %d\n", startGlobalCol, globalHE[startGlobalCol]._y, globalHE[startGlobalCol]._x);

					/*shift the vector*/
					vecEUp = _mm512_alignr_epi32(vecE, vecEUp, 15);
					vecHUp = _mm512_alignr_epi32(vecH, vecHUp, 15);

					/*update E value*/
					vecE = _mm512_max_epi32(_mm512_sub_epi32(vecEUp, vecGapExtend), _mm512_sub_epi32(vecHUp, vecGapOE));

					/*update F value*/
					vecF = _mm512_max_epi32(_mm512_sub_epi32(vecF, vecGapExtend), _mm512_sub_epi32(vecH, vecGapOE));

					/*update the H value*/
					/*get characters from sequence 2*/
					vecSeq2Symbols = _mm512_alignr_epi32(vecSeq2Symbols, _mm512_set1_epi32(*seq2Ptr), 15);
					//printVector(vecSeq1Symbols, "sequence 1\n");
					//printVector(vecSeq2Symbols, "sequence 2\n");

					/*get the substitution score*/
					vecSubScore = _mm512_mask_mov_epi32(vecMismatch, _mm512_cmpeq_epi32_mask(vecSeq1Symbols, vecSeq2Symbols), vecMatch);
					vecSubScore = _mm512_mask_mov_epi32(vecSubScore, _mm512_kor(_mm512_cmpeq_epi32_mask(vecSeq1Symbols, vecDummy), _mm512_cmpeq_epi32_mask(vecSeq2Symbols, vecDummy)), vecZero);
					//printVector(vecSubScore, "substitution score\n");

					vecH = _mm512_add_epi32(vecHUp2, vecSubScore);
					vecH = _mm512_max_epi32(_mm512_max_epi32(vecE, vecF), _mm512_max_epi32(vecH, vecZero));
					//printVector(vecHUp2, "vecHUp2: \n");
					//printVector(vecH, "Kernel two vecH: \n");

					/*calculate the maximum alignment score*/
					localVecMaxScore = _mm512_max_epi32(localVecMaxScore, vecH);

					/*save the H value for the next diagonal*/
					vecHUp2 = vecHUp;

					/*save the new H and E values*/
					_mm512_store_epi32(buffer, vecH);
					_mm512_store_epi32(buffer2, vecE);
					globalHE[startGlobalCol- 15] = int2(buffer[15], buffer2[15]);
					//printf("StartGlobalCol: %d E %d H %d\n", startGlobalCol - 15, globalHE[startGlobalCol - 15]._y, globalHE[startGlobalCol - 15]._x);
				}
				//printf("tileCol: %d beta: %d\n", tileCol, beta);
				//
				/*save the horizontal H and F values*/
				_mm512_store_epi32(buffer, vecH);
				globalHC[(blockCol + 1) * alpha + tileRow] = buffer[0];
				_mm512_store_epi32(buffer2, vecF);
				globalFC[(blockCol + 1) * alpha + tileRow] = buffer2[0];
				/*save to the MPI buffers*/
				if(blockCol + 1 == hNumBlocks) {
					mpiHC[startGlobalRow + tileRow] = buffer[0];
					mpiFC[startGlobalRow + tileRow] = buffer2[0];
					//printf("mpiHC %d mpiFC %d\n", buffer[0], buffer2[0]);
				}

				/*shift the diagonal value to the left*/
				for(int32_t i = 1; i <= 15; ++i, ++startGlobalCol) {

					/*calculate a mask*/
					vecMask = _mm512_cmpge_epi32_mask(vecIndices, _mm512_set1_epi32(i));
					//printf("startGlobalCol: %d vecMask: 0x%x\n", startGlobalCol, vecMask);
					/*shift the vector*/
					vecEUp = _mm512_alignr_epi32(vecE, vecZero, 15);
					vecHUp = _mm512_alignr_epi32(vecH, vecZero, 15);

					/*update E value*/
					vecE = _mm512_max_epi32(_mm512_sub_epi32(vecEUp, vecGapExtend), _mm512_sub_epi32(vecHUp, vecGapOE));

					/*update F value*/
					vecF = _mm512_max_epi32(_mm512_sub_epi32(vecF, vecGapExtend), _mm512_sub_epi32(vecH, vecGapOE));

					/*update the H value*/

					/*get characters from sequence 2*/
					vecSeq2Symbols = _mm512_alignr_epi32(vecSeq2Symbols, vecDummy, 15);
					//printVector(vecSeq1Symbols, "sequence 1\n");
					//printVector(vecSeq2Symbols, "sequence 2\n");

					/*get the substitution score*/
					vecSubScore = _mm512_mask_mov_epi32(vecMismatch, _mm512_mask_cmpeq_epi32_mask(vecMask, vecSeq1Symbols, vecSeq2Symbols), vecMatch);
					vecSubScore = _mm512_mask_mov_epi32(vecSubScore, _mm512_kor(_mm512_cmpeq_epi32_mask(vecSeq1Symbols, vecDummy), _mm512_cmpeq_epi32_mask(vecSeq2Symbols, vecDummy)), vecZero);
					//printVector(vecSubScore, "substitution score\n");

					vecH = _mm512_mask_add_epi32(vecH, vecMask, vecHUp2, vecSubScore);
					vecH = _mm512_mask_max_epi32(vecH, vecMask, _mm512_mask_max_epi32(vecZero, vecMask, vecE, vecF), _mm512_max_epi32(vecH, vecZero));
					//printVector(vecHUp2, "vecHUp2:\n");
					//printVector(vecH, "Kernel three vecH: \n");

					/*calculate the maximum alignment score*/
					localVecMaxScore = _mm512_max_epi32(localVecMaxScore, vecH);

					/*save the H value for the next diagonal*/
					vecHUp2 = vecHUp;

					/*save the new vertical H and Evalues*/
					_mm512_store_epi32(buffer, vecH);
					_mm512_store_epi32(buffer2, vecE);
					globalHE[startGlobalCol - 15] = int2(buffer[15], buffer2[15]);
					//printf("startGlobalCol: %d Saved E %d H %d\n", startGlobalCol - 15, globalHE[startGlobalCol - 15]._y, globalHE[startGlobalCol - 15]._x);

					/*save the horizontal H and F values*/
					row = (blockCol + 1) * alpha + tileRow + i;
					_mm512_store_epi32(buffer, vecH);
					globalHC[row] = buffer[i];
					_mm512_store_epi32(buffer2, vecF);
					globalFC[row] = buffer2[i];

					/*save to the MPI buffers*/
					if(blockCol + 1 == hNumBlocks) {
						row = startGlobalRow + tileRow + i;
						mpiHC[row] = buffer[i];
						mpiFC[row] = buffer2[i];
						//printf("mpiHC %d mpiFC %d\n", buffer[i], buffer2[i]);
					}
				}
			}
#if 0
			for(int i = 0; i < s2Length; ++i) {
				printf("i: %d E: %d H: %d\n", i, globalE[i], globalHup[i]);
			}
#endif

		} /*finish an external diagonal*/

		/*save the maximum alignment score*/
		vecMaxScores[tid] = localVecMaxScore;

		/*swap the HD buffers*/
		int32_t* __restrict__ tmp = globalHD;
		globalHD = globalHDC;
		globalHDC = tmp;

		tmp = globalF;
		globalF = globalFC;
		globalFC = tmp;

		tmp = globalH;
		globalH = globalHC;
		globalHC = tmp;
	}

	/*get the maximum alignment score*/
	maxScore = 0;
	for(int32_t i = 0; i < numThreadsUsed; ++i) {
		maxScore = max(maxScore, _mm512_reduce_max_epi32(vecMaxScores[i]));
	}

	/*release memory*/
	_mm_free(vecMaxScores);
	_mm_free(globalF);
	_mm_free(globalH);
	_mm_free(globalHC);
	_mm_free(globalFC);
	_mm_free(globalHD);
	_mm_free(globalHDC);
#else
	maxScore = -1; /*executed on the host*/
#endif
	return maxScore;
}

void AlignCoreTiling::_shiftHFBuffers(const int32_t* mpiHF,
		const int32_t* mpiHFC, const int32_t hfMpiVecSize,
		const int32_t realNumProcs) {

#ifdef MPI_PARALLEL
	const int32_t transSize = hfMpiVecSize * 2 + 1;
	/*for a single process*/
	if(realNumProcs == 1) {
		/*copy the data*/
		memcpy((void*)mpiHF, (void*)mpiHFC, transSize * sizeof(int32_t));
		return;
	}

	/*for multiple processes*/
	MPI_Request requests[2];
	const int32_t leftRank = _rank - 1;
	const int32_t rightRank = _rank + 1;
	if (_rank == 0) {
		/*send the data*/
		if (MPI_Send((void*)mpiHFC, transSize, MPI_INT, rightRank, _rank,
						MPI_COMM_WORLD) != MPI_SUCCESS) {
			Utils::exit(
					"MPI_Send failed for process %d at line %d in file %s\n",
					_rank, __LINE__, __FILE__);
		}
	} else if (_rank == realNumProcs - 1) {
		/*receive the data*/
		if (MPI_Recv((void*)mpiHF, transSize, MPI_INT, leftRank, leftRank,
						MPI_COMM_WORLD, MPI_STATUS_IGNORE) != MPI_SUCCESS) {
			Utils::exit(
					"MPI_Recv failed for process %d at line %d in file %s\n",
					_rank, __LINE__, __FILE__);
		}
	} else if (_rank > 0 && _rank < realNumProcs - 1) {
		/*send the data*/
		if (MPI_Isend((void*)mpiHFC, transSize, MPI_INT, rightRank, _rank,
						MPI_COMM_WORLD, &requests[0]) != MPI_SUCCESS) {
			Utils::exit(
					"MPI_Send failed for process %d at line %d in file %s\n",
					_rank, __LINE__, __FILE__);
		}

		/*receive the data*/
		if (MPI_Irecv((void*)mpiHF, transSize, MPI_INT, leftRank, leftRank,
						MPI_COMM_WORLD, &requests[1]) != MPI_SUCCESS) {
			Utils::exit(
					"MPI_Recv failed for process %d at line %d in file %s\n",
					_rank, __LINE__, __FILE__);
		}

		/*wait for the completion of both send and recv*/
		MPI_Wait(&requests[0], MPI_STATUS_IGNORE);
		MPI_Wait(&requests[1], MPI_STATUS_IGNORE);

	}
#endif
}
void AlignCoreTiling::_alignShared() {
	int32_t maxScore;
	double runtime = 0;
	int8_t* seq1 = _seq1;
	int8_t* seq2 = _seq2;

	/*get the start time*/
	runtime = omp_get_wtime();

	/*start to offload*/
#pragma offload target(mic:_micIndex) \
		in(seq1: length(_seq1Length) MIC_ALLOC_FREE) \
		in(seq2: length(_seq2Length) MIC_ALLOC_FREE)
	{

		maxScore = _swKernelShared(seq1, _seq1Length, seq2, _seq2Length, _gapOE,
				_gapExtend, _match, _mismatch, _numMicThreads, _alpha, _beta);
	}

	if (maxScore == -1) {
		Utils::exit("Errors have occurred in the offload region\n");
	}

	/*calculate the CPU time on the device*/
	runtime = omp_get_wtime() - runtime;

	Utils::log("Computation time on the Xeon Phi: %f seconds\n", runtime);

	/*compute GCUPS*/
	double product = _seq1Length;
	product *= _seq2Length;
	Utils::log("Performance in GCUPS: %f GCUPS\n",
			product / runtime / 1000000 / 1000);
	Utils::log("Optimal local alignment score: %d\n", maxScore);
}

/*scale and vectorize optimization one*/
__attribute__((target(mic)))           int32_t AlignCoreTiling::_swKernelShared(
		const int8_t* __restrict__ seq1, const int32_t s1Length,
		const int8_t* __restrict__ seq2, const int32_t s2Length,
		const int32_t gapOE, const int32_t gapExtend, const int32_t match,
		const int32_t mismatch, const int32_t numThreads, const int32_t alpha,
		const int32_t beta) {
	int32_t maxScore;

#ifdef __MIC__
	__m512i* vecMaxScores;
	int2* __restrict__ globalHE;
	int32_t* __restrict__ globalH;
	int32_t* __restrict__ globalF;
	int32_t* __restrict__ globalHC;
	int32_t* __restrict__ globalFC;
	int32_t* __restrict__ globalHD; /*each tile has an diagonal value for the left-upper cell*/
	int32_t* __restrict__ globalHDC; /*each tile has an diagonal value for the left-upper cell*/
	const register int32_t vNumBlocks = (s1Length + alpha - 1) / alpha; /*vertical*/
	const register int32_t hNumBlocks = (s2Length + beta - 1) / beta; /*horizontal*/
	const register int32_t hfVecSize = (hNumBlocks + 1) * alpha;

	/*allocate global buffers*/
	globalHE = (int2*)_mm_malloc(s2Length* sizeof(int2), 64);
#pragma simd
	for(int32_t i = 0; i < s2Length; ++i) {
		globalHE[i] = int2(0, 0);
	}

	globalH = (int32_t*)_mm_malloc(hfVecSize * sizeof(int32_t), 64);
	globalHC = (int32_t*)_mm_malloc(hfVecSize * sizeof(int32_t), 64);
	globalF = (int32_t*)_mm_malloc(hfVecSize * sizeof(int32_t), 64);
	globalFC = (int32_t*)_mm_malloc(hfVecSize * sizeof(int32_t), 64);
#pragma simd
	for(int32_t i = 0; i < hfVecSize; ++i) {
		globalH[i] = 0;
		globalHC[i] = 0;
		globalF[i] = 0;
		globalFC[i] = 0;
	}

	globalHD = (int32_t*)_mm_malloc((hNumBlocks + 1) * sizeof(int32_t), 64);
	globalHDC = (int32_t*)_mm_malloc((hNumBlocks + 1) * sizeof(int32_t), 64);
#pragma simd
	for(int32_t i = 0; i <= hNumBlocks; ++i) {
		globalHD[i] = 0;
		globalHDC[i] = 0;
	}

	/*initialize invariant values*/
	const register __m512i vecGapOE = _mm512_set1_epi32(gapOE);
	const register __m512i vecGapExtend = _mm512_set1_epi32(gapExtend);
	const register __m512i vecMatch = _mm512_set1_epi32(match);
	const register __m512i vecMismatch = _mm512_set1_epi32(mismatch);
	const register __m512i vecZero = _mm512_setzero_epi32();
	const register __m512i vecIndices = _mm512_set_epi32(15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0); /*from high to low*/
	const register __m512i vecDummy = _mm512_set1_epi32(DUMMY_NUCLEOTIDE);

	/*set the number of threads*/
	const int32_t numThreadsUsed = min(hNumBlocks, numThreads);
	omp_set_num_threads(numThreadsUsed);

	/*allocate maximum score buffer and initialize them*/
	vecMaxScores = (__m512i*)_mm_malloc(sizeof(__m512i) * numThreadsUsed, 64);
	for(int32_t i = 0; i < numThreadsUsed; ++i) {
		vecMaxScores[i] = _mm512_setzero_epi32();
	}

	/*sequence one (seq1) is placed vertically and sequence two (seq2) horizontally; indexed from zero*/
#pragma omp parallel firstprivate (globalH, globalHC, globalF, globalFC, globalHD, globalHDC) default(shared)
	for(int32_t outdiag = 0; outdiag < vNumBlocks + hNumBlocks - 1; outdiag++) { /*for each minor diagonal; indexed from zero*/
		/* startBlockRow and endBlockRow are indexed from 1;
		 * calculate the starting and end column indices by the function x + y = diagonal (0<=x<s2Length, 0<=y<s1Length);
		 * */
		const register int32_t startBlockRow = max(0, outdiag - hNumBlocks + 1); /*intersection point with line y = 0 or x = s1Length - 1*/
		const register int32_t endBlockRow = min(outdiag, vNumBlocks - 1); /*intersection point with line x = 0 or y = s2Length - 1*/

		//printf("startBlockRow: %d endBlockRow: %d\n", startBlockRow, endBlockRow);
		/*get the thread ID*/
		const register int32_t tid = omp_get_thread_num();
		register int32_t row, vTileSize, hTileSize, tileRow, tileCol, blockCol, startGlobalCol;
		int8_t* __restrict__ seq2Ptr;
		register __mmask16 vecMask;
		register __m512i vecSeq1Symbols, vecSeq2Symbols, vecSubScore;
		register __m512i vecF, vecH, vecE, vecEUp, vecHUp, vecHUp2;
		register __m512i localVecMaxScore = vecMaxScores[tid];
		int32_t buffer[16] __attribute__((aligned(64)));
		int32_t buffer2[16] __attribute__((aligned(64)));

		/*synchronize all threads*/
#pragma omp barrier

		/*start the core loop*/
#pragma omp for schedule(static, 1) nowait
		for(int32_t blockRow = startBlockRow; blockRow <= endBlockRow; ++blockRow) {
			/*for the current thread, we calculate the row range for its corresponding tile*/
			blockCol = outdiag - blockRow;
			//printf("tid: %d blockRow: %d blockCol: %d\n", tid, blockRow, blockCol);

			/*for each vector length*/
			hTileSize = min(beta, s2Length - blockCol * beta);
			vTileSize = min(alpha, s1Length - blockRow * alpha);

			/*write the diagonal value to the buffer of the right neighboring tile*/
			globalHDC[blockCol + 1] = globalHE[blockCol * beta + hTileSize - 1]._x;

			for(tileRow = 0; tileRow < vTileSize; tileRow += 16) {

				/*get characters from sequence 1*/
				row = blockRow * alpha + tileRow;
				vecSeq1Symbols = _mm512_extloadunpacklo_epi32(vecZero, seq1 + row, _MM_UPCONV_EPI32_SINT8, _MM_HINT_NONE);
				vecSeq1Symbols = _mm512_extloadunpackhi_epi32(vecSeq1Symbols, seq1 + row + 64, _MM_UPCONV_EPI32_SINT8, _MM_HINT_NONE);
				vecSeq2Symbols = vecDummy;

				/*load the left apron values*/
				row = blockCol * alpha + tileRow;
				vecH = _mm512_load_epi32(globalH + row);
				vecF = _mm512_load_epi32(globalF + row);
				//printVector(vecF, "F values:\n");
				vecE = vecZero;

				/*form the diagonal values*/
				vecHUp2 = tileRow ? _mm512_set1_epi32(globalH[row - 1]): _mm512_set1_epi32(globalHD[blockCol]);
				//printVector(vecHUp2, "Diagonal value\n");

				/*set sequence 2 pointer*/
				startGlobalCol = blockCol * beta;
				seq2Ptr = (int8_t*)&seq2[startGlobalCol];
				for(tileCol = 0; tileCol < 16; ++tileCol, ++startGlobalCol, ++seq2Ptr) {

					/*calculate a mask*/
					vecMask = _mm512_cmple_epi32_mask(vecIndices, _mm512_set1_epi32(tileCol));
					//printf("vecMask: 0x%x, tileCol: %d\n", vecMask, tileCol);

					/*load the upper H and E values*/
					vecHUp = _mm512_set1_epi32(globalHE[startGlobalCol]._x);
					vecEUp = _mm512_set1_epi32(globalHE[startGlobalCol]._y);
					//printf("Load startGlobalCol: %d E: %d Hup: %d\n", startGlobalCol, globalHE[startGlobalCol]._y, globalHE[startGlobalCol]._x);

					/*shift the vector: put the new E and H value at the lowest element of E and H*/
					vecEUp = _mm512_alignr_epi32(vecE, vecEUp, 15);
					vecHUp = _mm512_alignr_epi32(vecH, vecHUp, 15);

					/*update E value*/
					vecE = _mm512_max_epi32(_mm512_sub_epi32(vecEUp, vecGapExtend), _mm512_sub_epi32(vecHUp, vecGapOE));
					//printVector(vecE, "Kernel one vecE: \n");

					/*update F value*/
					vecF = _mm512_max_epi32(_mm512_sub_epi32(vecF, vecGapExtend), _mm512_sub_epi32(vecH, vecGapOE));
					//printVector(vecF, "Kernel one vecF: \n");

					/*update the H value*/
					/*get characters from sequence 2*/
					vecSeq2Symbols = _mm512_alignr_epi32(vecSeq2Symbols, _mm512_set1_epi32(*seq2Ptr), 15);
					//printVector(vecSeq1Symbols, "sequence 1: ");
					//printVector(vecSeq2Symbols, "sequence 2: ");

					/*get the substitution score*/
					vecSubScore = _mm512_mask_mov_epi32(vecMismatch, _mm512_mask_cmpeq_epi32_mask(vecMask, vecSeq1Symbols, vecSeq2Symbols), vecMatch);
					vecSubScore = _mm512_mask_mov_epi32(vecSubScore, _mm512_kor(_mm512_mask_cmpeq_epi32_mask(vecMask, vecSeq1Symbols, vecDummy), _mm512_mask_cmpeq_epi32_mask(vecMask, vecSeq2Symbols, vecDummy)), vecZero);
					//printVector(vecData1, "substitution score\n");

					//printVector(vecHUp2, "diagonal score\n");
					vecH = _mm512_mask_add_epi32(vecH, vecMask, vecHUp2, vecSubScore);
					vecH = _mm512_mask_max_epi32(vecZero, vecMask, _mm512_mask_max_epi32(vecZero, vecMask, vecE, vecF), _mm512_max_epi32(vecH, vecZero));
					/*if(tid == THREAD){
					 printVector(vecHUp2, "vecHUp2: ");
					 printVector(vecH, "Kernel one vecH: ");
					 printVector(vecE, "Kernel one vecE: ");
					 printVector(vecF, "Kernel one vecF: ");
					 }*/

					/*calculate the maximum alignment score*/
					localVecMaxScore = _mm512_mask_max_epi32(localVecMaxScore, vecMask, localVecMaxScore, vecH);

					/*save the H value for the next diagonal*/
					vecHUp2 = _mm512_mask_mov_epi32(_mm512_set1_epi32(globalH[row + tileCol]), vecMask, vecHUp);
				}
				/*save the H and E values for the first triangle*/
				_mm512_store_epi32(buffer, vecH);
				_mm512_store_epi32(buffer2, vecE);
				globalHE[startGlobalCol - 16] = int2(buffer[15], buffer2[15]);
				//printf("startGlobalCol: %d saved E %d H %d\n", startGlobalCol - 16, globalHE[startGlobalCol - 16]._y, globalHE[startGlobalCol - 16]._x);

				//printf("in the middle\n");
				for(; tileCol < hTileSize; ++tileCol, ++startGlobalCol, ++seq2Ptr) {

					/*load the upper H and E values*/
					vecHUp = _mm512_set1_epi32(globalHE[startGlobalCol]._x);
					vecEUp = _mm512_set1_epi32(globalHE[startGlobalCol]._y);
					//printf("Load startGlobalCol: %d E: %d Hup: %d\n", startGlobalCol, globalHE[startGlobalCol]._y, globalHE[startGlobalCol]._x);

					/*shift the vector*/
					vecEUp = _mm512_alignr_epi32(vecE, vecEUp, 15);
					vecHUp = _mm512_alignr_epi32(vecH, vecHUp, 15);

					/*update E value*/
					vecE = _mm512_max_epi32(_mm512_sub_epi32(vecEUp, vecGapExtend), _mm512_sub_epi32(vecHUp, vecGapOE));
					//printVector(vecE, "Kernel two vecE:\n");

					/*update F value*/
					vecF = _mm512_max_epi32(_mm512_sub_epi32(vecF, vecGapExtend), _mm512_sub_epi32(vecH, vecGapOE));
					//printVector(vecF, "Kernel two vecF:\n");

					/*update the H value*/
					/*get characters from sequence 2*/
					vecSeq2Symbols = _mm512_alignr_epi32(vecSeq2Symbols, _mm512_set1_epi32(*seq2Ptr), 15);
					//printVector(vecSeq1Symbols, "sequence 1\n");
					//printVector(vecSeq2Symbols, "sequence 2\n");

					/*get the substitution score*/
					vecSubScore = _mm512_mask_mov_epi32(vecMismatch, _mm512_cmpeq_epi32_mask(vecSeq1Symbols, vecSeq2Symbols), vecMatch);
					vecSubScore = _mm512_mask_mov_epi32(vecSubScore, _mm512_kor(_mm512_cmpeq_epi32_mask(vecSeq1Symbols, vecDummy), _mm512_cmpeq_epi32_mask(vecSeq2Symbols, vecDummy)), vecZero);
					//printVector(vecSubScore, "substitution score\n");

					vecH = _mm512_add_epi32(vecHUp2, vecSubScore);
					vecH = _mm512_max_epi32(_mm512_max_epi32(vecE, vecF), _mm512_max_epi32(vecH, vecZero));
					/*if(tid == THREAD){
					 printVector(vecHUp2, "vecHUp2: ");
					 printVector(vecH, "Kernel two vecH: ");
					 }*/

					/*calculate the maximum alignment score*/
					localVecMaxScore = _mm512_max_epi32(localVecMaxScore, vecH);

					/*save the H value for the next diagonal*/
					vecHUp2 = vecHUp;

					/*save the new H and E values*/
					_mm512_store_epi32(buffer, vecH);
					_mm512_store_epi32(buffer2, vecE);
					globalHE[startGlobalCol- 15] = int2(buffer[15], buffer2[15]);
					//printf("StartGlobalCol: %d E %d H %d\n", startGlobalCol - 15, globalHE[startGlobalCol - 15]._y, globalHE[startGlobalCol - 15]._x);
				}
				//printf("tileCol: %d beta: %d\n", tileCol, beta);
				//
				/*save the horizontal H and F values*/
				_mm512_store_epi32(buffer, vecH);
				globalHC[(blockCol + 1) * alpha + tileRow] = buffer[0];
				_mm512_store_epi32(buffer, vecF);
				globalFC[(blockCol + 1) * alpha + tileRow] = buffer[0];

				/*shift the diagonal value to the left*/
				for(int32_t i = 1; i <= 15; ++i, ++startGlobalCol) {

					/*calculate a mask*/
					vecMask = _mm512_cmpge_epi32_mask(vecIndices, _mm512_set1_epi32(i));
					//printf("startGlobalCol: %d vecMask: 0x%x\n", startGlobalCol, vecMask);
					/*shift the vector*/
					vecEUp = _mm512_alignr_epi32(vecE, vecZero, 15);
					vecHUp = _mm512_alignr_epi32(vecH, vecZero, 15);

					/*update E value*/
					vecE = _mm512_max_epi32(_mm512_sub_epi32(vecEUp, vecGapExtend), _mm512_sub_epi32(vecHUp, vecGapOE));
					//printVector(vecE, "Kernel three vecE:\n");

					/*update F value*/
					vecF = _mm512_max_epi32(_mm512_sub_epi32(vecF, vecGapExtend), _mm512_sub_epi32(vecH, vecGapOE));
					//printVector(vecF, "Kernel three vecF:\n");

					/*update the H value*/

					/*get characters from sequence 2*/
					vecSeq2Symbols = _mm512_alignr_epi32(vecSeq2Symbols, vecDummy, 15);
					//printVector(vecSeq1Symbols, "sequence 1: ");
					//printVector(vecSeq2Symbols, "sequence 2: ");

					/*get the substitution score*/
					vecSubScore = _mm512_mask_mov_epi32(vecMismatch, _mm512_mask_cmpeq_epi32_mask(vecMask, vecSeq1Symbols, vecSeq2Symbols), vecMatch);
					vecSubScore = _mm512_mask_mov_epi32(vecSubScore, _mm512_kor(_mm512_cmpeq_epi32_mask(vecSeq1Symbols, vecDummy), _mm512_cmpeq_epi32_mask(vecSeq2Symbols, vecDummy)), vecZero);
					//printVector(vecSubScore, "substitution score\n");

					vecH = _mm512_mask_add_epi32(vecH, vecMask, vecHUp2, vecSubScore);
					vecH = _mm512_mask_max_epi32(vecH, vecMask, _mm512_mask_max_epi32(vecZero, vecMask, vecE, vecF), _mm512_max_epi32(vecH, vecZero));
					/*if(tid == THREAD){
					 printVector(vecHUp2, "vecHUp2: ");
					 printVector(vecH, "Kernel three vecH: ");
					 printVector(vecE, "Kernel three vecE: ");
					 printVector(vecF, "Kernel three vecF: ");
					 }*/

					/*calculate the maximum alignment score*/
					localVecMaxScore = _mm512_max_epi32(localVecMaxScore, vecH);

					/*save the H value for the next diagonal*/
					vecHUp2 = vecHUp;

					/*save the new vertical H and Evalues*/
					_mm512_store_epi32(buffer, vecH);
					_mm512_store_epi32(buffer2, vecE);
					globalHE[startGlobalCol - 15] = int2(buffer[15], buffer2[15]);
					//printf("startGlobalCol: %d Saved E %d H %d\n", startGlobalCol - 15, globalHE[startGlobalCol - 15]._y, globalHE[startGlobalCol - 15]._x);

					/*save the horizontal H and F values*/
					row = (blockCol + 1) * alpha + tileRow + i;
					_mm512_store_epi32(buffer, vecH);
					globalHC[row] = buffer[i];
					_mm512_store_epi32(buffer, vecF);
					globalFC[row] = buffer[i];
				}
			}
		} /*finish an external diagonal*/

		/*save the maximum alignment score*/
		vecMaxScores[tid] = localVecMaxScore;

		/*swap the HD buffers*/
		int32_t* __restrict__ tmp = globalHD;
		globalHD = globalHDC;
		globalHDC = tmp;

		tmp = globalF;
		globalF = globalFC;
		globalFC = tmp;

		tmp = globalH;
		globalH = globalHC;
		globalHC = tmp;

	}

	/*get the maximum alignment score*/
	maxScore = 0;
	for(int32_t i = 0; i < numThreadsUsed; ++i) {
		maxScore = max(maxScore, _mm512_reduce_max_epi32(vecMaxScores[i]));
	}

	/*release memory*/
	_mm_free(vecMaxScores);
	_mm_free(globalHE);
	_mm_free(globalF);
	_mm_free(globalH);
	_mm_free(globalHC);
	_mm_free(globalFC);
	_mm_free(globalHD);
	_mm_free(globalHDC);
#else
	maxScore = -1; /*executed on the host*/
#endif
	return maxScore;
}

