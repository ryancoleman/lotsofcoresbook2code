/*
 * AlignCore.cpp
 *
 *  Created on: Jun 21, 2013
 *      Author: yongchao
 */

#include "AlignCore.h"

AlignCore::AlignCore(Align* align, int32_t micIndex) {
	_align = align;

	/*scoring matrix*/
	_match = align->getMatchScore();
	_mismatch = -align->getMismatchPenalty(); /*negative*/
	_gapOE = align->getGapOE(); /*positive*/
	_gapExtend = align->getGapExtend(); /*positive*/

	/*Xeon Phi coprocessor information*/
#ifdef MPI_PARALLEL
	MPI_Comm_size(MPI_COMM_WORLD, &_numProcs);
	MPI_Comm_rank(MPI_COMM_WORLD, &_rank);
#else
	_rank = 0;
	_numProcs = 1;
#endif
	_micIndex = micIndex;
	_numMicThreads = align->getNumMicThreads();

	/*tile sizes on each Xeon Phi*/
	_alpha = align->getVerticalTileSize();
	_beta = align->getHorizonalTileSize();
#ifdef MPI_PARALLEL
	/*tile size for mutliple Xeon Phis*/
	_mpiAlpha = align->getMpiVerticalTileSize();
#endif

	/*get the database information*/
	_seq1 = _align->getSeq(0);
	_seq1Length = _align->getSeqLength(0);
	_seq2 = _align->getSeq(1);
	_seq2Length = _align->getSeqLength(1);
}
