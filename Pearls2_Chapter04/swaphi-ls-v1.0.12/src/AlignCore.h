/*
 * AlignCore.h
 *
 *  Created on: Jun 21, 2013
 *      Author: yongchao
 */

#ifndef ALIGNCORE_H_
#define ALIGNCORE_H_
#include "Macros.h"
#include "Utils.h"
#include "Align.h"

struct int2 {
	__attribute__((target(mic))) int2() {
		_x = 0;
		_y = 0;
	}
	__attribute__((target(mic))) int2(int32_t x, int32_t y) {
		_x = x;
		_y = y;
	}
	int32_t _x;
	int32_t _y;
}__attribute__((packed));

class AlignCore {
public:
	AlignCore(Align* align, int32_t micIndex);
	virtual ~AlignCore() {
	}
	;

	/*alignment*/
	virtual void align() {
	}
protected:
	Align *_align;
	int32_t _micIndex;
	int32_t _numMicThreads;

	/*scoring scheme*/
	int32_t _match;
	int32_t _mismatch;
	int32_t _gapOE;
	int32_t _gapExtend;

	/*tile sizes on each Xeon Phi*/
	int32_t _alpha;
	int32_t _beta;
#ifdef MPI_PARALLEL
	/*global tile sizes for MPI*/
	int32_t _mpiAlpha;
	int32_t _mpiBeta;
#endif
	/*process rank and number of processes*/
	int32_t _rank;
	int32_t _numProcs;

	/*two sequences*/
	int8_t* _seq1; /*sequence 1*/
	int32_t _seq1Length; /*length of sequence 1*/
	int8_t* _seq2; /*sequence 2*/
	int32_t _seq2Length; /*length of sequence 2*/

	/*static functions*/
	inline __attribute__((target(mic)))
	double getPhiSysTime() {
		double dtime;
		struct timeval tv;

		gettimeofday(&tv, NULL);

		dtime = (double) tv.tv_sec;
		dtime += (double) (tv.tv_usec) / 1000000.0;
		return dtime;
	}
	inline __attribute__((target(mic)))
	void printVector(const __m512i vec, const char* args = NULL) {
#ifdef __MIC__
		int32_t scores[16];

		/*print out the message*/
		if(args) {
			printf("%s", args);
		}

		/*store the alignment score*/
		_mm512_packstorelo_epi32(scores, vec);
		_mm512_packstorehi_epi32(scores + 16, vec);

		/*frow low to high*/
		for(int i = 0; i < 16; ++i) {
			printf("%d ", scores[i]);
		}
		printf("\n");
#endif

	}
};

#endif /* ALIGNCORE_H_ */
