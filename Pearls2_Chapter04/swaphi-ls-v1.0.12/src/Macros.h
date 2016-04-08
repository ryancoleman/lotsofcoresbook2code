/*
 * Macros.h
 *
 *  Created on: Jun 19, 2013
 *      Author: yongchao
 */

#ifndef MACROS_H_
#define MACROS_H_
#define _LARGEFILE64_SOURCE
#include <sys/types.h>
#include <sys/time.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <stdarg.h>
#ifdef COMPRESSED_INPUT
#include <zlib.h>
#endif
#include <string>
#include <list>
#include <vector>
#include <omp.h>
#include <immintrin.h>

using namespace std;

/*program version*/
#define		SWAPHI_LS_VERSION		"v1.0.12"
#define		SWAPHI_LS_VERSION_DATE		"May 19, 2015"

/*user-defined macros*/
#define		FILE_FORMAT_FASTA		1
#define 	FILE_FORMAT_FASTQ		2

/*for device memroy allocation*/
#define MIC_REUSE	alloc_if(0) free_if(0)
#define MIC_ALLOC	alloc_if(1) free_if(0)
#define MIC_FREE	length(0) alloc_if(0) free_if(1)
#define MIC_ALLOC_FREE alloc_if(1) free_if(1)

/*parallelization*/
#define INTER_TASK_MODEL		0
#define INTRA_TASK_MODEL		1

/*subject sequence alignment for the kernles*/
#define SEQ_LENGTH_ALIGN_SHIFT		3
#define SEQ_LENGTH_ALIGN			(1 << SEQ_LENGTH_ALIGN_SHIFT)
#define SEQ_LENGTH_ALIGN_MASK		((1 << SEQ_LENGTH_ALIGN_SHIFT) - 1)
#define DUMMY_NUCLEOTIDE					5

/*workload*/
#define MAX_BYTES_PER_FETCH		(((uint64_t)128) << 20)
struct WorkloadEntry {
	WorkloadEntry() {
		_first = 0;
		_firstAddrOff = 0;
		_numChunks = 0;
	}
	WorkloadEntry(uint64_t first, uint64_t firstAddrOff, uint64_t numChunks) {
		_first = first;
		_firstAddrOff = firstAddrOff;
		_numChunks = numChunks;
	}

	uint64_t _first;
	uint64_t _firstAddrOff;
	uint64_t _numChunks;
};

#define max(a, b) ((a) > (b) ? (a) : (b))
#define min(a, b) ((a) < (b) ? (a) : (b))
#endif /* MACROS_H_ */
