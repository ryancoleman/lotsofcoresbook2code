/*
 * Align.h
 *
 *  Created on: Jun 19, 2013
 *      Author: yongchao
 */

#ifndef ALIGN_H_
#define ALIGN_H_

#include "Macros.h"
#include "Utils.h"
#include "Sequence.h"
#ifdef MPI_PARALLEL
#include <mpi.h>
#include <set>
#include <map>
#include <iostream>
using namespace std;
class RegistrationInfo
{
public:
	RegistrationInfo(int32_t proc, int32_t num) {
		_procs.insert(proc);
		_numXeonPhis = num;
	}
	void add(int32_t proc, int32_t numXeonPhis) {
		if(_procs.size() > 0 && _numXeonPhis != numXeonPhis) {
			Utils::exit("Xeon Phi devices are inconsitent between processes\n");
		}
		if(_procs.size() == 0) {
			_numXeonPhis = numXeonPhis;
		}
		_procs.insert(proc);
	}
	bool isValid() {
		return _procs.size() <= _numXeonPhis;
	}
	set<int32_t>& getProcs() {return _procs;}
	int32_t getNumXeonPhis() {return _numXeonPhis;}
private:
	set<int32_t> _procs;
	int32_t _numXeonPhis;
};
#endif

class Align {
public:
	Align();
	~Align();

	/*parse the parameters*/
	bool parseParams(int argc, char* argv[]);

	/*run the alignment*/
	bool run();

	inline int32_t getMicIndex() {
		return _micIndex;
	}
	inline int32_t getNumMics() {
		return _numMics;
	}
	inline uint32_t getNumMicThreads() {
		return _numMicThreads;
	}
	inline int32_t getMatchScore() {
		return _match;
	}
	inline int32_t getMismatchPenalty() {
		return _mismatch;
	}
	inline int32_t getGapExtend() {
		return _gapExtend;
	}
	inline int32_t getGapOE() {
		return _gapOpen + _gapExtend;
	}
	inline int32_t getVerticalTileSize() {
		return _tiles[0];
	}
	inline int32_t getHorizonalTileSize() {
		return _tiles[1];
	}
#ifdef MPI_PARALLEL
	inline int32_t getMpiVerticalTileSize() {
		return _mpiAlpha;
	}
#endif
	inline int32_t getTiling() {
		return _tiling;
	}
	inline int32_t getSeqLength(int index) {
		return _seqLengths[index];
	}
	inline int8_t* getSeq(int index) {
		return _seqs[index];
	}
private:
	/*parallelization*/
	int32_t _compMode;
	int32_t _micIndex; /*the index of the Xeon Phi*/
	int32_t _numMics; /*number of Xeon Phis*/
	uint32_t _numMicThreads; /*number of device threads*/
	int32_t _numSimChars; /*number of simulated characters*/
	int32_t _rank;
	int32_t _numProcs;
	int32_t _numProcsSharingHost;
	int32_t _order;

	/*penalty*/
	int32_t _match;
	int32_t _mismatch;
	int32_t _gapOpen;
	int32_t _gapExtend;
	int32_t _gapOE;

	/*tiling*/
	int32_t _tiles[2];
#ifdef MPI_PARALLEL
	int32_t _mpiAlpha; /*vertical tile size for MPI*/
#endif
	int32_t _tiling;

	/*sequence file*/
	string _queryFileName;
	string _subjectFileName;

	/*sequences*/
	int8_t* _seqs[2];
	int32_t _seqLengths[2];

private:
	void _setDefaults(); /*set default parameters*/
	/*print out parameters*/
	void _printUsage();
	/*load sequences*/
	bool _loadSequences();
	/*simulate sequences*/
	bool _simulateSequences();
#ifdef MPI_PARALLEL
	/*assign the device*/
	int32_t _assignXeonPhis();
#endif

	/*check the availabiltiy of device*/
	bool _checkDevice(int32_t index, int32_t &numThreads) {
		bool targetOk = false;
		int32_t nprocs = 0;
#pragma offload target (mic:index)
		{
#ifdef __MIC__
			targetOk = true;
			nprocs = omp_get_num_procs();
#else
			targetOk = false;
			nprocs = 0;
#endif
		}
		//Utils::log("Number of processor cores on device %d: %d\n", index, nprocs);
		numThreads = nprocs;
		return targetOk;
	}
};

#endif /* ALIGN_H_ */
