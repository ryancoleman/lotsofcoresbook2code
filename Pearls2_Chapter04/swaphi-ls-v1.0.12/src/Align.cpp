/*
 * Align.cpp
 *
 *  Created on: Jun 19, 2013
 *      Author: yongchao
 */

#include "Align.h"
#include "SeqFileParser.h"
#include "AlignCoreNaive.h"
#include "AlignCoreTiling.h"
#include <sys/mman.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

void Align::_printUsage() {
#ifdef MPI_PARALLEL
	Utils::log("MPI-SWAPHI-LS (%s) is a parallelized Smith-Waterman algorithm on Intel Xeon Phi clusters\nfor pairwise comparison of long DNA sequences (released on %s)\n", SWAPHI_LS_VERSION, SWAPHI_LS_VERSION_DATE);
	Utils::log("Usage: MPI-SWAPHI-LS [options]\n");
#else
	Utils::log(
			"SWAPHI-LS (%s) is a parallelized Smith-Waterman algorithm on the Intel Xeon Phi\nfor pairwise comparison of long DNA sequences (released on %s)\n",
			SWAPHI_LS_VERSION, SWAPHI_LS_VERSION_DATE);
	Utils::log("Usage: SWAPHI-LS [options]\n");
#endif
	Utils::log("Input:\n");
	Utils::log("\t-i <str> (query DNA sequence file [REQUIRED])\n");
	Utils::log("\t-j <str> (subject DNA sequence file [REQUIRED])\n");
	Utils::log(
			"\t-k <int> (place the longer sequence horizontally, default = %d)\n",
			_order);

	Utils::log("Scoring scheme:\n");
	Utils::log("\t-m <int> (match score, default = %d)\n", _match);
	Utils::log("\t-M <int> (mismatch penalty, default = %d)\n", _mismatch);
	Utils::log("\t-g <int> (gap opening penalty, default = %d)\n", _gapOpen);
	Utils::log("\t-e <int> (gap extension penalty, default = %d)\n",
			_gapExtend);

	Utils::log("Tiling:\n");
	Utils::log(
			"\t-a <int> (the vertical tile size within each device, default = %d [NO need to change])\n",
			_tiles[0]);
	Utils::log(
			"\t-b <int> (the horizonal tile size within each device, default = %d [0 means auto])\n",
			_tiles[1]);
#ifdef MPI_PARALLEL
	Utils::log(
			"\t-C <int> (the vertical block size between devices, default = %d)\n",
			_mpiAlpha);
#else
	Utils::log("\t-c <int> (enable the tiling, default = %d)\n", _tiling);
#endif

	Utils::log("Compute:\n");
	Utils::log(
			"\t-t <int> (number of threads per Xeon Phi, deafult = %d [0 means auto])\n",
			_numMicThreads);
#ifndef MPI_PARALLEL
	Utils::log("\t-x <int> (index of the Xeon Phi used, default = %d)\n",
			_micIndex);
#endif
	Utils::log("\t-n <int> (simulate #int characters, deafult = %d [SPEED TEST])\n",
			_numSimChars);
}

Align::Align() {
	/*set the default vaues*/
	_setDefaults();
}
Align::~Align() {
	/*release memory*/
	for (int i = 0; i < 2; ++i) {
		if (_seqs[i]) {
			_mm_free(_seqs[i]);
		}
	}
}

bool Align::run() {
	int32_t numThreads;
	AlignCore* core;

	/*create the object*/
#ifdef MPI_PARALLEL

	/*check the usability of the device*/
	//Utils::log("process: %d id: %d\n", _rank, _micIndex);
	if (_checkDevice(_micIndex, numThreads) == false) {
		Utils::log("Failed to offload on device %d\n", _micIndex);
		return false;
	}
	/*automatically set the number of device threads used*/
	if(_numMicThreads == 0) {
		_numMicThreads = numThreads;
	}

	/*create the object*/
	core = new AlignCoreTiling(this, _micIndex);
#else

	/*check the usability of the device*/
	if (_checkDevice(_micIndex, numThreads) == false) {
		Utils::log("Failed to offload on device %d\n", _micIndex);
		return false;
	}
	if (_numMicThreads == 0) {
		_numMicThreads = numThreads;
	}

	/*create the object*/
	if (_tiling) {
		core = new AlignCoreTiling(this, _micIndex);
	} else {
		core = new AlignCoreNaive(this, _micIndex);
	}
#endif

	/*run the kernel*/
	core->align();

	/*release the object*/
	delete core;

	return true;
}
bool Align::parseParams(int argc, char* argv[]) {
	int32_t ch;
	int32_t val;

	/*check parameters*/
	if (argc < 2) {
		_printUsage();
		return false;
	}

#ifdef MPI_PARALLEL
	/*get the number of processes*/
	MPI_Comm_size(MPI_COMM_WORLD, &_numProcs);

	/*get the rank of the process*/
	MPI_Comm_rank(MPI_COMM_WORLD, &_rank);
#endif

	/*check the number of Xeon Phis*/
#ifdef __INTEL_OFFLOAD
	_numMics = _Offload_number_of_devices();
#endif
	if (_numMics < 1) {
		Utils::log("Xeon Phi coprocessors are unavailable\n");
		return false;
	}
	/*print out the command line*/
	if(_rank == 0){
		Utils::log("Command line: ");
  	for (int i = 0; i < argc; ++i) {
			Utils::log("%s ", argv[i]);
  	}
		Utils::log("\n");
	}

	/*parse parameters*/
#ifndef MPI_PARALLEL
	while ((ch = getopt(argc, argv, "i:j:m:M:g:e:t:x:n:a:b:t:c:k:")) >= 0) {
#else
		while ((ch = getopt(argc, argv, "i:j:m:M:g:e:t:n:a:b:t:C:c:k:")) >= 0) {
#endif
		switch (ch) {
		case 'i':
			_queryFileName = optarg;
			break;
		case 'j':
			_subjectFileName = optarg;
			break;
		case 'k':
			val = atoi(optarg);
			if (val != 0) {
				val = 1;
			}
			_order = val;
			break;
		case 'm':
			val = atoi(optarg);
			if (val < 0) {
				val = 0;
			}
			_match = val;
			break;
		case 'M':
			val = atoi(optarg);
			if (val < 0) {
				val = 0;
			}
			_mismatch = val;
			break;
		case 'g':
			val = atoi(optarg);
			if (val < 0) {
				val = 0;
			}
			_gapOpen = val;
			break;
		case 'e':
			val = atoi(optarg);
			if (val < 0) {
				val = 0;
			}
			_gapExtend = val;
			break;
		case 'a':
			val = atoi(optarg);
			if (val < 1) {
				val = 1;
			}
			/*must be a multiple of 16*/
			_tiles[0] = (val + 15) >> 4 << 4;
			break;
		case 'b':
			val = atoi(optarg);
			if (val <= 0) {
				val = 0;
			}
			/*must be a multiple of 16*/
			_tiles[1] = (val + 15) >> 4 << 4;
			break;
#ifdef MPI_PARALLEL
			case 'C':
			val = atoi(optarg);
			if (val < 16) {
				val = 16;
			}
			/*must be a multiple of 16*/
			_mpiAlpha = (val + 15) >> 4 << 4;
			break;
#endif
		case 'c':
#ifdef MPI_PARALLEL
		/*do nothing for the MPI program*/
		break;
#else
			val = atoi(optarg);
			if (val != 0) {
				val = 1;
			}
			_tiling = val;
			break;
#endif
		case 't':
			val = atoi(optarg);
			if (val < 0) {
				val = 0;
			}
			_numMicThreads = val;
			break;
#ifndef MPI_PARALLEL
		case 'x':
			val = atoi(optarg);
			if (val >= _numMics) {
				val = _numMics - 1;
			}
			if (val < 0) {
				val = 0;
			}
			_micIndex = val;
			break;
#endif
		case 'n':
			val = atoi(optarg);
			if (val < 0) {
				val = 0;
			}
			_numSimChars = val;
			break;
		default:
			Utils::log("Unsupported option: %s\n", optarg);
			return false;
		}
	}

	/*check the input file*/
	if (_numSimChars == 0 && (_queryFileName.length() == 0 || _subjectFileName.length() == 0)) {
		Utils::log("Two input files must be given through options \"-i\" and \"-j\".\n");
		return false;
	}

#ifdef MPI_PARALLEL
	/*assign the Xeon Phis*/
	_assignXeonPhis();

#endif

	/*read sequences*/
	double stime = Utils::getSysTime();
	if ((_numSimChars ? _simulateSequences() : _loadSequences()) == false) {
		return false;
	}
	double etime = Utils::getSysTime();
	Utils::log("Sequence loading takes %f seconds\n", etime - stime);
	return true;
}

bool Align::_loadSequences() {
	int32_t numSeqs, alignedLength;
	Sequence seqs[2];
	SeqFileParser* parser;
	if(_rank == 0){
		Utils::log("Start loading sequences\n");
	}

	/*load the query sequence*/
	parser = new SeqFileParser(_queryFileName.c_str());
	if (!parser) {
		Utils::log("Failed to open the file: %s\n", _queryFileName.c_str());
		return false;
	}
	/*get the sequennces*/
	if(!parser->getSeq(seqs[0])){
		Utils::log("The query sequence file doest not contain any sequence\n");
		return false;
	}
	/*release the file parser*/
	delete parser;

	/*load the subject sequence*/
	parser = new SeqFileParser(_subjectFileName.c_str());
	if(!parser){
		Utils::log("Failed to open the file: %s\n", _subjectFileName.c_str());
		return false;
	}
	/*get the sequence*/
	if(!parser->getSeq(seqs[1])){
		Utils::log("The subject sequence files does not contain any sequence\n");
		return false;
	}
	/*release the file parser*/
	delete parser;

	/*print out the sequence information*/
	if(_rank == 0){
		Utils::log("Query: %s\n", seqs[0]._name);
		Utils::log("\tlength: %d\n", seqs[0]._length);
		Utils::log("Target: %s\n", seqs[1]._name);
		Utils::log("\tlength: %d\n", seqs[1]._length);
	}

	/*compare the lengths of the sequences*/
	Sequence* seqPtrs[2] = { &seqs[0], &seqs[1] };
	if (_order) {
		/*the longer sequence is placed horizontally*/
		if (seqPtrs[0]->_length > seqPtrs[1]->_length) {
			swap(seqPtrs[0], seqPtrs[1]);
		}
	} else {
		/*the shoter sequence is placed horizontally*/
		if (seqPtrs[0]->_length < seqPtrs[1]->_length) {
			swap(seqPtrs[0], seqPtrs[1]);
		}
	}

	/*for sequence one*/
	for (numSeqs = 0; numSeqs < 2; ++numSeqs) {
		/*aligned to 16*/
		alignedLength =
				_tiling ?
						(seqPtrs[numSeqs]->_length + 15) >> 4 << 4 :
						seqPtrs[numSeqs]->_length;

		/*allocate space*/
		_seqs[numSeqs] = (int8_t*) _mm_malloc(alignedLength, 64); /*aligned to 64 bytes*/

		/*store the sequence length*/
		_seqLengths[numSeqs] = alignedLength;

		/*store the symbols*/
		memcpy(_seqs[numSeqs], seqPtrs[numSeqs]->_bases,
				seqPtrs[numSeqs]->_length);

		/*padding the dummy symbols*/
		for (int32_t i = seqPtrs[numSeqs]->_length; i < alignedLength; ++i) {
			_seqs[numSeqs][i] = DUMMY_NUCLEOTIDE;
		}

		/*reverse the 2nd sequence if tiling is not used*/
		if (!_tiling && numSeqs == 1) {
			for (int32_t i = 0; i < alignedLength / 2; ++i) {
				swap(_seqs[numSeqs][i], _seqs[numSeqs][alignedLength - 1 - i]);
			}
		}
	}
	if(_rank == 0){
		Utils::log("Finish loading sequences\n");
	}
	return true;
}
bool Align::_simulateSequences() {
	int32_t alignedLength;

	Utils::log("Start simulating sequences\n");
	/*simulate sequence one*/
	alignedLength = _tiling ? (_numSimChars + 15) >> 4 << 4 : _numSimChars;

	_seqs[0] = (int8_t*) _mm_malloc(alignedLength, 64);
	Utils::log("Length of sequence 1 (placed vertically): %d\n", alignedLength);

	/*simulate the sequence*/
	srand48(11);
	for (int i = 0; i < alignedLength; ++i) {
		_seqs[0][i] = (int8_t) (drand48() * 4);
	}
	_seqLengths[0] = alignedLength;

	/*simulate sequence two*/
	alignedLength = _tiling ? (_numSimChars + 15) >> 4 << 4 : _numSimChars;
	Utils::log("Length of sequence 2 (placed horizontally): %d\n",
			alignedLength);
	_seqs[1] = (int8_t*) _mm_malloc(alignedLength, 64);

	/*simulate the sequence*/
	for (int i = 0; i < alignedLength; ++i) {
		_seqs[1][i] = (int8_t) (drand48() * 4);
	}
	_seqLengths[1] = alignedLength;

	/*reverse the sequence if tiling is not used*/
	if (!_tiling) {
		for (int32_t i = 0; i < alignedLength / 2; ++i) {
			swap(_seqs[1][i], _seqs[1][alignedLength - 1 - i]);
		}
	}

	Utils::log("Finish simulating sequences\n");
	return true;
}

#ifdef MPI_PARALLEL
int32_t Align::_assignXeonPhis()
{
	const int32_t hostNameLength = 4095;
	const int32_t maxBuffer = 4095;
	char hostName[4096];
	char buffer[4096];
	int32_t numXeonPhis;
	int32_t totalXeonPhis = 0;

	//get the host name
	if(gethostname(hostName, hostNameLength)) {
		fprintf(stderr, "Get host name failed");
	}

	if(_rank == 0) {
		map<string, RegistrationInfo*> hostMap;

		/*insert the Xeon Phi information on its own node*/
		string host((const char*)hostName);
		map<string, RegistrationInfo*>::iterator iter = hostMap.find(host);
		if(iter != hostMap.end()) {
			iter->second->add(_rank, _numMics);
		} else {
			RegistrationInfo* regInfo = new RegistrationInfo(_rank, _numMics);
			hostMap.insert(pair<string, RegistrationInfo*>(host, regInfo));
		}
		/*receive the Xeon Phi information from any other process*/
		for(int32_t rank = 1; rank < _numProcs; ++rank) {
			MPI_Status status;
			MPI_Recv(buffer, maxBuffer, MPI_BYTE, rank, 0, MPI_COMM_WORLD, &status);

			/*forming a string*/
			int bufferLength;
			MPI_Get_count(&status, MPI_BYTE, &bufferLength);
			buffer[bufferLength] = '\0';

			/*get the hostname and number of devices*/
			sscanf(buffer, "%s %d", hostName, &numXeonPhis);

			/*insert the Xeon Phi information*/
			string host((const char*)hostName);
			map<string, RegistrationInfo*>::iterator iter = hostMap.find(host);
			if(iter != hostMap.end()) {
				iter->second->add(rank, numXeonPhis);
			} else {
				RegistrationInfo* regInfo = new RegistrationInfo(rank, numXeonPhis);
				hostMap.insert(pair<string, RegistrationInfo*>(host, regInfo));
			}
		}
		//check the validity of the each host, and assign Xeon Phis to each process
		for(map<string, RegistrationInfo*>::iterator iter = hostMap.begin(); iter != hostMap.end(); ++iter) {
			if(!iter->second->isValid()) {
				Utils::log("The number of processes exceed the number of Xeon Phis in the node %s\n", iter->first.c_str());
				Utils::log("Processes: %ld devices: %d\n", iter->second->getProcs().size(), iter->second->getNumXeonPhis());
				exit(-1);
			}
			int32_t phiIndex = 0;
			string host = iter->first;
			set<int32_t> procs = iter->second->getProcs();
			for(set<int32_t>::iterator piter = procs.begin(); piter != procs.end(); ++piter) {
				if(*piter == 0) {
					//set the device ID
					_micIndex = phiIndex;
					_numProcsSharingHost = procs.size();
				} else {
					//send the device ID
					MPI_Send(&phiIndex, 1, MPI_INTEGER, *piter, 0, MPI_COMM_WORLD);

					//send the number of procs in the host
					int32_t procNum = procs.size();
					MPI_Send(&procNum, 1, MPI_INTEGER, *piter, 0, MPI_COMM_WORLD);
				}
				Utils::log("Process %d uses the Xeon Phi with ID %d in the host \"%s\"\n", *piter, phiIndex, host.c_str());
				++phiIndex;
			}
			totalXeonPhis += procs.size();
		}

		//release the host map
		for(map<string, RegistrationInfo*>::iterator iter = hostMap.begin(); iter != hostMap.end(); ) {
			delete iter->second;
			hostMap.erase(iter++);
		}
	} else {
		int32_t phiIndex;
		//send the Xeon Phi information in this host
		sprintf(buffer, "%s %d", hostName, _numMics);
		MPI_Send(buffer, strlen(buffer), MPI_BYTE, 0, 0, MPI_COMM_WORLD);

		//receive the assigned device ID
		MPI_Status status;
		MPI_Recv(&phiIndex, 1, MPI_INTEGER, 0, 0, MPI_COMM_WORLD, &status);
		//set the GPU device number
		_micIndex = phiIndex;

		//receive the number of procs in the host
		MPI_Recv(&_numProcsSharingHost, 1, MPI_INTEGER, 0, 0, MPI_COMM_WORLD, &status);
		totalXeonPhis = 1;
	}

	//broadcast the number of GPUs in the distribued system
	MPI_Bcast(&totalXeonPhis, 1, MPI_INTEGER, 0, MPI_COMM_WORLD);

	return totalXeonPhis;
}
#endif

void Align::_setDefaults() {
	/*set scoring matrix*/
	_match = 1;
	_mismatch = 3;
	_gapOpen = 5;
	_gapExtend = 2;

	/*number of threads per Xeon Phi*/
	_micIndex = 0;
	_numMics = 1;
	_numMicThreads = 0;
	_numSimChars = 0;
	_order = 1;

	/*tiling*/
	_tiles[0] = 16; /*vertical tile size*/
	_tiles[1] = 0; /*horizontal tile size*/
	_tiling = 1;
#ifdef MPI_PARALLEL
	_mpiAlpha = 1 << 17;
#endif
	_rank = 0;
	_numProcs = 1;

	/*sequences*/
	_seqs[0] = _seqs[1] = NULL;
	_seqLengths[0] = _seqLengths[1] = 0;
}
