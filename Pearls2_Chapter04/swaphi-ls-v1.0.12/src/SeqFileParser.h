/*
 * SeqFileParser.h
 *
 *  Created on: Dec 23, 2011
 *      Author: yongchao
 */

#ifndef SEQFILEPARSER_H_
#define SEQFILEPARSER_H_
#include "Macros.h"
#include "Sequence.h"
#include "Utils.h"
#include "MyFile.h"

class SeqFileParser {
public:
	/*public member functions*/
	SeqFileParser(const char* path, size_t BUFFER_SIZE = 4095);
	~SeqFileParser();

	//get the next sequence from the file
	inline uint32_t getSeq(Sequence& seq) {
		return (_format == FILE_FORMAT_FASTA) ?
				getFastaSeq(seq) : getFastqSeq(seq);
	}

private:
	void resizeBuffer(uint32_t nsize);
	uint32_t getFastaSeq(Sequence& seq);
	uint32_t getFastqSeq(Sequence& seq);

	/*buffered file operations*/
	inline int myfgetc(MyFilePt file) {
		/*check the end-of-file*/
		if (_fileBufferSentinel >= _fileBufferLength) {
			/*re-fill the buffer*/
			_fileBufferSentinel = 0;
			/*read file*/
			_fileBufferLength = myfread(_fileBuffer, 1, 4096, file);
			if (_fileBufferLength == 0) {
				/*reach the end of the file*/
				if (myfeof(file)) {
					return -1;
				} else {
					Utils::exit("File reading failed in function %s line %d\n",
							__FUNCTION__, __LINE__);
				}
			}
		}
		/*return the current character, and increase the sentinel position*/
		return _fileBuffer[_fileBufferSentinel++];
	}
	inline int myungetc(int ch, MyFilePt file) {
		if (_fileBufferSentinel >= 0) {
			_fileBuffer[--_fileBufferSentinel] = ch;
		} else {
			Utils::log("Two consecutive ungetc operations occurred\n");
			return -1; /*an error occurred, return end-of-file marker*/
		}
		return ch;
	}
	inline uint8_t _encode(char ch) {
		switch (ch) {
		case 'a':
		case 'A':
			return 0;
		case 'c':
		case 'C':
			return 1;
		case 'g':
		case 'G':
			return 2;
		case 't':
		case 'T':
			return 3;
		}
		return 4;
	}

	//buffer for file reading
	uint8_t* _buffer;
	uint32_t _length;
	uint32_t _size;

	//FASTA/FASTQ file handler
	MyFilePt _fp;
	uint8_t* _fileBufferR;
	uint8_t* _fileBuffer;
	int _fileBufferLength;
	int _fileBufferSentinel;
	int _format;

	friend class Sequence;
};

#endif /* SEQFILEPARSER_H_ */
