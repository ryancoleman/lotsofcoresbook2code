#include "SeqFileParser.h"
#include "Utils.h"
#include <zlib.h>

SeqFileParser::SeqFileParser(const char* path, size_t BUFFER_SIZE) {

	//allocate buffer for file reading
	_size = BUFFER_SIZE;
	_length = 0;
	_buffer = new uint8_t[_size + 1];
	if (_buffer == NULL) {
		Utils::exit("Memory allocation failed in file %s in line %d\n",
				__FUNCTION__, __LINE__);
	}

	/*create the file buffer*/
	_fileBufferSentinel = 0;
	_fileBufferLength = 0;
	_fileBufferR = new uint8_t[4096 + 8];
	if (_fileBufferR == NULL) {
		Utils::exit("Memory allocation failed in file %s in line %d\n",
				__FUNCTION__, __LINE__);
	}
	_fileBuffer = _fileBufferR + 8; /*make it aligned*/

	/*open the input file*/
	if (strcmp(path, "-") == 0) {
		_fp = myopenstdin("rb");
	} else {
		_fp = myfopen(path, "rb");
	}
	if (_fp == NULL) {
		Utils::cmpExit("Failed to open file: %s\n", path);
	}

	//detecting the file format in the first line
	int ch;
	while ((ch = myfgetc(_fp)) != -1 && ch != '>' && ch != '@' && ch != '\n')
		;
	if (ch == -1 || ch == '\n') {
		Utils::exit("Unrecognized file format\n");
	} else if (ch == '>') {
		_format = FILE_FORMAT_FASTA;
		myungetc(ch, _fp);
		//Utils::log("FASTA format identified\n");
	} else {
		_format = FILE_FORMAT_FASTQ;
		myungetc(ch, _fp);
		//Utils::log("FASTQ format identified\n");
	}
}
SeqFileParser::~SeqFileParser() {
	if (_buffer) {
		delete[] _buffer;
	}

	if (_fileBufferR) {
		delete[] _fileBufferR;
	}
	/*close the file*/
	myfclose(_fp);
}

uint32_t SeqFileParser::getFastaSeq(Sequence& seq) {
	int ch;

	//find the header
	while ((ch = myfgetc(_fp)) != -1 && ch != '>')
		;
	if (ch == -1) {
		return 0; //reach the end of file
	}
	//read the sequence name (only one line)
	_length = 0;
	while ((ch = myfgetc(_fp)) != -1 && ch != '\n') {
		if (_length >= _size) {
			resizeBuffer(_size + 256);
		}
		if (isspace(ch)) {
			ch = '\0';
		}
		_buffer[_length++] = ch;
	}
	if (ch == -1) {
		Utils::exit("Incomplete file\n");
	}
	_buffer[_length] = '\0';

	/*trim characters /[12]$ like BWA*/
	if (_length > 2 && _buffer[_length - 2] == '/'
			&& (_buffer[_length - 1] == '1' || _buffer[_length - 1] == '2')) {
		_length -= 2;
		_buffer[_length] = '\0';
	}

	//save the sequence name
	seq.setNameSize(_length + 1); /*adjust the name buffer size*/
	strcpy((char*) seq._name, (char*) _buffer);

	//read the sequence bases
	_length = 0;
	do {
		//filter out the blank lines
		while ((ch = myfgetc(_fp)) != -1 && (ch == '\r' || ch == '\n'))
			;
		if (ch == -1)
			break; //reach the end of file
		if (ch == '>') { //reaching another sequence
			myungetc(ch, _fp);
			break;
		}
		//save the current encoded base
		if (_length >= _size) {
			resizeBuffer(_size + 256);
		}
		_buffer[_length++] = _encode(ch);
	} while (1);

	//save the sequence length and its bases
	seq._length = _length;
	if (_length > 0) {
		/*adjust the sequence buffer size*/
		seq.setSequenceSize(seq._length);
		memcpy(seq._bases, _buffer, _length);
	}

	return _length;
}

uint32_t SeqFileParser::getFastqSeq(Sequence& seq) {
	int ch;

	//find the header
	while ((ch = myfgetc(_fp)) != -1 && ch != '@')
		;
	if (ch == -1)
		return 0; //reach the end of file

	//read the sequence name (only one line)
	_length = 0;
	while ((ch = myfgetc(_fp)) != -1 && ch != '\n') {
		if (_length >= _size) {
			resizeBuffer(_size + 256);
		}
		if (isspace(ch)) {
			ch = '\0';
		}
		_buffer[_length++] = ch;
	}
	if (ch == -1) {
		Utils::exit("Incomplete file\n");
	}
	_buffer[_length] = '\0';

	/*trim characters /[12]$ like BWA*/
	if (_length > 2 && _buffer[_length - 2] == '/'
			&& (_buffer[_length - 1] == '1' || _buffer[_length - 1] == '2')) {
		_length -= 2;
		_buffer[_length] = '\0';
	}

	//save the sequence name
	seq.setNameSize(_length + 1); /*adjust the name buffer size*/
	strcpy((char*) seq._name, (char*) _buffer);

	//read the sequence bases
	_length = 0;
	do {
		//filter out the blank lines
		while ((ch = myfgetc(_fp)) != -1 && (ch == '\r' || ch == '\n'))
			;
		if (ch == -1)
			Utils::exit("Incomplete FASTQ file\n");
		if (ch == '+')
			break; //the comment line

		//save the current encoded base
		if (_length >= _size) {
			resizeBuffer(_size + 256);
		}
		_buffer[_length++] = _encode(ch);

	} while (1);

	//save the sequence length and its bases
	seq._length = _length;

	/*adjust the sequence buffer size*/
	seq.setSequenceSize(seq._length);
	memcpy(seq._bases, _buffer, _length);

	//read the comment line (only one line)
	while ((ch = myfgetc(_fp)) != -1 && ch != '\n')
		;
	if (ch == -1)
		Utils::exit("Incomplete FASTQ file\n");

	//read the quality scores
	while ((ch = myfgetc(_fp)) != -1 && ch != '\n')
		;

	return seq._length;
}

void SeqFileParser::resizeBuffer(uint32_t nsize) {

	/*check the buffer size*/
	if (nsize <= _size) {
		return;
	}

	//allocate a new buffer
	_size = nsize * 2;	/*double the buffer size*/
	uint8_t* nbuffer = new uint8_t[_size];
	if (!nbuffer) {
		Utils::exit("Memory reallocation failed in file %s in line %d\n",
				__FUNCTION__, __LINE__);
	}
	//copy the old data
	memcpy(nbuffer, _buffer, _length);

	//release the old buffer
	delete[] _buffer;
	_buffer = nbuffer;
}

