/*
 * Sequence.cpp
 *
 *  Created on: Dec 23, 2011
 *      Author: yongchao
 */

#include "Sequence.h"
#include "SeqFileParser.h"
#include "Utils.h"

Sequence::Sequence() {
	_name = NULL;
	_bases = NULL;
	_length = 0;
}
Sequence::Sequence(const Sequence & s) {
	_length = s._length;
	if (_length == 0) {
		_name = NULL;
		_bases = NULL;
		return;
	}
	if (s._name) {
		_name = new char[strlen(s._name) + 1];
		if (_name == NULL) {
			Utils::exit("Memory allocation failed in function %s line %d\n",
					__FUNCTION__, __LINE__);
		}
		strcpy(_name, s._name);
	}
	if (s._bases) {
		_bases = (uint8_t*) _mm_malloc(_length, 64);
		if (_bases == NULL) {
			Utils::exit("Memory allocation failed in function %s line %d\n",
					__FUNCTION__, __LINE__);
		}
		memcpy(_bases, s._bases, _length);
	}
}
Sequence::~Sequence() {
	clear();
}
void Sequence::clear() {
	if (_name) {
		delete[] _name;
	}
	if (_bases) {
		_mm_free(_bases);
	}
	_name = NULL;
	_bases = NULL;
	_length = 0;
}
void Sequence::setNameSize(uint32_t size) {
	_name = new char[size];
	if (_name == NULL) {
		Utils::exit("Memory allocation failed in function %s line %d\n",
				__FUNCTION__, __LINE__);
	}
}
void Sequence::setSequenceSize(uint32_t size) {
	_bases = (uint8_t*) _mm_malloc(size, 64);
	if (_bases == NULL) {
		Utils::exit("Memory allocation failed in function %s line %d\n",
				__FUNCTION__, __LINE__);
	}
}
