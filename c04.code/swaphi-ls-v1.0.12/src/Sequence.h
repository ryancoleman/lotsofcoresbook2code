/*
 * Sequence.h
 *
 *  Created on: Dec 23, 2011
 *      Author: yongchao
 */

#ifndef SEQUENCE_H_
#define SEQUENCE_H_
#include "Macros.h"
#include "Utils.h"

struct Sequence {
	/*member functions*/
	Sequence();
	Sequence(const Sequence& s);
	~Sequence();
	void setNameSize(uint32_t size);
	void setSequenceSize(uint32_t size);
	void clear();

	/*member variables*/
	char* _name;
	uint8_t* _bases;
	uint32_t _length;
};

#endif /* SEQUENCE_H_ */
