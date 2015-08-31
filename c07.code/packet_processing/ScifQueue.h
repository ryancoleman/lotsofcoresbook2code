/*
 *
 * This file is provided under a dual BSD/GPLv2 license.  When using or
 * redistributing this file, you may do so under either license.
 *
 * GPL LICENSE SUMMARY
 *
 * Copyright(c) 2015 Intel Corporation.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of version 2 of the GNU General Public License as
 * published by the Free Software Foundation.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * BSD LICENSE
 *
 * Copyright(c) 2015 Intel Corporation.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 *  - Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *  - Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in
 *    the documentation and/or other materials provided with the
 *    distribution.
 *  - Neither the name of Intel Corporation nor the names of its
 *    contributors may be used to endorse or promote products derived
 *    from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 */
#ifndef __SCIFQUEUE_H__
#define __SCIFQUEUE_H__

#include <scif.h>

// Write Memory Barrier for serializing writes to MIC.
// Flushes WC buffers too.
#define wmb()  asm volatile ("sfence":::"memory");

const uint32_t PAGE_SIZE = 0x1000;
const uint32_t START_PAGE_OFFSET = 0x80000;
const uint32_t START_BUF_OFFSET = START_PAGE_OFFSET + 0x10000;

#define MAX_QUEUE_NAME 64
#define SCIFQ_MAX_MSG_SIZE 64

// This assumes that one slot will always be unused. So in someone specifies a queue of size 1 (-n),
// we will create a size of 2.
#define QUEUE_FULL(_head, _tail, _size) ((((_head) + 1) % (_size)) == _tail)
#define QUEUE_EMPTY(_head, _tail, _size) ((_head) == (_tail))

void microSleep(uint32_t micros);

// 
// General abstraction for a Queue (FIFO).
// Consits of a virtually contiguous buffer + a head and a tail. 
// In the classic "producer/consumer" example, elements are "enqueued" at head and it is 
// incremented by the producer. The consumer "dequeues" elements from tail and increments
// tail. Tail catches upto head. (think of a queue at the bank).
//
// There are two parameters of interest in this queue, the "capacity" of the queue and the size of each
// element in the queue defined my "msglen"
//
// buf  q_begin                     q_end   buf_end
// |      |                           |       |
// V      V                           V       V
// --------------------------------------------
// |      |      |       |      |     |       |
// | hdr  | elem |  elem |      |elem |  pad  |
// |      |      |       |      |     |       |
// --------------------------------------------
//
// There is also an additional "page" that is used to house the pointers that the remote side can used to 
// refresh its local copy (or write to for "head updates")
//

class ScifQueue 
{
	public:
		virtual bool pair(const char *partner) = 0;	
		virtual void unpair(const char *partner) = 0;
		
		inline uint32_t *get_pEQIndex() const { return pEQIndex; }
		inline uint32_t *get_pDQIndex() const { return pDQIndex; }
		inline uint32_t *get_pDQCount() const { return pDQCount; }
		
		inline uint32_t capacity() const { return Qsize; }
		inline uint32_t msg_length() const { return msglen; }

	protected:
		ScifQueue(const char *qname, uint32_t m_len, uint32_t q_size, 
					     scif_epd_t endpoint, uint32_t fixed = 0); 
		virtual ~ScifQueue();
		
		void barrier() const;				

		// Queue name and partner queue name
		char name[MAX_QUEUE_NAME];
		char partner[MAX_QUEUE_NAME];

		scif_epd_t epd;	
		uint32_t map_fixed;
	
		// size of each queue element (msglen), the total number of elements 
		// in the queue (Qsize) and the size of the buf allocated.
		uint32_t msglen;
		uint32_t Qsize;

		// the above two parameters define the "buf_size" and total number of 
		// PAGE_SIZE pages needed for the queue data-structure..
		uint32_t buf_size;

		// Buffer allocated for the actual queue (begin and end)
		char *buf;
		char *buf_end;

		// This is where the actual *queue* begins and ends (keeping capacity in mind).
		char *q_begin;
		char *q_end;

		// for dequeue counter
		char *page;

		// index of last push (enqueue index)
		uint32_t *pEQIndex;	

		// index of last pop (dequeue index)
		uint32_t *pDQIndex;	

		// total count of dequeued elements.
		uint32_t *pDQCount;	

		// "shadow" head/tail
		uint32_t EQIndex;
		uint32_t DQIndex;	

	private:
		ScifQueue(const ScifQueue&);
		ScifQueue& operator=(const ScifQueue&);
};

class ScifSendQueue : public ScifQueue
{
	public:
		ScifSendQueue(const char *qname, uint32_t m_len, uint32_t q_size, 
		                                 scif_epd_t endpoint, uint32_t fixed = 0);
		~ScifSendQueue();

		bool pair(const char *recv_qname);
		void unpair(const char *recv_qname);

		// push "n" items into the queue (the enqueue operation)
		uint32_t push(const char *data, uint32_t data_len, uint32_t n = 1); 

	private:
		// this is where things are enqueued to
		char *buf_mmap_addr;
		char *page_mmap_addr;
};

class ScifRecvQueue : public ScifQueue
{
	public:
		ScifRecvQueue(const char *qname, uint32_t m_len, uint32_t q_size, 
		                                 scif_epd_t endpoint, uint32_t fixed = 0);
		~ScifRecvQueue();

		bool pair(const char *send_qname);	
		void unpair(const char *send_qname);		
	
		// pop "n" items out of the queue (the dequeue operation).
		uint32_t pop(char *data, uint32_t data_len, uint32_t n = 1);

	private:
		// this is where things are dequeued from
		off_t buf_offset;
		off_t page_offset;
};

#endif
