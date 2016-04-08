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
// A Queue is a FIFO with a head and a tail.
// Each queue consists of a fixed capacity (number of elements or "Qsize") where each element
// is of a fixed size.
//

#include "ScifQueue.h"

#define MAX_QUEUE_NAME 64
#include <stdio.h>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <sys/mman.h>
#include <sys/time.h>

using namespace std;

#define EQINDEX_SANITY 0xC0DE0001
#define DQINDEX_SANITY 0xC0DE0002
#define DQCOUNT_SANITY 0xC0DE0003

// Read Timestamp Counter (rdtsc)
static __inline uint64_t rdtsc(void)
{
	uint32_t low, high;
	__asm __volatile("rdtsc" : "=a" (low), "=d" (high));
	return (low | ((uint64_t)high << 32));
}

// t2 - t1, assumes t2 is newer (i.e. there is no check for it)
long long difference(const timeval *t1, const timeval *t2)
{
  long long ret = 0;
  if (t1 != NULL) {
    ret = t1->tv_usec + (1000000LL * t1->tv_sec);
  }

  if (t2 != NULL) {
    ret -= (t2->tv_usec + (1000000LL * t2->tv_sec));
  }

  return ret;
}

// sleep for a bit. The loop method avoids yielding the CPU and a
// scheduling latency when the thread is ready to run eventually.
void microSleep (uint32_t micros)
{
#if 1
   timeval now,then;
   gettimeofday(&now,0);then=now;
   while (difference(&now,&then) <= micros)
   {
      gettimeofday(&now,0);
   }
#else
   usleep(micros);
#endif
}

//-----------------------------------------------------------------------
//           			SCIF Queue
//-----------------------------------------------------------------------

//
// Instantiate a SCIF Queue (made to be accessible from the peer endpoint (can be on 
// the same node or different nodes).
//
// Allocate page sized (and aligned) buffer. But wrap around exactly where the caller 
// wants to stop (Qsize)
//
// See the header file for details on what the queue structure looks like.

ScifQueue::ScifQueue(const char *qname, uint32_t m_len, uint32_t q_size, scif_epd_t endpoint, uint32_t fixed)
:epd(endpoint), map_fixed(fixed), msglen(m_len), Qsize(q_size) 
{
	cout << "ScifQueue::ScifQueue Creating Queue ... " << endl;

	// me
	memset(name, 0, MAX_QUEUE_NAME);
	strncpy(name, qname, MAX_QUEUE_NAME - 1); 

	// will be available at "pair" time.
	memset(partner, 0, MAX_QUEUE_NAME);

	// see QUEUE_FULL for why we bump up the Qsize by 1
	Qsize++;

	// figure out how many pages of memory we will need to size the queue to "Qsize"	
	// The buffer itself starts with a "header" of type long long. 

	// the extra uint64_t will hold "pEQIndex"
	uint32_t size = sizeof(uint64_t) + Qsize * msglen;
	uint32_t num_page = size/PAGE_SIZE;

	// page aligned.
	if (size % PAGE_SIZE) ++num_page;

	// this is how big our buffer needs to be. We will remember the real buffer for later. 
	buf_size = num_page * PAGE_SIZE;

	int err;
	if ((err = posix_memalign((void **)&buf, PAGE_SIZE, buf_size)) < 0) {
		cout << "ScifQueue::ScifQueue memory allocation failure for queue buffer" << endl;
		// error handling left out for simplicity
	}	

	// start with a zeroed queue.
	memset(buf, 0x0, buf_size);
	
	// A separate page for the "dequeue" counter
	if ((err = posix_memalign((void **)&page, PAGE_SIZE, PAGE_SIZE)) < 0) {
		cout << "ScifQueue::ScifQueue memory allocation failure for dequeue page" << endl;
		// error handling left out for simplicity
	}

	memset(page, 0x0, PAGE_SIZE);
	
	// top of "buf"
	buf_end = buf + buf_size;

	// real start and end of the queue 
	q_begin = buf + sizeof(uint64_t);
	q_end = q_begin + Qsize * msglen;

#ifdef VERBOSE_DEBUG
	cout << "Buffer: begin: " << static_cast<const void*>(buf) << " end: " << static_cast<const void*>(buf_end) << endl;
	cout << "Queue: begin: " << static_cast<const void*>(q_begin) << " end: " << static_cast<const void*>(q_end) << endl;
	cout << "Page: begin: " << static_cast<const void*>(page) << endl;
#endif

	// For now map these to local memory, we will eventually get these to point to remote memory.
	pEQIndex = (uint32_t *)buf; 
	pDQIndex = (uint32_t *)page;
	pDQCount = pDQIndex + sizeof(uint64_t);

	cout << "ScifQueue::ScifQueue Queue instantiated! " << endl;
}

ScifQueue::~ScifQueue()
{
	free(page);
	free(buf);
	cout << "ScifQueue::~ScifQueue" << endl;
}

// A barrier across two SCIF endpoints. 
// Essentially a node sends a message (the name of its partner) to its peer/partner 
// and receives a similar message. Then they do a comparison to make sure they know each other!
void ScifQueue::barrier() const
{
	int err;
	char message[SCIFQ_MAX_MSG_SIZE];

	memset(message, 0, SCIFQ_MAX_MSG_SIZE);

	if ((err = scif_send(epd, (char *)partner, SCIFQ_MAX_MSG_SIZE, 1)) <= 0) { 
		cout << "ScifQueue::barrier send failed" << endl;
		// error handling left out for simplicity
	}
 
	if ((err = scif_recv(epd, message, SCIFQ_MAX_MSG_SIZE, 1)) <= 0) { 
		cout << "ScifQueue::barrier: recv failed" << endl;
		// error handling left out for simplicity
	}
	
	if (strcmp(name, message) != 0) {
		cout << "ScifQueue::barrier names mismatch, name: " << string(name) 
						<< " message: " << string(message) << endl; 
		// error handling left out for simplicity
	}
}

//-----------------------------------------------------------------------
//           			SCIF Send Queue
//-----------------------------------------------------------------------
ScifSendQueue::ScifSendQueue(const char *qname, uint32_t m_len, uint32_t q_size, scif_epd_t endpoint, uint32_t fixed)
: ScifQueue(qname, m_len, q_size, endpoint, fixed)
{
}

ScifSendQueue::~ScifSendQueue()
{
}

//
// Send pair operation. Obtain the SCIF offsets for the queue buffer and the dequeue page via send/recv operations.
// Then mmap() them into the local process address space to get access to them in usermode (ring3).
//
bool ScifSendQueue::pair(const char *recv_qname)
{
	bool ret = true;
	off_t offset;

	strncpy(partner, recv_qname, MAX_QUEUE_NAME - 1);	

	// make sure the peer endpoint is here too.
	barrier();

	// Get SCIF offset (registration) for FIFO (queue) buffer
	scif_recv(epd, &offset, sizeof(offset), 1);

#ifdef VERBOSE_DEBUG
	cout << "ScifSendQueue::pair received offset 0x" << std::hex << offset << " for buf " << std::dec << endl;
#endif

	// Optimization #1:
	// Use of scif_mmap() to obtain a pointer to the remote queue. This allows local
	// code to access the buffer directly from ring3 (user mode). A kernel transition is 
	// completely avoided when reading/writing into the queue and manipulating its 
	// head and tail!
	//
	// Optimization #2:
	// scif_mmap() on the host maps remote memory (PCIe) in the address space of the current
	// process. This memory is mapped Write Combining (WC)
	
#ifdef VERBOSE_DEBUG
	cout << "ScifSendQueue::pair attempting to scif_mmap at address " << static_cast<const void*>(buf) << endl;
#endif
	if ((buf_mmap_addr = (char *)scif_mmap(buf, buf_size, SCIF_PROT_READ|SCIF_PROT_WRITE, SCIF_MAP_FIXED, 
							epd, offset)) == MAP_FAILED) {
		ret = false;
		// error handling left out for simplicity
	}

#ifdef VERBOSE_DEBUG
	cout << "ScifSendQueue::pair remote map @ address (queue) " << static_cast<const void*>(buf_mmap_addr) << endl;
#endif
	// make sure the peer endpoint is here too.
	barrier();	

	// Get SCIF offset (registration) for dequeue page next.
	scif_recv(epd, &offset, sizeof(offset), 1);

#ifdef VERBOSE_DEBUG
	cout << "ScifSendQueue::pair received offset 0x" << std::hex << offset << " for page " << std::dec << endl;

	cout << "ScifSendQueue::pair attempting to scif_mmap @ address (page) " << static_cast<const void*>(page) << endl;
#endif
	// use the offset we just got to scif_mmap() a pointer to the queue on the other side
	// NOTE: we're passing in SCIF_MAP_FIXED (i.e. exactly at this address).
	if ((page_mmap_addr = (char *)scif_mmap(page, PAGE_SIZE, SCIF_PROT_READ|SCIF_PROT_WRITE, SCIF_MAP_FIXED, 
							epd, offset)) == MAP_FAILED) {
		ret = false;
		// error handling left out for simplicity
	}

	// These pointers were pointing to the local buf and page. Reset them to
	// the remote side now.
	pEQIndex = (uint32_t *)buf_mmap_addr;
	pDQIndex = (uint32_t *)page_mmap_addr;
	pDQCount = pDQIndex + sizeof(uint64_t);

	q_begin = buf_mmap_addr + sizeof(uint64_t);

	// Now pEQIndex should point to the remote "buf" and pDQIndex should point to the remote "page"
	// Verify by reading what's there.
	
	// A quick sanity test should verify that they are infact pointing to the right remote pages!
	cout << "ScifSendQueue::pair -- Sanity Checking ... ";
	if ((*pEQIndex != EQINDEX_SANITY) || (*pDQIndex != DQINDEX_SANITY) || (*pDQCount != DQCOUNT_SANITY)) {
		cout << "ScifSendQueue:: Pair unexpected failure " << cout << endl;
		ret = -1;
	}
	cout << "Succeeded" << endl;

#ifdef VERBOSE_DEBUG
	cout << "ScifSendQueue:: Pair " << std::hex << *pEQIndex << " " << *pDQIndex << std::dec << endl;
#endif

	// Initialize local and remote indices to 0 -- the queues are empty!
	*pEQIndex = *pDQIndex = *pDQCount = 0;
	EQIndex = DQIndex = 0;

	// make sure the peer endpoint is here too.
  	barrier();

	if (ret) {
		cout << "ScifSendQueue::pair succeeded with token " << recv_qname << endl;
	}
	else {
		cout << "ScifSendQueue::pair failed with token " << recv_qname << endl;
	}

	return ret;
}

void ScifSendQueue::unpair(const char* recv_qname)
{
	int err;

	// Wait for the other side to get here.
	// Start by unmaping the page tables to remote memory.
	barrier();

	if ((err = scif_munmap(page_mmap_addr, PAGE_SIZE)) < 0) {
		cout << "ScifSendQueue::unpair failed for page" << endl;
		// error handling left out for simplicity
	}

	if ((err = scif_munmap(buf_mmap_addr, buf_size)) < 0) {
		cout << "ScifSendQueue::unpair failed for buf" << endl;
		// error handling left out for simplicity
	}
 
	barrier();

	cout << "ScifSendQueue::unpair waiting for remote side to unregister memory now ... ";

	barrier();
	cout << "done " << endl;

	cout << "ScifSendQueue::unpair unpaired with " << recv_qname << " done" << endl;
}

//
// Success: 
// Failure: returns 0
//
uint32_t ScifSendQueue::push(const char* data, uint32_t data_len, uint32_t n)
{
	// Only handle one message at this time to keep it simple.
	if (n != 1) return 0;

	// Optimization #4: "Shadow pointers". Notice that we only "refresh" our local EQIndex 
	// and DQIndex variables with what's really on the other side of the PCIe bus when 
	// it is necessary. This avoids constantly reading the indices for every push/pop 
	// operation.
	if (QUEUE_FULL(EQIndex, DQIndex, Qsize)) {

		// queue appears to be full, really? Perhaps I have stale data?
		// read across PCIe (case where endpoints are across the bus)
		// and refresh my local indices. In this case I'm really interested
		// in finding out if the remote side has advanced tail by popping
		// some messages.
		DQIndex = *pDQIndex;

		// Is it really full?
		if (QUEUE_FULL(EQIndex, DQIndex, Qsize)) 
			return 0;
		
	}

#ifdef VERBOSE_DEBUG
	cout << "ScifSendQueue::push EQIndex " << EQIndex << " DQIndex " << DQIndex << endl;
#endif

	// we only increment in units of "msglen"
	if (data_len <= msglen) {
		// Optimization #3:
		// Push the data to the other side, then serialize before updating the "head"
		// pointer. This is a posted transaction which is cheap.
		// Optimization #6:
		// Use icc/icpc with the right optimizations to ensure that this memcpy will be vectorized.
		memcpy(q_begin + (EQIndex * msglen) , data, data_len);
#ifdef HOST
		// Optimization #2: Memory mapped WC on the host will require a fence (flush) 
		// via a serialization instruction
		wmb();
#endif
		// Optimization #4: Use of shadown pointers.
		EQIndex = (EQIndex + 1) % Qsize;

		// Optimization #3: 
		// Let the other side know that it has one more message. Note that this is 
		// a write operation across PCIe i.e. a posted transaction which is cheap. 
		*pEQIndex = EQIndex;
#ifdef HOST
		// Optimization #2: Memory mapped WC on the host will require a fence (flush) 
		// via a serialization instruction
		wmb();
#endif

#ifdef VERBOSE_DEBUG
	cout << "ScifSendQueue::push EQIndex " << EQIndex << " DQIndex " << DQIndex << endl;
#endif
		return n;
	}

	return 0; // trying to copy a single data element that is >= msglen
}

//-----------------------------------------------------------------------
//           			SCIF Recv Queue
//-----------------------------------------------------------------------

ScifRecvQueue::ScifRecvQueue(const char *qname, uint32_t m_len, uint32_t q_size, scif_epd_t endpoint, uint32_t fixed)
: ScifQueue(qname, m_len, q_size, endpoint, fixed)
{
	// We will write these values in the buffer/page so that the Sender can 
	// verify this after it has scif_mmaped to these pages from the remote side.

#ifdef VERBOSE_DEBUG
	cout << "ScifRecvQueue::ScifRecvQueue: Setting SANITY values for check on pair operation " << endl;
#endif
	*pEQIndex = EQINDEX_SANITY;
	*pDQIndex = DQINDEX_SANITY;
	*pDQCount = DQCOUNT_SANITY;
}

ScifRecvQueue::~ScifRecvQueue()
{
}

bool ScifRecvQueue::pair(const char *send_qname)
{
	bool ret = true;
	off_t offset = 0;

	strcpy(partner, send_qname);

	// make sure the peer is here too.
	barrier();

	// Optimization #1: 
	// scif_mmap() from the other side requires an offset which can be generated by scif_register().
	// Register local memory and generate an offset, then send it to the other side.
	if ((buf_offset = scif_register(epd, buf, buf_size, offset, SCIF_PROT_READ|SCIF_PROT_WRITE, 
							map_fixed ? SCIF_MAP_FIXED : 0)) < 0) {
		ret = false;
	}

#ifdef VERBOSE_DEBUG
	cout << "ScifRecvQueue::pair sending offset 0x" << std::hex << buf_offset << " for buf " << std::dec << endl;
#endif
	scif_send(epd, &buf_offset, sizeof(buf_offset), 1);

	// make sure the peer is here too
	barrier();
		
	// Optimization #1: 
	// scif_mmap() from the other side requires an offset which can be generated by scif_register().
	// Register local memory and generate an offset, then send it to the other side.
	if ((page_offset = scif_register(epd, page, PAGE_SIZE, offset, SCIF_PROT_READ|SCIF_PROT_WRITE, 
							map_fixed? SCIF_MAP_FIXED : 0)) < 0) {
		ret = false;
	}

#ifdef VERBOSE_DEBUG
	cout << "ScifRecvQueue::pair sending offset 0x" << std::hex << page_offset << " for page " << std::dec << endl;
#endif
	scif_send(epd, &page_offset, sizeof(page_offset), 1);

	// final sync
	barrier();

	if (ret) {
		cout << "ScifRecvQueue::pair done with token " << send_qname << endl;
	}
	else {
		cout << "ScifRecvQueue::pair failure with token " << send_qname << endl;
	}

	return ret;
}

void ScifRecvQueue::unpair(const char* send_qname)
{
	int err;

	barrier();
	
	cout << "ScifRecvQueue::unpair wait for remote side to scif_munmap first " << endl;

	barrier();

	if ((err = scif_unregister(epd, page_offset, PAGE_SIZE)) < 0) {
		cout << "scif_unregister failed for page" << endl;
	}

	if ((err = scif_unregister(epd, buf_offset, buf_size)) < 0) {
		cout << "scif_unregister failed for buf" << endl;
	}

	barrier();

	cout << "ScifRecvQueue::unpair with " << send_qname << " done" << endl;
}

uint32_t ScifRecvQueue::pop(char* data, unsigned data_len, uint32_t n)
{
	if (n != 1) return 0;


	// In the case of the Recv Queue, pEQIndex and pDQIndex should be pointers
	// to local memory and this should be not go across the PCIe bus.
	if (QUEUE_EMPTY(*pEQIndex, *pDQIndex, Qsize)) 
		return 0;

#ifdef VERBOSE_DEBUG
	cout << "ScifRecvQueue::pop  EQIndex " << *pEQIndex << " DQIndex " << *pDQIndex << endl;
#endif
	if (data_len <= msglen) {

		// Optimization #3:
		// Push the data to the other side, then serialize before updating the "head"
		// pointer. This is a posted transaction which is cheap.
		// Optimization #6:
		// Use icc/icpc with the right optimizations to ensure that this memcpy will be vectorized.
		memcpy(data, q_begin + (*pDQIndex * msglen), data_len);
#ifdef HOST
		// Optimization #2: Memory mapped WC on the host will require a fence (flush) 
		// via a serialization instruction
		wmb();
#endif
		// Optimization #4: 
		// The remote side will read this value only when it thinks it has a stale copy.
		// Alternatively, we can push this to the remote side too!
		*pDQIndex = DQIndex = (DQIndex + 1) % Qsize;
		*pDQCount += 1;
#ifdef HOST
		// Optimization #2: Memory mapped WC on the host will require a fence (flush) 
		// via a serialization instruction
		wmb();
#endif
		return 1;
	}

	return 0; 	
}
