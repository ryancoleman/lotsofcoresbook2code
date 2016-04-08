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

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <unistd.h>
#include <fcntl.h>
#include <string.h>
#include <sys/mman.h>
#include <sys/ioctl.h>
#include <scif.h>
#include <iostream>
#include "ScifNode.h"
#include "ScifQueue.h"
#include <sys/io.h>

using namespace std;

extern char *optarg;
extern int optind, opterr, optopt;

void r3_cli(void)
{
	// disable interrupts from ring3 
	iopl(3);
	__asm__ __volatile__ ("cli":::"memory");
}
 
void r3_sti(void)
{
	// enable interrupts from ring3 
	__asm__ __volatile__ ("sti":::"memory");
	iopl(0);
}

void usage(char *p)
{
	cerr << p << " -s <message size> -c <queue capacity> -n <total messages>" << endl;
	cerr << "Arguments: " << endl;
	cerr << "\t-s message size - size of each message sent between host and coprocessor " << endl;
	cerr << "\t-c queue capacity - number of elements the queue should be sized to " << endl;
	cerr << "\t-n total messages - total number of messages that should be sent from the host " << endl;
}

int main(int argc, char *argv[])
{
	int err;
	uint32_t msg_size, capacity, num_msg;
	cpu_set_t cpuset;
	char *data;

	// Optimization #6: Thread affinity. 
	CPU_ZERO(&cpuset);
	CPU_SET(42, &cpuset);
	if (sched_setaffinity(0, sizeof(cpuset), &cpuset) < 0) {
		cerr << "failed to set affinity of main thread" << endl;
		exit(1);
	}

	if (argc != 7) {
		usage(argv[0]);
		exit(1);
	}

	while ((err = getopt(argc, argv, "s:c:n:")) != -1) {
		switch(err) {
			case 's':
				msg_size = atoi(optarg);
				break;
			case 'c':
				capacity = atoi(optarg);
				break;
			case 'n':
				num_msg = atoi(optarg);
				break;
			case '?':
				usage(argv[0]);
				exit(1);
		}
	}
	
	cout << "Coprocessor setup as server" << endl;

	// I'm MIC0!
	ScifServer Server("mic0", SCIF_DEVICE_MIC);

	// accepting connections from "host"
	const uint32_t local_port = 2050;
	Server.accept("host", local_port, 2);
 
	ScifRecvQueue RecvQueue("mic_recv_queue", (uint32_t)msg_size, (uint32_t)capacity, Server.pop_epd());
	RecvQueue.pair("host_send_queue");

	ScifSendQueue SendQueue("mic_send_queue", (uint32_t)msg_size, (uint32_t)capacity, Server.pop_epd());
	SendQueue.pair("host_recv_queue");

	uint32_t num_page = (capacity * msg_size)/PAGE_SIZE; 
	
	if ((err = posix_memalign((void **)&data, 0x1000, num_page * PAGE_SIZE)) < 0) {
		cerr << "posix_memalign failed for data" << endl;
		exit(1);
	}

	uint32_t tot_recv = 0;
	uint32_t num_recv = 0, num_send;

	// While there are more messages to get
	while (tot_recv < num_msg) {
		num_send = 0;
		num_recv = RecvQueue.pop((char*)data, msg_size, 1);
		if (num_recv) {
			while (!num_send)  
				num_send = SendQueue.push((char*)data, msg_size, num_recv);
			tot_recv += num_recv;
		}
	}

	cout << "Coprocessor processed  " << tot_recv << " messages! " << endl;	

	uint32_t *p = SendQueue.get_pDQCount();

#ifdef VERBOSE_DEBUG
	cout << "Waiting for Coprocessor to send everything " << endl;
	cout << "Currently DQCount is " << *p << endl;
#endif
	// wait until the host consumes all packets
	while (*p != (uint32_t)num_msg);

	return 0;
}
