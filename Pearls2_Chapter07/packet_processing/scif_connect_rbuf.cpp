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
#include <sys/time.h>
#include <scif.h>
#include <iostream>
#include <vector>
#include <numeric>
#include <algorithm>
#include <fstream>

#include "ScifNode.h"
#include "ScifQueue.h"

#define MAX_N 100000000

using namespace std;
extern char *optarg;
extern int optind, opterr, optopt;


// track time for send/recv of data to the coprocessor
uint64_t sendTime[MAX_N] = {0};
uint64_t recvTime[MAX_N] = {0};
uint64_t capacity;

//uint64_t* recvTime1, *recvTime2, *recvTime3;
uint64_t *micRDTSC1, *micRDTSC2;
uint32_t num_msg, msg_size;

static inline uint64_t get_time_stamp()
{
	uint64_t time;
	timeval tv;

	gettimeofday(&tv, NULL);
	
	// convert to micro-seconds
	time = tv.tv_sec * (uint64_t)1000000 + tv.tv_usec;
	return (time);
}

// 
// The receive thread:
//
char *recv_data;
volatile int recv_not_ready = 1;
void *receive(void *param)
{
	ScifRecvQueue* p = (ScifRecvQueue *)param;
	unsigned int num_page = (capacity * msg_size)/PAGE_SIZE;
	int err;
	uint64_t time;
	int tot_recv = 0;
	int num_recv = 0;

	if ((err = posix_memalign((void **)&recv_data, 0x1000, num_page * PAGE_SIZE)) < 0) {
		cout << "receive thread: posix_memalign() failed " << endl;
		exit(1);
	}

	memset(recv_data, 0xa, num_page * PAGE_SIZE);
	err = mlock(recv_data, num_page * PAGE_SIZE);
	if (err) {
		cout << "recieve thread: mlock failed" << endl;
		// error handling left out for simplicity
	}

	// Let the main thread start sending packets to the card now.
	recv_not_ready = 2;
	while (tot_recv < num_msg) {
		// get one packet
		num_recv = p->pop(recv_data, msg_size, 1);
		if (num_recv == 1)  {
			recvTime[tot_recv] = get_time_stamp();
			tot_recv += num_recv;
		}
	}

	cout << "Host side recieve thread processed " << tot_recv << " packets" << endl;

	return (void *)0;
}

void usage(char *p)
{
	cout << p << " -s <message size> -c <queue capacity> -t <time interval> -n <total messages>" << endl;
	cout << "Arguments: " << endl;
	cout << "\t-s message size - size of each message (4KB max) sent between host and coprocessor (1000)" << endl;
	cout << "\t-c queue capacity - number of elements the queue should be sized to (100)" << endl;
	cout << "\t-t time interval - inter packet gap (IGP) (10 microseconds) " << endl;
	cout << "\t-n total messages - total number of messages that should be sent from the host (100)" << endl;
}

int main(int argc, char *argv[])
{
	int err, time;
	cpu_set_t cpuset;
	char *data;

	msg_size = 1000;
	capacity = 100;
	time = 10;
	num_msg = 100;

	// Optimization #6: Single thread affinity.
	CPU_ZERO(&cpuset);
	CPU_SET(6, &cpuset);

	if (sched_setaffinity(0, sizeof(cpuset), &cpuset) < 0) {
		cout << "failed to set affinity of main thread" << endl;
		// error handling left out for simplicity
	}

	if (argc != 9) {
		usage(argv[0]);
		exit(1);
	}

	// no default values on purpose!
	while ((err = getopt(argc, argv, "s:c:t:n:")) != -1) {
		switch(err) {
			case 's':
				msg_size = atoi(optarg);
				break;
			case 'c':
				capacity = atoi(optarg);
				break;
			case 't':
				time = atoi(optarg);
				break;
			case 'n':
				num_msg = atoi(optarg);
				break;
			case '?':
				usage(argv[0]);
				exit(1);
		}
	}

	cout << "Host setup as client" << endl;
	cout << "msg_size " << msg_size << " capacity " << capacity << " time " << time 
					<< " num_msg " << num_msg << endl;

	// I'm the host!
	ScifClient Host("host", SCIF_DEVICE_HOST);

	// connect to mic0 - this sample assumes we have atleast mic0
	const unsigned int server_port = 2050;
	Host.connect("mic0", server_port, 2);

	// Create a Send Queue
	ScifSendQueue SendQueue("host_send_queue", (unsigned int)msg_size, (unsigned int)capacity, Host.pop_epd());
	SendQueue.pair("mic_recv_queue");

	// Create a Recv Queue
	ScifRecvQueue RecvQueue("host_recv_queue", (unsigned int)msg_size, (unsigned int)capacity, Host.pop_epd());
	RecvQueue.pair("mic_send_queue"); 

	cout << "Number of messages: " << num_msg << endl;
	
	// Optimization #6: Single thread affinity
	// Create a single recieve thread for processing packets coming back from the coprocessor.
	CPU_ZERO(&cpuset);
	CPU_SET(8, &cpuset);
	
	pthread_t thread;
	pthread_attr_t attr;
	struct sched_param param;

	pthread_attr_init(&attr);
	param.sched_priority = 20;

	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
	pthread_attr_setschedpolicy(&attr, SCHED_FIFO);
	pthread_attr_setinheritsched(&attr, PTHREAD_EXPLICIT_SCHED);
	pthread_attr_setschedparam(&attr, &param);

	if (pthread_attr_setaffinity_np(&attr, sizeof(cpuset), &cpuset) != 0) {
		cout << "unable to setaffinity for recv thread" << endl;
		// error handling left out for simplicity
	}

	memset(recvTime, 0xa, MAX_N * sizeof(uint64_t));
	err = mlock(recvTime, MAX_N * sizeof(uint64_t));
	if (err) {
		cout << "recvTime: mlock failed" << endl;
			
	}

	memset(sendTime, 0xa, MAX_N * sizeof(uint64_t));
	err = mlock(sendTime, MAX_N * sizeof(uint64_t));
	if (err) {
		cout << "SendTime: mlock failed" << endl;
	}

	if ((err = posix_memalign((void **)&data, 0x1000, PAGE_SIZE)) < 0) {
		cout << "posix_memalign failed for data" << endl;
		exit(1);
	}

	memset(data, 0xa, PAGE_SIZE);

	// lock it so we don't page fault. This will require sufficient privileges to work.
	err = mlock(data, PAGE_SIZE);
	if (err) {
		cout << "data: mlock failed" << endl;
		// error handling left out for simplicity
	}

	pthread_create(&thread, &attr, &receive, (void*)&RecvQueue);

	// wait until recv thread has spun up!
	while (recv_not_ready == 1);

	// Number of random numbers in each msg. The first slot of
	// every message will have the timing information.
	int num_rn = msg_size/sizeof(uint64_t) - 1;
	
	uint64_t tdata;
	int tot_send = 0; 
	int num_send = 0;

#ifdef VERBOSE_DEBUG
	cout << "Starting packet transmissions to coprocessor " << endl;
#endif
	// Start queueing up packets for the coprocessor.
	while (tot_send < num_msg) {
		// store the send time.
		tdata = get_time_stamp();

		// record the send time for packet 'i'
		sendTime[tot_send] = *(uint64_t *)data = tdata;

		for (int j = 1; j < num_rn; ++j) {
			// some random data for this example.
			*((uint64_t *)data + j) = 0xDEADC0DE + tot_send;
		}		

		num_send = SendQueue.push((char *)data, msg_size);
		if (num_send == 0) 
			continue;

		if (num_send != 1) {
			cout << "Number of sent packets != 1\n" << endl;
		}

		// generate inter-packet gap.
		if (time > 0)
			microSleep(time);

		tot_send += num_send;
	}

#ifdef VERBOSE_DEBUG
	cout << "Waiting for receive thread to join " << endl;
#endif
	pthread_join(thread, 0);

	bool flag = true;

	//----------------------------------------------------------------------------
	// 				Statistics
	//----------------------------------------------------------------------------
	if (flag) {
		cout << "Total RTT from first send to last receive for " << num_msg << " messages: " 
				<< recvTime[num_msg - 1] - sendTime[0] << " microsecs" << endl;

		vector<double> rtt;
		vector<double> send_int;

		for (int i = 0; i < num_msg; ++i) {
			
			// RTT[i] 
			rtt.push_back(recvTime[i] - sendTime[i]);
			if (i != num_msg - 1) {
				send_int.push_back(sendTime[i + 1] - sendTime[i]);	
			}
		}

		// This is so we don't see a ton of data. Only the first 20 and last 20 please.
		if (num_msg >= 40) {
			cout << "RTT for the first 20 messages: ";
			for (int i = 0; i < 20; ++i) 
				cout << rtt[i] << " ";

			cout << endl;
			cout << "RTT for the last 20 messages: ";
			for (int i = 0; i < 20; ++i)
				cout << rtt[num_msg - 20 + i] << " ";
			cout << endl;
		}

		// Mean send and recv time.
		double tot = 0;
		for (size_t i = 0; i < rtt.size(); ++i) {
			tot += rtt[i];
		}

		double rtt_mean = tot / rtt.size();

		tot = 0;
		for (size_t i = 0; i < send_int.size(); ++i) {
			tot += send_int[i];
		}

		double send_int_mean = tot / send_int.size();

		// sort the time vectors
		sort(rtt.begin(), rtt.end());
		sort(send_int.begin(), send_int.end());

		double percentiles[] = {50, 70, 75, 80, 85, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 99.9, 99.99, 99.999, 99.9999};
		int psize = sizeof(percentiles)/sizeof(double);

		const int num_print = 2;
		string title[num_print] = {"Mean RTT: ", "Mean time interval between sends: "};
		double mean[num_print] = {rtt_mean, send_int_mean};
		vector<double>* pv[num_print] = {&rtt, &send_int};	

		for (int i = 0; i < num_print; ++i) {
			cout << title[i] << mean[i] << " microsecs" << endl;
			for (int j = 0; j < psize; j++) {
				printf("%9.4f", percentiles[j]);
				printf("%%");
			}
			cout << endl;

			for (int j = 0; j < psize; j++) {
				size_t pos = (size_t)(0.01 * percentiles[j] * pv[i]->size());
				if (pos >= pv[i]->size()) pos = pv[i]->size() - 1;
					printf("%10d", (int)(*pv[i])[pos]);
			}
			cout << endl;
		}
	}

	return 0;
}
