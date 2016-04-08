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

#include <scif.h>
#include "scif_tutorial.h"

#define START_OFFSET 0x80000
#define PAGE_SIZE 0x1000
#define PORT_NO 2050
#define BACKLOG 5

/* Barrier: Nothing more than receive a "control_msg" and send a "control_msg" back.
 * What is send/recv'd isn't important 
 */
#define BARRIER(newepd, string) { \
	printf("%s\n", string); \
	if ((err = scif_send(newepd, &control_msg, sizeof(control_msg), 1)) <= 0) { \
		printf("scif_send failed with err %d\n", get_curr_status()); \
		fflush(stdout); \
		goto __end; \
	} \
	if ((err = scif_recv(newepd, &control_msg, sizeof(control_msg), 1)) <= 0) { \
		printf("scif_recv failed with err %d\n", get_curr_status()); \
		fflush(stdout); \
		goto __end; \
	} \
}

typedef struct window_info {
	void *self_addr;
	off_t offset;
} win_t;

void usage(char *p)
{
	printf("%s -n <4K page count>\n", p);
	printf("Argunments\n");
	printf("\tcount of 4K page: Number of 4K pages to scif_mmap() \n");
}


int main(int argc, char *argv[])
{
	scif_epd_t epd, newepd;
	struct scif_portID portID;
	off_t suggested_offset;
	win_t buffer;
	char ack_msg[32]={0};
	char *temp, *end_addr;
   	volatile char *curr_addr;
	int j, err, conn_port, num_buf, self=0;
	int control_msg, msg_size, map_fixed, mark;

	if (argc != 3) {
		usage(argv[0]);
		exit(1);
	}

	num_buf = atoi(argv[2]);

	msg_size = (num_buf + 1) * PAGE_SIZE;
	if ((msg_size <= 0) || (msg_size > INT_MAX)) {
		printf("not valid msg size");
		exit(1);
	}

	/* open end pt */
	if ((epd = scif_open()) == SCIF_OPEN_FAILED) {
		printf("scif_open failed with err %d\n", get_curr_status());
		exit(1);
	}

	/* bind end pt to specified port */
	if ((conn_port = scif_bind(epd, PORT_NO)) < 0) {
		printf("scif_bind failed with error %d\n", get_curr_status());
		exit(2);
	}

	/* marks an end pt as listening end pt and queues up a maximum of BACKLOG
	 * no: of incoming connection requests
	 */
        printf("Listening for incoming connection(s) ... ");
	if (scif_listen(epd, BACKLOG) != 0) {
		printf("scif_listen failed with error %d\n", get_curr_status());
		exit(1);
	}

        printf("Accepting connection ... ");
	/* accepts a conn request by creating a new end pt that connects to peer */
	if (scif_accept(epd, &portID, &newepd, SCIF_ACCEPT_SYNC) != 0) {
		printf("scif_accept failed with error %d\n", get_curr_status());
		exit(1);
	}
        printf("done \n");

	 /* addresses in VAS & RAS must be multiple of page size */
	err = posix_memalign(&buffer.self_addr, 0x1000, msg_size);
	if (err != 0) {
		printf("posix_memalign failed with error : %d\n", err);
		goto __end;
	}
	memset(buffer.self_addr, 0x0, msg_size);

	err = posix_memalign((void **)&temp, 0x1000, PAGE_SIZE);
	if (err != 0) {
		printf("posix_memalign failed with error : %d\n", err);
		goto __end;
	}
	memset(temp, 0x0, PAGE_SIZE);

	printf("Registering memory of size %d bytes ... ", msg_size);
	/* scif_register : marks a memory region for remote access */
	if ((buffer.offset = scif_register (newepd,
					    buffer.self_addr,
					    msg_size,
					    START_OFFSET,
					    SCIF_PROT_READ | SCIF_PROT_WRITE,
					    SCIF_MAP_FIXED)) < 0) {
		printf("scif_register failed with error : %d\n", get_curr_status());
		goto __end;
	}

	printf("Done (offset 0x%lx)\n", (unsigned long)buffer.offset);

    	/* See note in the other code, but this is where this code is telling the other
 	 * side of the connection that it is done with registration and is ready for action 
 	 */
	BARRIER(newepd, "Card: I'm done with registration!");

	/* Now wait for the other side to write some data and let us know when it is done */
	scif_recv(newepd, ack_msg, sizeof(ack_msg), 1);

	/* Reversing the written data in PAGE_SIZE chunks */
	curr_addr = (volatile char *)buffer.self_addr;
	/* Last chunk of RAS window is reserved for value update by scif_fence_signal */
	end_addr = (char *)buffer.self_addr + msg_size - (2 * PAGE_SIZE);

        /* Quick check: Now that the other side has told us that it is done with moving bytes, 
 	 * does the first byte of the first page have 0x1 and does the first byte of the last page 
 	 * have num_buf? That's our expected data right?!
         */
	if (*curr_addr != 0x1 && *end_addr != num_buf){
		printf("Hmmm... scif_writeto didn't do its job! \n");
	} else {
		printf("Data received correctly! \n");
	}

	/* This is very stupid on purpose -- swapping page - really? */
	for (j = 0; j < (int)num_buf/2; j++) {
		memcpy(temp, curr_addr, PAGE_SIZE);
		memcpy(curr_addr, end_addr, PAGE_SIZE);
		memcpy(end_addr, temp, PAGE_SIZE);

		curr_addr = curr_addr + PAGE_SIZE;
		end_addr = end_addr - PAGE_SIZE;
	}

	/* tell the other side to suck the data back upto the host now */
	BARRIER(newepd, "Card: Data reverse done");

	curr_addr = (volatile char *)buffer.self_addr + num_buf * PAGE_SIZE; 

        /* --------------------------------------------------------------------
	 * Okay now we'll test the scif_mmap() capabilities -- the "low latency" 
         * way of doing things
         * --------------------------------------------------------------------
         */
	while (*(volatile unsigned int *)curr_addr != 0xC0DE);   /* wait for this to become "go" */
	printf("Flag transitioned to \"go\" \n");
        printf("Host wrote: %s \n", (char *)buffer.self_addr);

        /* card will now write something to page 1 */
        strcpy((char *)buffer.self_addr + PAGE_SIZE, "I'm okay, How are you doing, Host?");
	
	/* set the flag to "go" for the host */
	*(unsigned int *)curr_addr = 0xDEED;   /* wait for this to become "go" */
	
	/* Receive ACK msg from peer when its done moving data. This is so we don't unregister
         * the local buffers too soon 
         */
	printf("Waiting for final message from host \n");
	scif_recv(newepd, ack_msg, sizeof(ack_msg), 1);
	printf("Received %s \n", ack_msg);
	

	/* scif_unregister : closes the registered window */
	if ((err = scif_unregister(newepd, buffer.offset, msg_size)) < 0) {
		printf("scif_unregister failed with error : %d\n", get_curr_status());
		goto __end;
	}

	BARRIER(newepd, "Card: Window unregistered");
	errno = 0;

__end:
	if (buffer.self_addr != NULL) 
		free(buffer.self_addr);

	if (temp != NULL)
		free(temp);

	scif_close(newepd);
	scif_close(epd);

	if (errno == 0)
		printf("Success\n");
	else
		printf("Failed\n");

	return errno;
}

