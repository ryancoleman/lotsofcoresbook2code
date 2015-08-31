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
#include "scif.h"
#include "scif_tutorial.h"

#define START_OFFSET 0x80000
#define PAGE_SIZE 0x1000
#define PEER_PORT 2050

/* Note that all this barrier is doing is sending a "control_msg" to the other side of this connection
 * and waiting (scif_recv) for a "control_msg" from the remote side
 */
#define BARRIER(epd, string) { \
	printf("%s\n", string); \
	if ((err = scif_send(epd, &control_msg, sizeof(control_msg), 1)) <= 0) { \
		printf("scif_send failed with err %d\n", get_curr_status()); \
		fflush(stdout); \
		goto __end; \
	} \
	if ((err = scif_recv(epd, &control_msg, sizeof(control_msg), 1)) <= 0) { \
		printf("scif_recv failed with err %d\n", get_curr_status()); \
		fflush(stdout); \
		goto __end; \
	} \
}

typedef struct window_info {
	void *self_addr;
	off_t offset;
}win_t;

void usage(char *p)
{
	printf("%s -n <4K page count> -c <cpu vs. dma 1/0> -r <remote node>\n", p);
	printf("Argunments\n");
	printf("\tcount of 4K page: Number of 4K pages to scif_mmap() \n");
	printf("\tcpu vs. DMA: Method of data transfer \n");
	printf("\tremote node: Remote node to connect to \n");
}

int main(int argc, char *argv[])
{
	scif_epd_t epd;
	struct scif_portID portID;
	off_t curr_offset, signal_offset;
	win_t buffer;

	volatile char *curr_addr;
	char ack_msg[32];
	int j, err, data, num_buf, control_msg, tries = 20;
	int conn_port, msg_size, map_fixed, use_cpu;
	char *remotep;

	if (argc != 7) {
		usage(argv[0]);
		exit(1);
	}

	num_buf = atoi(argv[2]); /* how many 4K pages? */
	use_cpu = atoi(argv[4]);  /* CPU or DMA for transfer? */
	portID.node = atoi(argv[6]);  /* remote node ID */

	msg_size = (num_buf + 1) * PAGE_SIZE; /* total buffer size, one extra for sync */
	if (msg_size <= 0 || msg_size > INT_MAX) {
		printf("not valid msg size");
		exit(1);
	}

	printf("Params: msg_size %d bytes, using %s remote node %d \n", 
				num_buf * PAGE_SIZE, (use_cpu == 0) ? "DMA" : "CPU", portID.node);

	portID.port = PEER_PORT;

  	/* The code here is similar to BSD sockets. You want to connect to the remote side before
  	 * doing anything else. 
  	 * 1. scif_open
  	 * 2. scif_bind to a port
  	 * 3. scif_connect. Will block until the connection is accepted on the other side
  	 *
  	 * I've added printf's just so you know what is going on .. you can take them out.
  	 */
	if ((epd = scif_open()) == SCIF_OPEN_FAILED) {
		perror("scif_open failed");
		exit(1);
	}

	/* bind end point to available port, generated dynamically */
	if ((conn_port = scif_bind(epd, 0)) < 0) {
		perror("scif_bind failed");
		exit(1);
	}

	/* Initiate a connection to remote node, when successful returns the
	 * peer portID. Try to connect about 20 times, and give up. The other side isn't 
	 * accepting -- did you start the other side?
	 */
        printf("Estabilishing connection ... ");
__retry:
	if (scif_connect(epd, &portID) < 0) {
		if ((get_curr_status() == ECONNREFUSED) && (tries > 0)) {
			printf("connection to node %d failed : trial %d\n", portID.node, tries);
			tries--;
			sleep(1);
			goto __retry;
		}
		perror("scif_connect failed");
		exit(1);
	}

	printf("done \n");

	/* Allocate a buffer. Addresses in Virtual Address Space (VAS) and 
         * Registered Address Space (RAS) must be multiple of PAGE_SIZE (0x1000)
         */
	err = posix_memalign(&buffer.self_addr, 0x1000, msg_size);
	if (err != 0) {
		printf("posix_memalign failed with error : %d\n", err);
		goto __end;
	}

	/* Register the allocated memory. This will return a "offset". This offset is useful
         * for all kinds of things as we will see.
         *
         * Note1: Note the SCIF_PROT flags (read/write) and whether or not we want to ask for a 
         * fixed offset (MAP_FIXED) so we don't have to tell the other side about it.
         */
        printf("Registering memory of size %d bytes ...", msg_size);
	if ((buffer.offset = scif_register (epd,
					    buffer.self_addr,
					    msg_size,
					    START_OFFSET,
					    SCIF_PROT_READ | SCIF_PROT_WRITE,
					    SCIF_MAP_FIXED)) < 0) {
		perror("scif_register failed");
		goto __end;
	}

	printf("Done (offset 0x%lx)\n", (unsigned long)buffer.offset);

        /* Now here is something you have to pay attention to.
  	 * We need to make sure the other side has allocated its memory and registerd it too.
  	 * We use the BARRIER mechanism to wait for an "ack" from the other side. There is nothing
  	 * more to it!
  	 */
	BARRIER(epd, "Host: I'm done with registration!");

        /* What the loop below does is the following:
 	 * 1. memset a page worth of buffer to 0x1, 0x2, 0x3 .. 
 	 * 2. scif_writeto the remote registered address space (i.e. remote buffer). All this does
 	 *    is use either the CPU or the DMA engine to program the transfer. If the DMA engine
 	 *    is used, you have to make sure the transfer is complete by fencing (See below).
 	 */
	data = 0x1;
	curr_offset = buffer.offset;
	curr_addr = (volatile char *)buffer.self_addr;

        printf("Moving bytes now ... \n");

	for (j = 0; j < num_buf; j++) {
		memset(curr_addr, data, PAGE_SIZE);

		if ((err = scif_writeto(epd,
					curr_offset, /* local RAS offset */
					PAGE_SIZE,
					curr_offset, /* remote RAS offset */
					use_cpu ? SCIF_RMA_USECPU : 0))) {
			perror("scif_writeto failed");
			goto __end;
		}

		/* prints the starting value of PAGE_SIZE chunk to verify data */
		printf("%p[0] : %lx\n", curr_addr, (unsigned long)((int *)curr_addr)[0]);

		/* next page .. with the next pattern (see step 1 above) */
		data = data + 1;
		curr_addr = curr_addr + PAGE_SIZE;
		curr_offset = curr_offset + PAGE_SIZE;
	}

        /* NOTE: The other side is waiting (via scif_recv) for data transfer to be complete now, 
 	 * until we let it go (via send or write a "status flag". for simplicity now, we'll send
	 * a message. I'll show how to use scif_mmap() to write a flag shortly
	 */ 

	/* Remember we had allocated one extra page, we'll use it for "signaling" DMA completions.
         * We need to tell scif_fence_signal what offset to write the final message when DMA is done.
         * This has to be an offset. We use the last page's offset for this
         */
	signal_offset = buffer.offset + (num_buf * PAGE_SIZE);
	curr_addr = (volatile char *)buffer.self_addr + (num_buf * PAGE_SIZE);

        printf("Fencing DMA Write to card ... ");
	if ((err = scif_fence_signal(epd,
				signal_offset,	/* local offset */
				0xabcd,		/* value to be written to local offset */
				signal_offset,	/* remote offset */
				0xabcd,		/* value to be written to remote offset */
				SCIF_FENCE_INIT_SELF | SCIF_SIGNAL_LOCAL)) != 0) {
		perror("scif_fence_signal failed"); 
		goto __end;
	}
	printf("Done \n");

	/* Wait for DMA to be done, really.
 	 * Note that I've used the offset to tell HW what to write at the page (the offset describes 
 	 * it). But I'll use the virtual address of the "status" page to spin for completion.
 	 * Also note, this is a spin in user mode - I'm not spinning or blocking in the kernel. 
 	 * Talk about latency!
 	 */

	/* Now this is the thing we talked about in the call Elmer/Rama.
 	 * Notice how I called scif_writeto() above, then scif_fence_signal and then waited on the 
 	 * result. If I had something on the host (computation) that I could do in the meantime,
 	 * while the DMA engine was moving bytes, I could have done it in between the 
 	 * scif_fence_signal and the spin in user mode.
 	 *
 	 * If on the other hand I don't have anything to do, I can simplify this. I can pass in the
 	 * SCIF_RMA_SYNC flag into scif_writeto (see scif.h for more) and in that case scif_writeto
 	 * will block until the DMA is completely done.
 	 *
 	 * Hope this makes sense.
 	 */
	printf("Spinning for data transfer completion ...");
	while (*(unsigned short *)curr_addr != 0xabcd);
	printf("Done\n");
	
        
	/* Send a msg to peer about RMA completion -- the data is on the card now */
	strcpy(ack_msg, "DMA done");
	scif_send(epd, ack_msg, sizeof(ack_msg), 1);

	BARRIER(epd, "Host: Waiting for card to reverse the data");

        /* Now we use the readfrom API to read data back from the other side */
	curr_offset = buffer.offset;
        printf("Reading data from remote side ... ");
	for (j = 0; j < num_buf; j++) {
		if ((err = scif_readfrom(epd,
					 curr_offset, /* local RAS offset */
					 PAGE_SIZE,
					 curr_offset, /* remote RAS offset */
					 use_cpu ? SCIF_RMA_USECPU : 0))) {
			perror("scif_readfrom failed");
			goto __end;
		}
		curr_offset = curr_offset + PAGE_SIZE;
	}
        printf("Done \n");

        printf("Fencing DMA Read from card ... ");
	/* Fencing against scif_readfrom operations initiated from local endpt*/
	if ((err = scif_fence_signal(epd,
				     signal_offset,	/* local offset */
				     0xabcd,		/* value to be written to local offset */
				     signal_offset,	/* remote offset */
			             0xabcd,		/* value to be written to remote offset */
			             SCIF_FENCE_INIT_SELF | SCIF_SIGNAL_LOCAL)) != 0) {
		perror("scif_fence_signal failed");
		goto __end;
	}
	printf("Done \n");

	printf("Spinning for data transfer completion ...");
	while (*(volatile unsigned short *)curr_addr != 0xabcd);
	printf("Done\n");

	/* prints the starting value of PAGE_SIZE chunk to verify data */
	curr_addr = (volatile char *)buffer.self_addr;
	for (j = 0; j < num_buf; j++) {
		printf("curr_addr[0] : %lx\n", (unsigned long)((int *)curr_addr)[0]);
		curr_addr = (char *)curr_addr + PAGE_SIZE;
	}

	/* ---------------------------------------------------------------------------------------
 	 * Now onto using scif_mmap(). Since we have n pages registered + 1 page for status here's
 	 * how we'll just use them. First the host will write a string "Hello Card, How are you?"
 	 * to the first page *directly* on the card. To do this, I'll need a pointer to a buffer
 	 * on the card. Then I'll trigger the status page by writing a flag (0xC0DE) to the "status"
 	 * page (which is the 1 extra page I had). The card, meanwhile, will spin on the flag to 
 	 * become 0xC0DE (from whatever it was) and when it does, will print the string 
 	 * Then the card will write a string and flip a flag (0xDEED) and I'll wait on it. 
 	 *
 	 * And so on...
 	 * ---------------------------------------------------------------------------------------
 	 */
        
        /* Step1: get a pointer to the card side (registered) memory. I'll reuse the offset I got 
         * earlier
         */
	printf("Getting pointer to buffer on card ... ");
       	remotep = scif_mmap(NULL, msg_size, SCIF_PROT_READ|SCIF_PROT_WRITE, 0, epd, START_OFFSET);
	if (remotep == SCIF_MMAP_FAILED) {
		perror("scif_mmap failed");
		goto __end;
	}
	printf("Done\n");

	
	/* Step2: write something to page 0 */
	strcpy(remotep, "Hello Card, How are you?");

	/* Step3: flip a flag on the card */
	curr_addr = remotep + num_buf * PAGE_SIZE;

	__asm __volatile("sfence"::);

	/* "go" */
	*(volatile unsigned int *)curr_addr  = 0xC0DE;  
	
	/* Since card memory is mapped WC, we need to make
         * sure the above write is flushed from the WC buffer in the CPU
         * Any serializing instruction will do, I'm using sfence.
         */	
	__asm __volatile("sfence"::);

	printf("Waiting for card to 0xDEED \n");
	/* Step4: Wait for the card to signal you back */
	while(*(volatile unsigned int *)curr_addr != 0xDEED);

	/* remember the card to wrote to page 1 -- just for fun */
        printf("Received from Card: %s \n", remotep + PAGE_SIZE); 

	/* Good bye, card! */
	strcpy(ack_msg, "Good Bye");
	scif_send(epd, ack_msg, sizeof(ack_msg), 1);

	/* scif_unregister : closes the registered window */
	if ((err = scif_unregister(epd, buffer.offset, msg_size)) < 0) {
		printf("scif_unregister failed with error : %d\n", get_curr_status());
		goto __end;
	}
	BARRIER(epd, "Host: Window unregistered");
	errno = 0;

__end:
	if (buffer.self_addr != NULL)
		free(buffer.self_addr);

	scif_close(epd);

	if (errno == 0)
		printf("Success\n");
	else
		printf("Failed\n");

	return errno;
}
