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

/* 
 * scif_accept_register.c : demonstrates basic implementation fo registering windows
 * and mmap calls
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <unistd.h>
#include <fcntl.h>
#include <string.h>
#include <limits.h>
#include <sys/mman.h>
#include <sys/ioctl.h>
#include <scif.h>

#define START_OFFSET 0x80000
#define PAGE_SIZE 0x1000
#define PORT_NO 2050
#define BACKLOG 5
#define LOCAL_PATTERN 0xbc
#define REMOTE_PATTERN 0xadadadad

#define BARRIER(newepd, string) { \
	printf("%s\n", string); \
	if ((err = scif_send(newepd, &control_msg, sizeof(control_msg), 1)) <= 0) { \
		printf("scif_send failed with err %d\n", errno); \
		fflush(stdout); \
		goto __end; \
	} \
	if ((err = scif_recv(newepd, &control_msg, sizeof(control_msg), 1)) <= 0) { \
		printf("scif_recv failed with err %d\n", errno); \
		fflush(stdout); \
		goto __end; \
	} \
	printf("==============================================================\n"); \
}

typedef struct window_info {
	void *peer_addr;
	void *self_addr;
	off_t offset;
}win_t;

void usage(char *p)
{
	printf("%s -n <4K page count> -m <mmap chooser 0/1/2>\n", p);
	printf("Argunments\n");
	printf("\tcount of 4K page: Number of 4K pages to scif_mmap() \n");
	printf("\tmmap chooser: method to obtain the offset via scif_mmap \n");
}

int main(int argc, char *argv[])
{
	scif_epd_t epd, newepd;
	struct scif_portID portID;
	off_t suggested_offset = 0, peer_offset;
	win_t buffer;

	int j, err, conn_port;
	int control_msg, msg_size, map_manager;

	if (argc != 5) {
		usage(argv[0]);
		exit(1);
	}

	/* number of 4K pages * PAGE_SIZE */
	msg_size = atoi(argv[2]) * PAGE_SIZE;
	if (msg_size <= 0 || msg_size > INT_MAX) {
		printf("invalid message size");
		exit(1);
	}

	/* method of specifying the SCIF mmap offset */
	map_manager = atoi(argv[4]);

	/* open end pt */
	if ((epd = scif_open()) == SCIF_OPEN_FAILED) {
		perror("scif_open failed");
		exit(1);
	}

	/* bind end pt to specified port */
	if ((conn_port = scif_bind(epd, PORT_NO)) < 0) {
		perror("scif_bind failed");
		exit(1);
	}

	printf("bind success to port %d\n", conn_port);

	/* Marks an end pt as listening end pt and queues up a maximum of BACKLOG
	 * no: of incoming connection requests
	 */
	if (scif_listen(epd, BACKLOG) != 0) {
		perror("scif_listen failed");
		exit(1);
	}

	/* accepts a conn request by creating a new end pt that connects to peer */
	if (scif_accept(epd, &portID, &newepd, SCIF_ACCEPT_SYNC) != 0) {
		perror("scif_accept failed");
		exit(1);
	}
	printf("connection accepted!\n");

	/* addresses in VAS & RAS must be multiple of page size */
	err = posix_memalign(&buffer.self_addr, 0x1000, msg_size);
	if (err != 0) {
		perror("posix_memalign failed");
		goto __end;
	}

	memset(buffer.self_addr, LOCAL_PATTERN, msg_size);
	printf("Buffer allocated at 0x%lx\n", (unsigned long)buffer.self_addr);

	/* 
 	 * map_manager=0 : SCIF_MAP_FIXED not set, scif manages RAS; implementaion
	 * defined suitable offsets are chosen for mapping len bytes
	 *
	 * map_manager=1 : SCIF_MAP_FIXED set, user managed; specified fixed offset are
	 * used. For scif_register doesnt replace existing registrations & returns error.
	 * However, scif_mmap replaces existing mappings
	 *
	 * map_manager=2 : SCIF_MAP_FIXED set, OS managed; offset is same as virtual addr.
	 * This relieves the app from the need to track association between VA and RA as
	 * they are same.
	 */
	if (map_manager == 0)
		suggested_offset = 0;
	else if (map_manager == 1)
		suggested_offset = START_OFFSET;
	else if (map_manager == 2)
		suggested_offset = (off_t) buffer.self_addr;

	/* scif_register : marks a memory region for remote access starting at offset po,
	 * a function of suggested_offset & msg_size which backs the VAS starting at
	 * buffer.self_addr. Successful registration returns po, offset where mapping
	 * is placed
	 */
	if ((buffer.offset = scif_register (newepd,
					buffer.self_addr,
					msg_size,
					suggested_offset,
					SCIF_PROT_READ | SCIF_PROT_WRITE,
					map_manager? SCIF_MAP_FIXED : 0)) < 0) {
		perror("scif_register failed with error");
		goto __end;
	}

	printf("Registered buffer at offset 0x%lx\n", (unsigned long)buffer.offset);

	BARRIER(newepd, "Memory registration done \n");

	/* if we didn't choose to specify a fixed offset, send the offset to the other side 
 	 * for scif_mmap() 
 	 */
	if ((map_manager == 0) || (map_manager == 2)) {
		/* receive peer's window offset to mmap */
		if ((err = scif_recv(newepd, &peer_offset, sizeof(peer_offset), 1)) < 0) {
			perror("scif_recv failed");
			goto __end;
		}

		/* send registered window offset to peer to allow mmaping */
		if ((err = scif_send(newepd, &buffer.offset, sizeof(buffer.offset), 1)) < 0) {
			perror("scif_send failed");
			goto __end;
		}
	}
	else
		peer_offset = buffer.offset;

	/* scif_mmap : maps pages in VAS starting at peer_addr to remote window
	 * starting at buffer.offset where peer_addr is a function of buffer.self_addr
	 * & msg_size. successful mapping returns peer_addr, the address where
	 * mapping is placed.
	 */
	printf("Obtained peer offset: 0x%lx\n", (unsigned long)peer_offset);
	if ((buffer.peer_addr = scif_mmap(buffer.self_addr,
					msg_size,
					SCIF_PROT_READ | SCIF_PROT_WRITE,
					map_manager? SCIF_MAP_FIXED : 0,
					newepd,
					peer_offset)) == MAP_FAILED) {
		perror("scif_mmap failed with error");
		goto __end;
	}

	printf("SCIF mmap'd buffers at address 0x%lx\n", (unsigned long)buffer.peer_addr);

	/* verify that the page is correct and readable via the locally obtained mapping. 
 	 * Of course, this assumes that the remote side allocated memory and wrote the pattern
 	 * below before letting you map and read it!
 	 */
	for (j = 0; j < (msg_size/sizeof(int)); j++) {
		if (((int*)buffer.peer_addr)[j] != REMOTE_PATTERN) {
			printf("data mismatch: peer_addr[%d] %x\n",
					j, ((int*)buffer.peer_addr)[j]);
			errno = -1;
			goto __end;
		}
	}

	BARRIER(newepd, "pattern verified");

	/* scif_unmap :  removes mapping to a remote window */
	if ((err = scif_munmap(buffer.peer_addr, msg_size)) < 0) {
		perror("scif_munmap failed with error");
		goto __end;
	}

	/* scif_unregister : closes the registered window. window is not
	 * deleted until all references to it are removed
	 */
	if ((err = scif_unregister(newepd, buffer.offset, msg_size)) < 0) {
		perror("scif_unregister failed with error");
		goto __end;
	}

	BARRIER(newepd, "unmap and unregister");
	errno = 0;

__end:
	/* cleanup */
	if (buffer.self_addr != NULL)
		free(buffer.self_addr);

	scif_close(newepd);
	scif_close(epd);

	if (errno == 0)
		printf("Success\n");
	else
		printf("Failed \n");

	return errno;
}
