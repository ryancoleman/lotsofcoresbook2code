# to test w/ icc
CARD_CC=icc -mmic
HOST_CC=icc

# to test w/ gcc
#CARD_CC=/usr/linux-k1om-4.7/bin/x86_64-k1om-linux-gcc
#HOST_CC=gcc

# everything
all: scif_accept_rma_fence.o scif_connect_rma_fence.o

# Host side
scif_connect_rma_fence.o: scif_connect_rma_fence.c scif_tutorial.h scif.h
	$(HOST_CC) scif_connect_rma_fence.c -lscif -L/usr/lib64 -o scif_connect_rma_host

# MIC side
scif_accept_rma_fence.o: scif_accept_rma_fence.c scif_tutorial.h scif.h
	$(CARD_CC) scif_accept_rma_fence.c -lscif -o scif_accept_rma_mic


clean:
	rm -f scif_connect_rma_host scif_accept_rma_mic

