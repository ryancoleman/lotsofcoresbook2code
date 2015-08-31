#!/bin/sh

# Example for low latency optimizations.
#
#DEBUG_FLAGS=-DVERBOSE_DEBUG

all: scif_accept_rbuf_mic scif_connect_rbuf_host

# card side compiled with icc.
scif_accept_rbuf_mic: ScifNode.cpp ScifNode.h ScifQueue.cpp ScifQueue.h scif_accept_rbuf.cpp 
	icc -mmic -O3 -I /usr/include  -Wall ScifQueue.cpp ScifNode.cpp  scif_accept_rbuf.cpp -lscif -lpthread -o scif_accept_rbuf_mic $(DEBUG_FLAGS)

# host side compiled with icpc
scif_connect_rbuf_host: ScifNode.cpp ScifNode.h ScifQueue.cpp ScifQueue.h scif_connect_rbuf.cpp
	icpc -g -O3 -I /usr/include  -Wall -DHOST ScifQueue.cpp ScifNode.cpp scif_connect_rbuf.cpp -lscif -lpthread -o scif_connect_rbuf_host $(DEBUG_FLAGS)

clean:
	rm -f scif_connect_rbuf_host scif_accept_rbuf_mic
