#
# Simple example to demonstrate a SCIF connection.
#

CARD_CC=icc -mmic
HOST_CC=icc

# Use gcc instead?
#CARD_CC=/usr/linux-k1om-4.7/bin/x86_64-k1om-linux-gcc
#HOST_CC=gcc

# everything
all: scif_accept_register scif_connect_register

# Host side
scif_accept_register: scif_accept_register.c scif_tutorial.h scif.h
	$(HOST_CC) scif_accept_register.c -lscif -L/usr/lib64 -o scif_accept_register

# MIC side
scif_connect_register: scif_connect_register.c scif_tutorial.h scif.h
	$(CARD_CC) scif_connect_register.c -lscif -o scif_connect_register 

clean:
	rm -f scif_connect_register scif_accept_register

