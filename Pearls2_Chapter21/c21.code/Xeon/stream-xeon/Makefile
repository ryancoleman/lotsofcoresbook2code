#CC = gcc
#CFLAGS = -O2
#
#FF = g77
#FFLAGS = -O2

CC=icc

CFLAGS=-mcmodel medium -shared-intel -O3 -xhost -DSTREAM_ARRAY_SIZE=120000000 -DOFFSET=0 -openmp -opt-streaming-stores always
#CFLAGS=-mcmodel medium -shared-intel -O3 -xAVX -DSTREAM_ARRAY_SIZE=120000000 -DOFFSET=0 -openmp -opt-streaming-stores always

FF=ifort
FFLAGS=-O2

.PHONY: all stream clean

all: stream

stream: stream_c.exe

stream_f.exe: stream.f mysecond.o
	$(CC) $(CFLAGS) -c mysecond.c
	$(FF) $(FFLAGS) -c stream.f
	$(FF) $(FFLAGS) stream.o mysecond.o -o stream_f.exe

stream_c.exe: stream.c
	$(CC) $(CFLAGS) stream.c -o stream_c.exe

clean:
	rm -f stream_f.exe stream_c.exe *.o
