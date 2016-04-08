CC=icc
CFLAGS=-mmic -O3 -openmp -DSTREAM_ARRAY_SIZE=64000000 -opt-streaming-cache-evict=0 -opt-streaming-stores always
LIBHUGETLB_DIR=libhugetlbfs-2.16
TARGET=stream
DEPS=stream.c

INTEL_LIBS=libiomp5.so
SHARED_LIBS=

ENV=KMP_AFFINITY=scatter OMP_NUM_THREADS=60 LD_LIBRARY_PATH=/tmp


ifeq ($(HUGE_PAGES), yes)
CFLAGS+=-Wl,-T $(LIBHUGETLB_DIR)/ldscripts/elf_x86_64.xBDT -L $(LIBHUGETLB_DIR)/obj64
TARGET:=$(TARGET).huge
DEPS+=$(LIBHUGETLB_DIR)/obj64/libhugetlbfs.so
SHARED_LIBS+=$(LIBHUGETLB_DIR)/obj64/libhugetlbfs.so
INTEL_LIBS+=libsvml.so libintlc.so.5 libintlc.so libintlc.so libimf.so libirng.so
endif

ifeq ($(PREFETCH), no)
CFLAGS+=-no-opt-prefetch
TARGET:=$(TARGET).nopf
else
CFLAGS+=-opt-prefetch-distance=64,8
TARGET:=$(TARGET).pf
endif


$(TARGET):$(DEPS)
	$(CC) $(CFLAGS) -o $@ $<

$(LIBHUGETLB_DIR)/obj64/libhugetlbfs.so:
	make -C $(LIBHUGETLB_DIR) ARCH=x86_64 CC64='icc -mmic' libs BUILDTYPE=NATIVEONLY

run:$(TARGET)
ifneq ($(RUNONLY), yes)
	scp $(addprefix $(MKLROOT)/../compiler/lib/mic/, $(INTEL_LIBS)) $(SHARED_LIBS) $< mic0:/tmp
ifeq ($(HUGE_PAGES), yes)	
	ssh mic0 'echo 623 > /proc/sys/vm/nr_hugepages'
	ssh mic0 'mkdir -p /mnt/hugetlbfs'
	ssh mic0 'mount -t hugetlbfs none /mnt/hugetlbfs'
endif
endif
	ssh mic0 $(ENV) /tmp/$<

clean:
	rm -f stream.*pf
	make -C $(LIBHUGETLB_DIR) clean
