# Intel
export CPATH=$CPATH:/usr/include

export CC  := mpiicc 
export CFLAGS := -mkl=sequential -g -O3 -xHost -mcmodel=medium -openmp -align -ansi-alias -restrict -fno-alias -std=gnu99
export CFLAGS += -DL1_CUTOFF=3 -I../resources/iniparser/src -I/usr/include
export LIBS := -L../resources/iniparser -lgsl -liniparser
