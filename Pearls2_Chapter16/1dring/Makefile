CC=mpiicc

BIN=mpi3shm_1Dring.v2.1

all: host mic 

host: mpi3shm_1Dring.v2.1.c
	$(CC) -o $(BIN)   $(FLAGS) $<


mic: mpi3shm_1Dring.v2.1.c
	$(CC) -mmic -o $(BIN).mic  $(FLAGS) $<


run-host:
		@mpirun -n 2 ./$(BIN)

run-mic:
		@I_MPI_MIC=1 mpirun -n 2 -host $(shell hostname -s)-mic0 -genv I_MPI_MIC_POSTFIX .mic ./$(BIN)
		
run-symm:
		@I_MPI_MIC=1 mpirun -n 4 -ppn 2 -hosts $(shell hostname -s),$(shell hostname -s)-mic0 -genv I_MPI_MIC_POSTFIX .mic ./$(BIN)

clean:
		rm -f $(BIN) $(BIN).mic

