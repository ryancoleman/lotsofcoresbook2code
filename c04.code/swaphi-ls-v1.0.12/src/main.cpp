/*
 * main.cpp
 *
 *  Created on: Jun 19, 2013
 *      Author: yongchao
 */
#include "Align.h"

int main(int argc, char* argv[]) {
	Align* align;

#ifdef MPI_PARALLEL
	MPI_Init(&argc, &argv); /*initialize the MPI runtime environment*/
#endif

	/*create alignment object*/
	align = new Align();

	/*parse the parameter*/
	if (align->parseParams(argc, argv) == false) {

		/*release the object*/
		delete align;

#ifdef MPI_PARALLEL
		MPI_Finalize();
#endif
		return -1;
	}

	/*run the kernel*/
	align->run();

	/*release the object*/
	delete align;

#ifdef MPI_PARALLEL
	MPI_Finalize(); /*finalize the MPI*/
#endif

	return 0;
}

