#if defined(__APPLE__)
#if _GLIBCXX_ATOMIC_BUILTINS == 1
#undef _GLIBCXX_ATOMIC_BUILTINS
#endif // _GLIBCXX_ATOMIC_BUILTINS
#endif // __APPLE__

#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <stdio.h>
#include <sys/mman.h>
#include <errno.h>
#include "omp.h"
#include "math.h"
#include "offload.h"
#include "Timer.h"
#include "MICStencil.cpp"

#define LINESIZE	64
#define ALLOC   	alloc_if(1)
#define FREE    	free_if(1)
#define RETAIN  	free_if(0)
#define REUSE   	alloc_if(0)

////////////////////////////////////////////////////////////////
// TODO: Tune Threads, Partitions according to card's parameters
////////////////////////////////////////////////////////////////

template <class T> __declspec(target(mic)) void
MICStencil<T>::stencilOP(  T* pIn, 
							unsigned int nIters,
							unsigned int uDimWithHalo, 
							unsigned int uImgElements,
							unsigned int uHaloWidth)
{
	T wcenter = 0.25;
	T wdiag = 0.05;
	T wcardinal = 0.15; 
	
	unsigned int uRowPartitions = sysconf(_SC_NPROCESSORS_ONLN) / 4 - 1;
	unsigned int uColPartitions = 4;	// Threads per core for KNC

	unsigned int uRowTileSize	= (uDimWithHalo - 2 * uHaloWidth) / uRowPartitions;
	unsigned int uColTileSize	= (uDimWithHalo - 2 * uHaloWidth) / uColPartitions;

	uRowTileSize = ((uDimWithHalo - 2 * uHaloWidth) % uRowPartitions > 0) ? (uRowTileSize + 1) : (uRowTileSize);

	// printf("Dividing rows int %d partitions and columns into %d partitions.\n", uRowPartitions, uColPartitions);
	// printf("Tile size rows %d, Tile size columns %d.\n", uRowTileSize, uColTileSize);
	// printf("Halo width %d, dimension with halo %d.\n", uHaloWidth, uDimWithHalo);
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	// Should use the "Halo Val" when filling the memory space

	T *pTmp	 = (T*)pIn;
	T *pCrnt = (T*)memset((T*)_mm_malloc(uImgElements * sizeof(T), LINESIZE), 0, uImgElements * sizeof(T));

//printf("%d %d\n",mlock(pTmp,uImgElements * sizeof(T)),errno);
//printf("%d %d\n",mlock(pCrnt,uImgElements * sizeof(T)),errno);
	
	#pragma omp parallel firstprivate(pTmp, pCrnt, uRowTileSize, uColTileSize, uHaloWidth, uDimWithHalo)
	{
		unsigned int uThreadId  = omp_get_thread_num();

		unsigned int uRowTileId	= uThreadId / uColPartitions;
		unsigned int uColTileId	= uThreadId % uColPartitions;

		unsigned int uStartLine	= uRowTileId * uRowTileSize + uHaloWidth;
		unsigned int uStartCol	= uColTileId * uColTileSize + uHaloWidth;

		unsigned int uEndLine	= uStartLine + uRowTileSize;
		uEndLine = (uEndLine > (uDimWithHalo - uHaloWidth)) ? uDimWithHalo - uHaloWidth : uEndLine;

		unsigned int uEndCol	= uStartCol  + uColTileSize;
		uEndCol  = (uEndCol  > (uDimWithHalo - uHaloWidth)) ? uDimWithHalo - uHaloWidth : uEndCol;

		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		T	cardinal0	= 0.0;
		T	diagonal0	= 0.0;
		T	center0		= 0.0;

		unsigned int cntIterations, i, j;

		for (cntIterations = 0; cntIterations < nIters; cntIterations ++)
		{
			// Do Stencil Operation
			for (i = uStartLine; i < uEndLine; i++)
			{
				T * pCenter 	= &pTmp [ i * uDimWithHalo];
				T * pTop 		= pCenter - uDimWithHalo;
				T * pBottom 	= pCenter + uDimWithHalo;
				T * pOut 		= &pCrnt[ i * uDimWithHalo];

				__assume_aligned(pCenter, 	64);
				__assume_aligned(pTop, 		64);
				__assume_aligned(pBottom, 	64);
				__assume_aligned(pOut, 		64);

				#pragma simd vectorlengthfor(float)
				for (j = uStartCol; j < uEndCol; j++)
				{
					cardinal0	= pCenter[j - 1] + pCenter[j + 1] + pTop[j] + pBottom[j];
					diagonal0	= pTop[j - 1] + pTop[j + 1] + pBottom[j - 1] + pBottom[j + 1];
					center0		= pCenter[j];

					pOut[j]		= wcardinal * cardinal0 + wdiag * diagonal0 + wcenter * center0;
				}
			}

			#pragma omp barrier
			;

			// Switch pointers
			T* pAux	= pTmp;
			pTmp 	= pCrnt;
			pCrnt	= pAux;
		} // End For

	} // End Parallel

//printf("%d\n",munlock(pTmp,uImgElements * sizeof(T)));
//printf("%d\n",munlock(pCrnt,uImgElements * sizeof(T)));

	_mm_free(pCrnt);
}

template <class T> void
MICStencil<T>::operator()( Matrix2D<T>& mtx, unsigned int nIters )
{
	unsigned int uDimWithHalo	= mtx.GetNumRows();
	unsigned int uHaloWidth		= LINESIZE / sizeof(T);
	unsigned int uImgElements	= uDimWithHalo * uDimWithHalo;

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	T* p = mtx.GetFlatData();
    __declspec(target(mic),	align(LINESIZE)) T* pIn = mtx.GetFlatData();

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    __declspec(target(mic), align(sizeof(T)))	T wcenter	= this->wCenter;
    __declspec(target(mic), align(sizeof(T)))	T wdiag		= this->wDiagonal;
    __declspec(target(mic), align(sizeof(T)))	T wcardinal	= this->wCardinal;

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	// double tTotal = curr_second();

	#pragma offload target(mic) in(pIn:length(uImgElements) ALLOC RETAIN)
	{
		// Just copy pIn to compute the copy transfer time
	}

	// double tKernel = curr_second();

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	#pragma offload target(mic) in(pIn:length(uImgElements) REUSE RETAIN)	\
								in(uImgElements) in(uDimWithHalo)			\
								in(wcenter) in(wdiag) in(wcardinal)
	{
		//printf("Off: %f %f %f %d %d %d %d\n", wcenter, wdiag, wcardinal, nIters, uDimWithHalo, uImgElements, uHaloWidth);
		stencilOP( pIn, 
				nIters,
				uDimWithHalo, 
				uImgElements,
				uHaloWidth);
	}//*/
			
	/*#pragma offload target(mic) in(pIn:length(uImgElements) REUSE RETAIN)	\
								in(uImgElements) in(uDimWithHalo)			\
								in(wcenter) in(wdiag) in(wcardinal)
    {
		unsigned int uRowPartitions = sysconf(_SC_NPROCESSORS_ONLN) / 4 - 1;
		unsigned int uColPartitions = 4;	// Threads per core for KNC

		unsigned int uRowTileSize	= (uDimWithHalo - 2 * uHaloWidth) / uRowPartitions;
		unsigned int uColTileSize	= (uDimWithHalo - 2 * uHaloWidth) / uColPartitions;

		uRowTileSize = ((uDimWithHalo - 2 * uHaloWidth) % uRowPartitions > 0) ? (uRowTileSize + 1) : (uRowTileSize);

		// printf("Dividing rows int %d partitions and columns into %d partitions.\n", uRowPartitions, uColPartitions);
		// printf("Tile size rows %d, Tile size columns %d.\n", uRowTileSize, uColTileSize);
		// printf("Halo width %d, dimension with halo %d.\n", uHaloWidth, uDimWithHalo);
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		// Should use the "Halo Val" when filling the memory space

        T *pTmp	 = (T*)pIn;
        T *pCrnt = (T*)memset((T*)_mm_malloc(uImgElements * sizeof(T), LINESIZE), 0, uImgElements * sizeof(T));
		
        #pragma omp parallel firstprivate(pTmp, pCrnt, uRowTileSize, uColTileSize, uHaloWidth, uDimWithHalo)
        {
			unsigned int uThreadId  = omp_get_thread_num();

			unsigned int uRowTileId	= uThreadId / uColPartitions;
			unsigned int uColTileId	= uThreadId % uColPartitions;

			unsigned int uStartLine	= uRowTileId * uRowTileSize + uHaloWidth;
			unsigned int uStartCol	= uColTileId * uColTileSize + uHaloWidth;

			unsigned int uEndLine	= uStartLine + uRowTileSize;
			uEndLine = (uEndLine > (uDimWithHalo - uHaloWidth)) ? uDimWithHalo - uHaloWidth : uEndLine;

			unsigned int uEndCol	= uStartCol  + uColTileSize;
			uEndCol  = (uEndCol  > (uDimWithHalo - uHaloWidth)) ? uDimWithHalo - uHaloWidth : uEndCol;

        	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

			T	cardinal0	= 0.0;
			T	diagonal0	= 0.0;
			T	center0		= 0.0;

			unsigned int cntIterations, i, j;

        	for (cntIterations = 0; cntIterations < nIters; cntIterations ++)
        	{
				// Do Stencil Operation
				for (i = uStartLine; i < uEndLine; i++)
				{
					T * pCenter 	= &pTmp [ i * uDimWithHalo];
					T * pTop 		= pCenter - uDimWithHalo;
					T * pBottom 	= pCenter + uDimWithHalo;
					T * pOut 		= &pCrnt[ i * uDimWithHalo];

					__assume_aligned(pCenter, 	64);
					__assume_aligned(pTop, 		64);
					__assume_aligned(pBottom, 	64);
					__assume_aligned(pOut, 		64);

					#pragma simd vectorlengthfor(float)
					for (j = uStartCol; j < uEndCol; j++)
					{
						cardinal0	= pCenter[j - 1] + pCenter[j + 1] + pTop[j] + pBottom[j];
						diagonal0	= pTop[j - 1] + pTop[j + 1] + pBottom[j - 1] + pBottom[j + 1];
						center0		= pCenter[j];

						pOut[j]		= wcardinal * cardinal0 + wdiag * diagonal0 + wcenter * center0;
					}
				}

				#pragma omp barrier
				;

				// Switch pointers
				T* pAux	= pTmp;
				pTmp 	= pCrnt;
				pCrnt	= pAux;
        	} // End For

        } // End Parallel

        _mm_free(pCrnt);
    } // End Offload */

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 
    // tKernel = curr_second() - tKernel;

	#pragma offload target(mic) out(pIn:length(uImgElements) REUSE FREE)
	{
		// Just copy back pIn
	}
	//printf("Final %p %f\n", pIn, pIn[128200]);

	/*
	tTotal = curr_second() - tTotal;

	double N		= (double) (uDimWithHalo - uHaloWidth);
	double fFLOP 	= 1000.0 * N * N * 11.0;
	double kFLOPS	= (fFLOP / (tKernel * 1e9));
	double tFLOPS	= (fFLOP / (tTotal * 1e9));

	printf("\n Kernel GFLOPS (no   PCIe) = %.3f", kFLOPS);
	printf("\n Total  GFLOPS (with PCIe) = %.3f\n", tFLOPS);
	*/
}

void
EnsureStencilInstantiation( void )
{
    MICStencil<float> csf( 0, 0, 0, 0 );
    Matrix2D<float> mf( 2, 2 );
    csf( mf, 0);

    MICStencil<double> csd( 0, 0, 0, 0 );
    Matrix2D<double> md( 2, 2 );
    csd( md, 0);
}
