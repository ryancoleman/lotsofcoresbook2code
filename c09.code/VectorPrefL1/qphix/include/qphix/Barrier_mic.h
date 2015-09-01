/*
Copyright (c) 2014, Intel Corporation

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

   * Redistributions of source code must retain the above copyright notice,
     this list of conditions and the following disclaimer.
   * Redistributions in binary form must reproduce the above copyright
     notice, this list of conditions and the following disclaimer in the
     documentation and/or other materials provided with the distribution.
   * Neither the name of Intel Corporation nor the names of its contributors
     may be used to endorse or promote products derived from this software
     without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

/* file: Barrier_mic.h
 * Purpose to provide a fast barrier on Xeon Phi Knights Corner Co-processor
 * Supplied by Intel Parallel Computing Labs, Bangalore
 */

/* NB: This is a derived work from the original source code. 
 * The some of the constants have been changed from #define-s to 
 * user settable values in the constructor, and the memory management
 * has been made dynamic, in the constructor and destructor.
 * -- Balint Joo, Jefferson Lab */

#ifndef QPHIX_BARRIER_MIC_H
#define QPHIX_BARRIER_MIC_H

#include <immintrin.h>
#include <cmath>
#include <cstdio>
#include <cstdlib>


using namespace std;

namespace QPhiX { 

#define MAX_THREADS_PER_CORE 4

class Barrier {
	typedef struct {
		volatile int counter;
	} atomic_t;

#define atomic_set(v,i) (((v)->counter) = (i))

	static inline int atomic_dec_and_test( atomic_t *v )
	{
		return !(__sync_sub_and_fetch(&v->counter, 1));
	}


	typedef struct {
		int volatile sense;
	} htbar_t;


	typedef struct {
		int volatile htbar_sense;
		htbar_t htbar[MAX_THREADS_PER_CORE];
		int volatile *myflags[2];
		int **partnerflags[2];
		int parity;
		int sense;
	} kingbar_t;

	int num_kingbar_threads, num_kingkar_nodes, lg_num_kingbar_nodes, num_threads_per_core;
	atomic_t kingbar_threads_not_alloced;
	volatile int kingbar_alloc_done;
	kingbar_t **bar;
	int n_threads_per_core;

 public:
 Barrier(int num_threads, int threads_per_core) : n_threads_per_core(threads_per_core)
	{
		int num_nodes = num_threads / n_threads_per_core;
		bar = (kingbar_t **)_mm_malloc(num_nodes * sizeof(kingbar_t *), 64);
		for(int i=0; i < num_nodes; i++) bar[i] = 0x0;

		num_kingbar_threads = num_threads;
		lg_num_kingbar_nodes = (int)ceil(log2(num_nodes));
		atomic_set(&kingbar_threads_not_alloced, num_kingbar_threads);
		kingbar_alloc_done = 0;
	}

	~Barrier() {
	  
	  // bar is an array of length num_nodes
	  // bar[nid] needs to be cleaned up
	  int num_nodes = num_kingbar_threads / n_threads_per_core;
	  for(int nid=0; nid < num_nodes; nid++) { 
	    kingbar_t *node = bar[nid];
	    if(lg_num_kingbar_nodes > 0) {
	      for (int i = 0;  i < 2;  i++) {
		_mm_free((void *)(node->myflags[i]));
		_mm_free((void *)(node->partnerflags[i]));
	      }
	    }
	    _mm_free((void *)node);
	  }
	  
	  _mm_free((void *)bar);
	}


	void init(int tid)
	{
		kingbar_t *node;
		int j, m, n, nid = tid / n_threads_per_core, num_nodes = num_kingbar_threads / n_threads_per_core;

    // only the first of the 4 hardware threads per core sets up a node
		if ((tid % n_threads_per_core) == 0) {
			node = (kingbar_t *)_mm_malloc(sizeof(kingbar_t), 64);

			node->htbar_sense = 1;
			for(int i = 0; i < n_threads_per_core; i++) node->htbar[i].sense = 1;
			if(lg_num_kingbar_nodes > 0) {
				for (int i = 0;  i < 2;  i++) {
					node->myflags[i] = (int *)_mm_malloc(lg_num_kingbar_nodes * sizeof(int), 64);
					node->partnerflags[i] = (int **)_mm_malloc(lg_num_kingbar_nodes * sizeof(int *), 64);
				}
			}
			node->parity = 0;
			node->sense = 1;

			bar[nid] = node;
		}

    // barrier to let all the nodes get installed into bar
		if (atomic_dec_and_test(&kingbar_threads_not_alloced)) {
			atomic_set(&kingbar_threads_not_alloced, num_kingbar_threads);
			kingbar_alloc_done = 1;
		} else {
			while(!kingbar_alloc_done);
		}

    // only the first of the 4 hardware threads per core completes setup
		if ((tid % n_threads_per_core) == 0) {
			for (int i = 0;  i < lg_num_kingbar_nodes;  i++) {
				node->myflags[0][i] = node->myflags[1][i] = 0;
				j = (nid + (1 << i)) % num_nodes;
				node->partnerflags[0][i] = (int *)&bar[j]->myflags[0][i];
				node->partnerflags[1][i] = (int *)&bar[j]->myflags[1][i];
			}
		}

    // barrier to let setup finish
		if (atomic_dec_and_test(&kingbar_threads_not_alloced)) {
			atomic_set(&kingbar_threads_not_alloced, num_kingbar_threads);
			kingbar_alloc_done = 0;
		} else {
			while(kingbar_alloc_done);
		}
		//printf("%s:%d: reached %d\n", __func__, __LINE__, tid);
	}



	void wait(int tid)
	{
		int i, nid = tid / n_threads_per_core, tofs = tid % n_threads_per_core, num_nodes = num_kingbar_threads / n_threads_per_core;
		kingbar_t *node = bar[nid];

    // first sync using a kbar inside the core
		node->htbar[tofs].sense = !node->htbar[tofs].sense;
		if (tofs == 0) {
			for (i = 1;  i < n_threads_per_core;  i++) {
				while (node->htbar[i].sense == node->htbar_sense) _mm_delay_64(100);
			}
		}
		else {
			while (node->htbar_sense != node->htbar[tofs].sense)_mm_delay_64(100);
		}

    // now, the first of the 4 hardware threads per core syncs with the others
		if (tofs == 0) {
			for (i = 0;  i < lg_num_kingbar_nodes;  i++) {
            //printf("thread %d setting %p, waiting on %p (sense %d)\n", tid, node->partnerflags[node->parity][i], &node->myflags[node->parity][i], node->sense);
				*node->partnerflags[node->parity][i] = node->sense;
				while (node->myflags[node->parity][i] != node->sense)_mm_delay_64(100);
			}
			if (node->parity == 1)
				node->sense = !node->sense;
			node->parity = 1 - node->parity;

        // wake up the waiting threads in this core
			node->htbar_sense = node->htbar[tofs].sense;
		}
	}
};

#undef MAX_THREADS_PER_CORE
}; // End Namespace

#endif
