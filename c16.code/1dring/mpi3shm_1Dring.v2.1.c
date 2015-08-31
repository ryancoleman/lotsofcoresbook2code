#include "mpi.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

/*
 * mpi3shm_1Dring.c can be used for the basic functionality testing of MPI-3 shared memory
 * in multi-node enviroment (such as Xeon and Xeon Phi based clusters).
 *
 * Each rank exhanges hello world info (rank, total number of ranks and node name) 
 * with its 2 neighbours (partners) in 1D ring topology under periodic boundary conditions.
 *
 * mpi3shm_1Dring.c serves as a prototype for MPI-3 shm addition to MPPtest halo code. 
 *
 * functions:
 *   get_n_partners   -- gets number of intra- and inter- node partners
 *   print_n_partners -- prints number of intra- and inter- node partners
 *   translate_ranks  -- defines global rank  -> shm commmunicator rank mapping
 *   get_partners_ptrs -- returns pointers to mem windows
 *   main
 *
 *   Example of compilation and usage:
 *   mpiicc -o mpi3shm_1Dring mpi3shm_1Dring.c
 *   mpiicc -mmic -o mpi3shm_1Dring.mic mpi3shm_1Dring.c
 *   export I_MPI_MIC=enable
 *   export I_MPI_MIC_POSTFIX=.mic
 *   mpirun -l -bootstrap ssh -n 112 -machinefile hostfile  ./mpi3shm_1Dring
 *   where hostfile:
 *   esg065:24
 *   esg065-mic0:32
 *   esg066:24
 *   esg066-mic0:32
*/

  int const n_partners=2; /* size of partners array in 1D-ring topology*/
  int verbose = 1; /* Switch on/off debugging printfs; can be controlled by env 1DRING_VERBOSE */


  /* count number of intra and inter node partners */
  void get_n_partners (int rank, int partners[], int partners_map[],
                       int *n_node_partners, int *n_inter_partners)
  { 
   int j, partner;

    for (j=0; j<n_partners; j++)
    {
       partner = partners[j]; /* partner is in the world notation */ 
       if (partner != MPI_PROC_NULL) {
       /* If partner has a valid mapping in shm communicator then it is on the same node */
          partners_map[j] == MPI_UNDEFINED ? (*n_inter_partners)++ : (*n_node_partners)++; 
       }
    }
  }



  /* print number of intra and inter node partners */
  void print_n_partners (int rank, int partners[], int partners_map[],
                         int n_node_partners, int n_inter_partners)
  {
    int j, partner;
    char tmp_str_intra[n_partners*16]; 
    char tmp_str_inter[n_partners*10];
    int pos_intra = 0, pos_inter = 0;

      for (j=0; j<n_partners; j++)
      {
        partner = partners[j]; /* partner is in the world notation */
        if (partner != MPI_PROC_NULL)
        {
           if (partners_map[j] != MPI_UNDEFINED) /* partner j is on the same node  */
                pos_intra += sprintf (&tmp_str_intra[pos_intra], ", %d (%d)", partner, partners_map[j]);
           else
                pos_inter += sprintf (&tmp_str_inter[pos_inter], ", %d", partner);
        }
      }

      if (n_inter_partners)
           printf ("i'm rank %d with %d internode partner%c%s \n", 
                  rank, n_inter_partners, n_inter_partners >1?'s':' ', tmp_str_inter);
       
      if (n_node_partners) 
           printf ("i'm rank %d with %d intranode partner%c%s\n", 
                  rank, n_node_partners, n_node_partners > 1?'s':' ', tmp_str_intra);
  }



  /* defines global rank  -> shmcomm rank mapping;
     output: partners_map is array of ranks in shmcomm  */
  void translate_ranks(MPI_Comm shmcomm, int partners[], int partners_map[])
  {
    MPI_Group world_group, shared_group;
    
    /* create MPI groups for global communicator and shm communicator */
    MPI_Comm_group (MPI_COMM_WORLD, &world_group); 
    MPI_Comm_group (shmcomm, &shared_group);

    MPI_Group_translate_ranks (world_group, n_partners, partners, shared_group, partners_map); 
  }



  /* returns pointers to mem windows partners_ptrs */
  void get_partners_ptrs(MPI_Win win, int partners[], int partners_map[], int **partners_ptrs )
  {
   int j, partner, dsp_unit;
   MPI_Aint sz;

    for (j=0; j<n_partners; j++) 
    {
      partners_ptrs[j] = NULL;
      partner = partners[j];
      if((partner != MPI_PROC_NULL) && (partners_map[j] != MPI_UNDEFINED))
       /* MPI_Win_shared_query queries the process-local address for memory segments created with MPI_Win_allocate_shared.
         This function can return different process-local addresses for the same physical memory on different processes.  */
          MPI_Win_shared_query (win, partners_map[j], &sz, &dsp_unit, &partners_ptrs[j]); /* returns partners_ptrs */
    }
  }



  /*        MAIN    */
  int
  main (int argc, char *argv[])
  {
    /* to be used for hello world exchanges */
    int rank, numtasks, namelen;
    char name[MPI_MAX_PROCESSOR_NAME];

    /* related to MPI-3 shm*/
    MPI_Comm shmcomm; /* shm communicator  */ 
    MPI_Win win;      /* shm window object */ 
    int *mem;         /* shm memory to be allocated on each node */
    int i0, i1; int* i2; /* for reading back from shm */

    /* current rank exchanges hello world info with partners */
    int partners[n_partners];
    int *partners_map;   /* mapping in shm communicator */
    int **partners_ptrs; /* ptrs to shared mem window for each partner*/
    int j, partner, alloc_len;
    int n_node_partners=0, n_inter_partners=0;

     /* non-blocking inter-node */
    MPI_Request *reqs, *rq;
    int rbuf[n_partners]; /* recv buffer */
    int req_num = 2;      /* each inter-node echange needs a pair of MPI_Irecv and MPI_Isend */

    if (getenv("1DRING_VERBOSE")) verbose = 1; /* Switch on/off printfs thru env */

    MPI_Init (&argc, &argv); 
    MPI_Comm_size (MPI_COMM_WORLD, &numtasks);
    MPI_Comm_rank (MPI_COMM_WORLD, &rank);
    MPI_Get_processor_name (name, &namelen);
    /* if (verbose) printf ("Hello world from COMM_WORLD: rank %d of %d is running on %s\n", rank, numtasks, name); */
    

    /* The 1D ring is defined in partners array. It can be easily expanded to the higher order stencils.
       The current rank has 2 neighbours: previous and next, i.e., prev-rank-next */
    partners[0] = rank-1; /* prev */
    partners[1] = rank+1; /* next */
    /* We will use periodic boundary conditions here */
    if (rank == 0)  partners[0] = numtasks - 1;
    if (rank == (numtasks - 1))  partners[1] = 0;

    /* MPI-3 SHM collective creates shm communicator  */
    MPI_Comm_split_type (MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &shmcomm); 

   /*  mapping: global rank  -> shmcomm rank is in partners_map */
    partners_map = (int*)malloc(n_partners*sizeof(int)); /* allocate partners_map */
    translate_ranks(shmcomm, partners, partners_map);

   /* number of inter and intra node partners */
    get_n_partners (rank, partners, partners_map,  &n_node_partners, &n_inter_partners); 
    if (verbose) print_n_partners (rank, partners, partners_map,  n_node_partners, n_inter_partners); 

    
    alloc_len = 2*sizeof(int) + namelen+1; /* the size of hello world info: 2 int and string; +1 for '\n' */
    if (n_node_partners > 0)
    {
     /* allocate shared memory windows on each node for intra-node partners  */
     MPI_Win_allocate_shared (alloc_len, 1, MPI_INFO_NULL, shmcomm, /* inputs to MPI-3 SHM collective */
                              &mem, &win);  /* outputs: mem - initial address of window; win - window object */

     /* pointers to mem windows */
     partners_ptrs = (int **)malloc(n_partners*sizeof(int*));  
     get_partners_ptrs (win, partners, partners_map,  partners_ptrs );

    }
    else
    {
       mem = (int *)malloc(alloc_len); 
    }    

    /* allocate MPI Request resources for inter-node comm. */
    if(n_inter_partners > 0)
    {
        reqs = (MPI_Request*)malloc(req_num*n_inter_partners*sizeof(MPI_Request));
        rq = reqs;
    }

    /* start halo exchange */

    if (n_node_partners > 0)
    {    
        /* Entering MPI-3 RMA access epoch required for MPI-3 shm */
        MPI_Win_lock_all (MPI_MODE_NOCHECK, win); 
        /*  alternatively, MPI_Win_lock_all, MPI_Win_sync and MPI_Barrier can be replaced with
        2 MPI_Win_fence calls surrounding update of shared memory.  */
        /* MPI_Win_fence(0, win); */ /* -- alternative */
    }    

    /* update MPI-3 shared memory (or local memory in case of lack of node partners) 
     * by writing hello_world info into mem */
    mem[0] = rank; 
    mem[1] = numtasks;
    memcpy(mem+2, name, namelen);

    if (n_node_partners > 0)
    {    
        /* MPI_Win_fence (0, win); */ /* -- alternative end */

        MPI_Win_sync (win);     /* memory fence to sync node exchanges */ 
        MPI_Barrier (shmcomm);  /* time barrier to make sure all ranks have updated their info */
    }


    for (j=0; j<n_partners; j++) 
    {
        if(partners_map[j] != MPI_UNDEFINED) /* partner j is on the same node  */ 
        {
            i0 = partners_ptrs[j][0]; /* load from MPI-3/SHM ops! */
            i1 = partners_ptrs[j][1];
            i2 = partners_ptrs[j]+2; 
            if(verbose) printf ("load MPI/SHM values from neighbour => rank %d, numtasks %d on %s\n", i0, i1, i2);
        }
        else /* inter-node non-blocking MPI-1 */
        {
            MPI_Irecv (&rbuf[j], 1, MPI_INT, partners[j], 1 , MPI_COMM_WORLD, rq++);
            MPI_Isend (&rank, 1, MPI_INT, partners[j], 1 , MPI_COMM_WORLD, rq++);
        }
    } 

   /* sync inter-node exchanges and print out receive buffer rbuf*/
   if(n_inter_partners > 0)
   {
       MPI_Waitall (req_num*n_inter_partners, reqs, MPI_STATUS_IGNORE); 
       if(verbose){
           for (j =0; j< n_partners;j++)
               if (partners_map[j] == MPI_UNDEFINED) printf("Recieved from my inter-node partner %d\n", rbuf[j]);
       }
   }    
   
   if (n_node_partners > 0)
   {    
       MPI_Win_unlock_all (win);  /* close RMA epoch */
       /* free resources */
       MPI_Win_free (&win);
       free (partners_ptrs); 
   } 
   
   if (n_inter_partners) free (reqs);
   
   free (partners_map);  
   
   MPI_Finalize ();

   return (0);
 }


